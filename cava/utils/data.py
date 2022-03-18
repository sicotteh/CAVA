#!/usr/bin/env python3


# Classes providing interfaces with annotation databases and the reference genome
#######################################################################################################################

import os
import sys

from . import conseq
from . import core
from . import csn

#import time

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + '/pysamdir')
import pysam


#######################################################################################################################

# Class representing the Ensembl (transcript) dataset ( can be any database)
class Ensembl(object):
    # Constructor
    def __init__(self, options, genelist, transcriptlist, codon_usage, reference):
        self.options = options
        # Openning tabix file representing the Ensembl database
        self.contigs = dict()
        try:
            self.tabixfile = pysam.Tabixfile(options.args['ensembl'])
        except:
            try:
                self.tabixfile = pysam.TabixFile(options.args['ensembl'])
            except:
                sys.stderr.write("CAVA: ERROR: error trying to open Tabix file for "+options.args['ensembl']+"\n")
        for chrom in self.tabixfile.contigs:
            if chrom in reference.reflens:
                self.contigs[chrom] = reference.reflens[chrom]
        self.proteinSeqs = dict()
        self.exonSeqs = dict()
        self.exoncache_hit = dict()
        self.cache_num = 0
        self.genelist = genelist
        self.transcriptlist = transcriptlist
        self.codon_usage = codon_usage
        # Transcript to Protein Map for HGVSp protein
        # copy it over to "self" in order to maintain the calling function signature of Record.annotate() called by run() (main.py)
        self.transcript2protein = options.transcript2protein
        self.nvar = 0  # counter for most recent transcript lookup

        # 2E4 genes, 200K isoforms in 3E9, so 11E5 (50K) bins should get 0-1 genes and
        #
        #  Cache fully parsed transcripts
        self.lasttranscript = None
        self.transcript_cache = dict()
        self.CACHESIZE = 10  # This will support multi-transcript queries. If use a lot of alternative splicing transcripts.. need to up that.
        self.transcript_nvar = dict()
        # Cache transcript positions
        self.transcript_bins = None
        self.chrom = None


        self.nbins = 0
        self.binsize = 50000
        self.chrom = None
        if 'loadalltranscripts' in self.options.args and self.options.args['loadalltranscripts'] is True:
            self.loadalltranscripts = True
        else:
            self.loadalltranscripts = False


# Get the list of transcript lines overlapping
# returns either an iterator over a tabix file .. or a list (that can be iterated over)
    def fetch_overlapping_transcripts(self,chrom,startpos0,endpos1): # Give 0-base coordinate for start and 1-base for stop
        # If current chromosome is not loaded, then
        #      get tabix iterator figure out length, figu.. and load all chromosomes lines
        # Check self.tabixfile.l
        if chrom not in self.contigs:
            return list()
        if self.loadalltranscripts is False:
            # This fetch consumes most of the runtime for CAVA (78%) as long
            # caching should be faster.
            return self.tabixfile.fetch(reference=chrom, start = startpos0, end = endpos1)
        if self.chrom is None or chrom != self.chrom:
            # Flush cache and load all transcripts.
            self.chrom = chrom
            self.transcript_bins = [None]*(1+int((1+self.contigs[chrom])/self.binsize))
            hits = self.tabixfile.fetch(reference = chrom)
            for line in hits:
                linedat = line.split("\t",8)
                transcriptid = linedat[0]
                transcriptStart = int(linedat[6])
                transcriptEnd = int(linedat[7])
                binstart = int(transcriptStart/self.binsize)
                binend  = int(transcriptEnd/self.binsize)
                for bin in range(binstart,binend+1):
                    if self.transcript_bins[bin] is None:
                        self.transcript_bins[bin] = [[transcriptid,transcriptStart,transcriptEnd,line]]
                    else:
                        self.transcript_bins[bin].append([transcriptid,transcriptStart,transcriptEnd,line])
        lines = list()
        got_transcript = dict()
        binstart = int((startpos0+1)/ self.binsize)
        binend = int(endpos1/ self.binsize)
        for bin in range(binstart, binend + 1):
            bin_list = self.transcript_bins[bin]
            if bin_list is not None:
                for bin_content in bin_list: # [transcriptid, transcriptStart, transcriptEnd, line])
                    if (startpos0+1) >= bin_content[1] and endpos1<=bin_content[2]:
                        transcriptid = bin_content[0]
                        if transcriptid not in got_transcript:
                            got_transcript[transcriptid]=1
                            lines.append(bin_content[3])
        return lines





#
# Loading transcripts and the exons is very costly (need to read all exons), so caching them save a lot of disk access
#

    def find_transcript_in_cache_or_in_file(self,line):
        transcriptid = line.split("\t")[0]
        self.nvar +=1
        self.transcript_nvar[transcriptid] = self.nvar
        if transcriptid in self.transcript_cache:
            transcript = self.transcript_cache[transcriptid]
        else:
            # with the reference sequence being cached, fetching a whole transcripts and exons takes 0.03-0.09 ms
            #        ... rather than 150-200 ms if the sequence was not cached.
            transcript = core.Transcript(line)
            self.transcript_cache[transcriptid] = transcript
            if len(self.transcript_cache)>self.CACHESIZE:
                vals = list(self.transcript_nvar.values())
                minval = min(vals)
                which_minval = vals.index(minval)
                rm_tr = ''+list(self.transcript_nvar.keys())[which_minval]
                self.transcript_cache.pop(rm_tr)
                self.transcript_nvar.pop(rm_tr)

        return transcript


    # Find transcripts overlapping with a variant

    def findTranscripts(self, variant):
        ret = dict()
        retOUT = dict()

        # Checking chromosome name
        goodchrom = core.convert_chrom(variant.chrom,self.tabixfile.contigs)
        if goodchrom is None:
            return ret, retOUT

        # Defining variant end points.
        # HS: Tabix uses 0-based indexing for start/pos
        if not variant.is_insertion:
            start = variant.pos -1
            end = variant.pos + len(variant.ref)
        else:  # for insertion, position got shifted to be after the insertion point .. shift back
            start = variant.pos - 2
            end = variant.pos -1
        if start<0:
            start = 0
        if end <= start:
            end = start + 1

        # Checking both end points of the variant
        reg1 = goodchrom + ':' + str(start) + '-' + str(start)
        reg2 = goodchrom + ':' + str(end) + '-' + str(end)

        if not variant.is_substitution:
            # HS notes:
            # tabix index is loaded in memory, and the data is buffered in 37K (2^16) blocks (from reading the pysam and htslib)
            #   the 'fetch' take about 0.004 - 0.049 ms .. MUCH faster than a single disk seek & read (10-15ms)
            # the iteration take 0.8 ms for the two .. must be cached .. but slow.
            # but uses 70% of the time of CAVA.
            # Loading transcript tables in memory and a pyranges implmementaion would help.
            #  ... so so it is clearly cached .. and no need for further optimization.
            #st_time = time.perf_counter_ns()
            #    .. so basically a sliding windown caching
            # Using pyranges takes 8 ms per pyrange creation.. and 33 ms for an overlap operation with all transcripts.
            # ... but if refactored to find the overlap of 1000's of SNPs/intervals at a time
            #  .. would be faster... since finding overlap with 1000 snps only takes 82ms, so 1.6ms for 2 overlap
            #
            # so pyranges would not be faster .. even if batches.
            #
            # The only way to get faster is to build an im-,memory genome-binned
            # Create per chromosomes "bins" (100K bins,so about 30,000 elements dictionary .. poointing to  a list of transcripts.
            # .. cache last bin request.
            hits1 = self.fetch_overlapping_transcripts(goodchrom,start,start+1)  #self.tabixfile.fetch(region=reg1)
            hits2 = self.fetch_overlapping_transcripts(goodchrom,end-1,end)  # self.tabixfile.fetch(region=reg2)
            #end_time0 = time.perf_counter_ns()
            #sys.stdout.write("Time for two tabix transcripts fetch=" + str(end_time - st_time) + "\n")

            hitdict1 = dict()
            hitdict2 = dict()
            for line in hits1:
                transcript = self.find_transcript_in_cache_or_in_file(line)
                if not (transcript.transcriptStart + 1 <= start <= transcript.transcriptEnd): continue
                # if not strand == transcript.strand: continue
                hitdict1[transcript.TRANSCRIPT] = transcript
            for line in hits2:
                transcript = self.find_transcript_in_cache_or_in_file(line)
                if not (transcript.transcriptStart + 1 <= end <= transcript.transcriptEnd): continue
                #  if not strand == transcript.strand: continue
                hitdict2[transcript.TRANSCRIPT] = transcript
            #end_time = time.perf_counter_ns()
            #sys.stdout.write("Two tabix fetch & iterate=" + str(end_time0 - st_time) + " "+ str(end_time - st_time) + "\n")

            for key, transcript in hitdict1.items():
                if len(self.genelist) > 0 and transcript.geneSymbol not in self.genelist: continue
                if len(self.transcriptlist) > 0 and transcript.TRANSCRIPT not in self.transcriptlist: continue

                if key in list(hitdict2.keys()):
                    ret[key] = transcript
                else:
                    if not variant.is_insertion: retOUT[key] = transcript  # partial overlap downstream of transcrupt

            if not variant.is_insertion:
                for key, transcript in hitdict2.items():  # check for partial overlap upstream of transcript
                    if len(self.genelist) > 0 and transcript.geneSymbol not in self.genelist: continue
                    if len(self.transcriptlist) > 0 and transcript.TRANSCRIPT not in self.transcriptlist: continue

                    if not key in list(hitdict1.keys()):
                        retOUT[key] = transcript

        else:  # Variant is Substitution
            hits1 = self.fetch_overlapping_transcripts(goodchrom,start,end)  #self.tabixfile.fetch(region=reg2)
            for line in hits1:
                transcript = self.find_transcript_in_cache_or_in_file(line)

                if len(self.genelist) > 0 and transcript.geneSymbol not in self.genelist: continue
                if len(self.transcriptlist) > 0 and transcript.TRANSCRIPT not in self.transcriptlist: continue

                if not (transcript.transcriptStart + 1 <= end <= transcript.transcriptEnd): continue
                ret[transcript.TRANSCRIPT] = transcript

        return ret, retOUT  # retOUT not populated for substitution

    # Check if a is between x and y
    def inrange(self, x, y, a):
        return x <= a <= y or y <= a <= x

    # Parse CSN coordinates
    def getIntronBases(self, x):
        idx = x.find('-')
        if idx < 1:
            idx = x.find('+')
            if idx < 1: return None
            return int(x[idx:])
        return int(x[idx:])

    # Check if variant is duplication overlapping SS boundary
    def isDupOverlappingSSBoundary(self, csnval, ssrange=8):

        if '_p' in csnval:
            [cpart, _] = csnval.split('_p')
        else:
            cpart = csnval

        idx = cpart.find('dup')
        if idx == -1: return False

        cpart = cpart[2:idx]

        if '_' in cpart:
            [x, y] = cpart.split('_')
        else:
            x = self.getIntronBases(cpart)
            if x is None: return False
            return x == ssrange or x == -ssrange

        x = self.getIntronBases(x)
        y = self.getIntronBases(y)
        if x is None or y is None: return False

        return self.inrange(x, y, ssrange) or self.inrange(x, y, -ssrange)

    # Correct CLASS annotations for duplications overlapping SS boundary
    def correctClasses(self, csn, class_plus, class_minus):
        if self.isDupOverlappingSSBoundary(csn, ssrange=int(self.options.args['ssrange'])):
            if class_plus == 'SS' and class_minus == 'INT': return 'INT', 'INT'
            if class_plus == 'INT' and class_minus == 'SS': return 'INT', 'INT'
        return class_plus, class_minus

    # Correct SO annotations for duplications overlapping SS boundary
    def correctSOs(self, csn, so_plus, so_minus):
        if self.isDupOverlappingSSBoundary(csn):
            if so_plus == 'intron_variant|splice_region_variant' and so_minus == 'intron_variant': return 'intron_variant', 'intron_variant'
            if so_plus == 'intron_variant' and so_minus == 'intron_variant|splice_region_variant': return 'intron_variant', 'intron_variant'
        return so_plus, so_minus

    # Annotating a variant based on Ensembl data

    def annotate(self, variant, reference, impactdir):
        # Create left-aligned and right-aligned versions of the variant
        if variant.is_deletion or variant.is_insertion:  # optimization, MNP or substitutions cannot be shifted
            variant_plus = variant.alignOnPlusStrand(reference)
            variant_minus = variant.alignOnMinusStrand(reference)
        else:
            variant_plus = variant
            variant_minus = variant

        # Checking if variant alignment makes any difference
        if variant_plus.pos == variant_minus.pos:
            difference = False
        else:
            difference = True

        # Initializing annotation strings
        TRANSCRIPTstring = ''
        GENEstring = ''
        GENEIDstring = ''
        LOCstring = ''
        CSNstring = ''
        CLASSstring = ''
        ALTFLAGstring = ''
        TRINFOstring = ''
        ALTANNstring = ''
        ALTCLASSstring = ''
        SOstring = ''
        ALTSOstring = ''
        IMPACTstring = ''
        PROTPOSstring = ''
        PROTREFstring = ''
        PROTALTstring = ''

        # Collecting transcripts that overlap with the variant
        transcripts_plus, transcriptsOUT_plus = self.findTranscripts(variant_plus)
        transcripts_minus, transcriptsOUT_minus = self.findTranscripts(variant_minus)

        if variant.is_deletion or variant.is_complex:
            # If variant is Deletion(or Del+repl),
            #         being partial at 3' end has predictable functional impact
            #         Deletion 5' end should be annotated because CAP site is not always correct)
            #         we want the CSN produced
            transcripts_allplus = set(list(transcripts_plus.keys()) + list(transcriptsOUT_plus.keys()))
            transcripts_allminus = set(list(transcripts_minus.keys()) + list(transcriptsOUT_minus.keys()))
            transcripts = set(list(transcripts_plus.keys()) + list(transcripts_minus.keys()) + list(transcriptsOUT_plus.keys()) + list(transcriptsOUT_minus.keys()))
            transcriptsOUT = set()
        else:
            transcripts_allplus = set(list(transcripts_plus.keys()))
            transcripts_allminus = set(list(transcripts_minus.keys()))
            transcripts = set(list(transcripts_plus.keys()) + list(transcripts_minus.keys()))
            transcriptsOUT = set(list(transcriptsOUT_plus.keys()) + list(transcriptsOUT_minus.keys()))

        transcripts = sorted(list(transcripts))
        transcriptsOUT = sorted(list(transcriptsOUT))
        transcripts_allplus = sorted(transcripts_allplus)
        transcripts_allminus = sorted(transcripts_allminus)

        # Annotating with transcripts that only partially overlap with the variant
        for TRANSCRIPT in transcriptsOUT:
            # if variant can be shifted (insertion or deletion), then it is possible
            # for a transcript to be both in the transcripts list and the transcriptsOUT lists
            if TRANSCRIPT in transcripts: continue
            # Only 1 end mapped in transcript
            #  either statement is true for both shifted variants pos .. or in transcriptt for one shifted version and not in for the other shifted version

            if TRANSCRIPT in list(transcriptsOUT_plus.keys()):
                transcript = transcriptsOUT_plus[TRANSCRIPT]
            else:
                transcript = transcriptsOUT_minus[TRANSCRIPT]

            if len(TRANSCRIPTstring) > 0:
                TRANSCRIPTstring += ':'
                GENEstring += ':'
                GENEIDstring += ':'
                TRINFOstring += ':'
                LOCstring += ':'
                CSNstring += ':'
                CLASSstring += ':'
                ALTFLAGstring += ':'
                ALTANNstring += ':'
                ALTCLASSstring += ':'
                SOstring += ':'
                ALTSOstring += ':'
                IMPACTstring += ':'
                PROTPOSstring += ':'
                PROTREFstring += ':'
                PROTALTstring += ':'

            TRANSCRIPTstring += TRANSCRIPT
            GENEstring += transcript.geneSymbol
            GENEIDstring += transcript.geneID
            TRINFOstring += transcript.TRINFO

            # if TRANSCRIPT in list(transcriptsOUT_plus.keys()):
            #         LOCstring += 'OUT'
            #     else:
            #         LOCstring += '.'
            # else:
            #     if TRANSCRIPT in list(transcriptsOUT_minus.keys()):
            #         LOCstring += 'OUT'
            #     else:
            #         LOCstring += '.'
            LOCstring += 'OUT'

            CSNstring += '.'
            CLASSstring += '.'
            ALTFLAGstring += '.'
            ALTANNstring += '.'
            ALTCLASSstring += '.'
            SOstring += '.'
            ALTSOstring += '.'
            IMPACTstring += '.'
            PROTPOSstring += '.'
            PROTREFstring += '.'
            PROTALTstring += '.'

        # Iterating through the list of transcripts with  variant fully inside transcript
        #  or Del/MNP with one end in transcript
        for TRANSCRIPT in transcripts:
            if TRANSCRIPT in list(transcripts_plus.keys()):
                transcript = transcripts_plus[TRANSCRIPT]
            elif TRANSCRIPT in list(transcripts_minus.keys()):
                transcript = transcripts_minus[TRANSCRIPT]
            elif TRANSCRIPT in list(transcriptsOUT_plus.keys()):
                transcript = transcriptsOUT_plus[TRANSCRIPT]
            elif TRANSCRIPT in list(transcriptsOUT_minus.keys()):
                transcript = transcriptsOUT_minus[TRANSCRIPT]

            # Separating annotations by different transcripts with colon
            if len(TRANSCRIPTstring) > 0:
                TRANSCRIPTstring += ':'
                GENEstring += ':'
                GENEIDstring += ':'
                TRINFOstring += ':'
                LOCstring += ':'
                CSNstring += ':'
                CLASSstring += ':'
                ALTFLAGstring += ':'
                ALTANNstring += ':'
                ALTCLASSstring += ':'
                SOstring += ':'
                ALTSOstring += ':'
                IMPACTstring += ':'
                PROTPOSstring += ':'
                PROTREFstring += ':'
                PROTALTstring += ':'

            # Creating the TRANSCRIPT, GENE and TRINFO annotations
            TRANSCRIPTstring += TRANSCRIPT
            GENEstring += transcript.geneSymbol
            GENEIDstring += transcript.geneID
            TRINFOstring += transcript.TRINFO

            # Creating the LOC annotation
            if TRANSCRIPT in list(transcripts_plus.keys()):
                loc_plus = transcript.whereIsThisVariant(variant_plus)
            elif TRANSCRIPT in list(transcriptsOUT_plus.keys()):
                loc_plus = 'OUT'
#                loc_plus = transcript.whereIsThisVariant(variant_plus)
            else:
                loc_plus = '.'

            if difference:
                if TRANSCRIPT in list(transcripts_minus.keys()):
                    loc_minus = transcript.whereIsThisVariant(variant_minus)
                elif TRANSCRIPT in list(transcriptsOUT_minus.keys()):
                    loc_minus = 'OUT'
 #                   loc_minus = transcript.whereIsThisVariant(variant_minus)
                else:
                    loc_minus = '.'
            else:
                loc_minus = loc_plus

            # Creating reference and mutated protein sequence
            #   . means outside transcript .. which can (now) happen for deletion or MNP
            notexonic_plus = (
                    ('5UTR' in loc_plus) or ('3UTR' in loc_plus) or ('-' in loc_plus) or ('In' in loc_plus) or (
                    loc_plus == 'OUT') or (loc_plus == '.'))
            if difference:
                notexonic_minus = (
                        ('5UTR' in loc_minus) or ('3UTR' in loc_minus) or ('-' in loc_minus) or (
                        'In' in loc_minus) or (loc_minus == 'OUT') or (loc_minus == '.'))  # Variants overlapping the ends of the protein
            else:
                notexonic_minus = notexonic_plus
            exonseqs = None
            if notexonic_plus and notexonic_minus:
                protein = ''
            else:
                self.cache_num += 1
                self.exoncache_hit[transcript.TRANSCRIPT] = self.cache_num
                if not transcript.TRANSCRIPT in list(self.proteinSeqs.keys()):
                    protein, exonseqs = transcript.getProteinSequence(reference, None, None, self.codon_usage)
                    self.proteinSeqs[transcript.TRANSCRIPT] = protein
                    self.exonSeqs[transcript.TRANSCRIPT] = exonseqs
                    transcript.exonseqs = exonseqs
                    if len(self.proteinSeqs) > self.CACHESIZE:  # Cache of proteins and exons data
                            vals = list(self.exoncache_hit.values())
                            minval = min(vals)
                            which_minval = vals.index(minval)
                            rm_tr = '' + list(self.exoncache_hit.keys())[which_minval]
                            self.exonSeqs.pop(rm_tr)
                            self.exoncache_hit.pop(rm_tr)
                            self.proteinSeqs.pop(rm_tr)
                else:
                    protein = self.proteinSeqs[transcript.TRANSCRIPT]
                    exonseqs = self.exonSeqs[transcript.TRANSCRIPT]


            if notexonic_plus or protein == '':
                mutprotein_plus = ''
            else:
                mutprotein_plus, exonseqsalt_plus = transcript.getProteinSequence(reference, variant_plus, exonseqs, self.codon_usage)

            if difference:
                if notexonic_minus or protein == '':
                    mutprotein_minus = ''
                else:
                    mutprotein_minus, exonseqsalt_minus = transcript.getProteinSequence(reference, variant_minus, exonseqs, self.codon_usage)
            else:
                mutprotein_minus = mutprotein_plus

            # Creating the CSN annotations both for left and right aligned variant
            if TRANSCRIPT in transcripts_allplus:
                csn_plus, protchange_plus = csn.getAnnotation(variant_plus, transcript, reference, protein, mutprotein_plus)
                csn_plus_str = csn_plus.getAsString()
            else:
                csn_plus_str, protchange_plus = '.', ('.', '.', '.')

            if difference:
                if TRANSCRIPT in transcripts_allminus:
                    csn_minus, protchange_minus = csn.getAnnotation(variant_minus, transcript, reference, protein, mutprotein_minus)
                    csn_minus_str = csn_minus.getAsString()
                else:
                    csn_minus_str, protchange_minus = '.', ('.', '.', '.')
            else:
                csn_minus_str, protchange_minus = csn_plus_str, protchange_plus

            # CLASS, SO and IMPACT

            so_plus = ''
            so_minus = ''
            class_plus = ''
            class_minus = ''
            impact_plus = ''
            impact_minus = ''

            if not impactdir is None or self.options.args['ontology'].upper() in ['CLASS', 'BOTH']:
                # Creating the CLASS annotations both for left and right aligned variant
                if TRANSCRIPT in transcripts_allplus:
                    class_plus = conseq.getClassAnnotation(variant_plus, transcript, protein, mutprotein_plus, loc_plus,
                                                           int(self.options.args['ssrange']))
                else:
                    class_plus = '.'

                if difference:
                    if TRANSCRIPT in transcripts_allminus:
                        class_minus = conseq.getClassAnnotation(variant_minus, transcript, protein, mutprotein_minus,
                                                                loc_minus, int(self.options.args['ssrange']))
                    else:
                        class_minus = '.'
                else:
                    class_minus = class_plus

            # Determining the IMPACT flag
            if not impactdir is None:
                if TRANSCRIPT in transcripts_allplus:
                    if class_plus in list(impactdir.keys()):
                        impact_plus = impactdir[class_plus]
                    else:
                        impact_plus = 'None'
                else:
                    impact_plus = '.'

                if TRANSCRIPT in transcripts_allminus:
                    if class_minus in list(impactdir.keys()):
                        impact_minus = impactdir[class_minus]
                    else:
                        impact_minus = 'None'
                else:
                    impact_minus = '.'

            if self.options.args['ontology'].upper() in ['SO', 'BOTH']:
                # Creating the SO annotations both for left and right aligned variant

                if TRANSCRIPT in transcripts_allplus:
                    so_plus = conseq.getSequenceOntologyAnnotation(variant_plus, transcript, protein, mutprotein_plus,
                                                                   loc_plus)
                else:
                    so_plus = '.'

                if difference:
                    if TRANSCRIPT in transcripts_allminus:
                        so_minus = conseq.getSequenceOntologyAnnotation(variant_minus, transcript, protein,
                                                                        mutprotein_minus, loc_minus)
                    else:
                        so_minus = '.'
                else:
                    so_minus = so_plus

            # Deciding which is the correct CSN and CLASS annotation
            if transcript.strand == 1:
                class_plus, class_minus = self.correctClasses(csn_plus_str, class_plus, class_minus)
                so_plus, so_minus = self.correctSOs(csn_plus_str, so_plus, so_minus)
                if class_plus == "ESS" or so_plus.startswith("splice"):
                    csn_plus_arr = csn_plus_str.split("_p.")
                    if len(csn_plus_arr)==1 or csn_plus_arr[1] in ["=","(=)","="]:
                        csn_plus_str = csn_plus_arr[0]  + "_p.?"
                if class_minus == "ESS" or so_minus.startswith("splice"):
                    csn_minus_arr = csn_minus_str.split("_p.")
                    if len(csn_minus_arr) == 1 or csn_minus_arr[1] in ["=", "(=)","="]:
                        csn_minus_str = csn_minus_arr[0] + "_p.?"
                CSNstring += csn_plus_str
                CLASSstring += class_plus
                ALTANN = csn_minus_str
                altCLASS = class_minus
                SOstring += so_plus
                altSO = so_minus
                LOCstring += loc_plus
                IMPACTstring += impact_plus
                PROTPOSstring += protchange_plus[0]
                PROTREFstring += protchange_plus[1]
                PROTALTstring += protchange_plus[2]
            else:
                class_plus, class_minus = self.correctClasses(csn_minus_str, class_plus, class_minus)
                so_plus, so_minus = self.correctSOs(csn_minus_str, so_plus, so_minus)
                if class_plus == "ESS" or so_plus.startswith("splice"):
                    csn_plus_arr = csn_plus_str.split("_p.")
                    if len(csn_plus_arr)==1 or csn_plus_arr[1] in ["=","(=)"]:
                        csn_plus_str = csn_plus_arr[0]  + "_p.?"
                if class_minus == "ESS" or so_minus.startswith("splice"):
                    csn_minus_arr = csn_minus_str.split("_p.")
                    if len(csn_minus_arr) == 1 or csn_minus_arr[1] in ["=", "(=)"]:
                        csn_minus_str = csn_minus_arr[0] + "_p.?"
                CSNstring += csn_minus_str
                CLASSstring += class_minus
                ALTANN = csn_plus_str
                altCLASS = class_plus
                SOstring += so_minus
                altSO = so_plus
                LOCstring += loc_minus
                IMPACTstring += impact_minus
                PROTPOSstring += protchange_minus[0]
                PROTREFstring += protchange_minus[1]
                PROTALTstring += protchange_minus[2]

            if self.options.args['givealt']:
                # Creating the ALTANN annotation
                if not csn_plus_str == csn_minus_str:
                    ALTANNstring += ALTANN
                else:
                    ALTANNstring += '.'

                # Creating the ALTCLASS annotation
                if not class_plus == class_minus:
                    ALTCLASSstring += altCLASS
                else:
                    ALTCLASSstring += '.'

                # Creating the ALTSO annotations
                if not so_plus == so_minus:
                    ALTSOstring += altSO
                else:
                    ALTSOstring += '.'

            if (not self.options.args['givealt']) or self.options.args['givealtflag']:
                # Creating the ALTFLAG annotation

                if self.options.args['ontology'].upper() == 'CLASS':
                    if not class_plus == class_minus:
                        ALTFLAGstring += 'AnnAndClass'
                    else:
                        if not csn_plus_str == csn_minus_str:
                            ALTFLAGstring += 'AnnNotClass'
                        else:
                            ALTFLAGstring += 'None'

                if self.options.args['ontology'].upper() == 'SO':
                    if not so_plus == so_minus:
                        ALTFLAGstring += 'AnnAndSO'
                    else:
                        if not csn_plus_str == csn_minus_str:
                            ALTFLAGstring += 'AnnNotSO'
                        else:
                            ALTFLAGstring += 'None'

                if self.options.args['ontology'].upper() == 'BOTH':
                    if not class_plus == class_minus:
                        if not so_plus == so_minus:
                            ALTFLAGstring += 'AnnAndClassAndSO'
                        else:
                            ALTFLAGstring += 'AnnAndClassNotSO'
                    else:
                        if csn_plus_str == csn_minus_str:
                            ALTFLAGstring += 'None'
                        else:
                            if not so_plus == so_minus:
                                ALTFLAGstring += 'AnnAndSONotClass'
                            else:
                                ALTFLAGstring += 'AnnNotClassNotSO'

        # Adding annotations to the variant
        variant.addFlag('TRANSCRIPT', TRANSCRIPTstring)
        variant.addFlag('GENE', GENEstring)
        variant.addFlag('GENEID', GENEIDstring)
        variant.addFlag('TRINFO', TRINFOstring)
        variant.addFlag('LOC', LOCstring)
        variant.addFlag('CSN', CSNstring)
        variant.addFlag('PROTPOS', PROTPOSstring)
        variant.addFlag('PROTREF', PROTREFstring)
        variant.addFlag('PROTALT', PROTALTstring)

        if self.options.args['ontology'].upper() in ['CLASS', 'BOTH']: variant.addFlag('CLASS', CLASSstring)
        if self.options.args['ontology'].upper() in ['SO', 'BOTH']: variant.addFlag('SO', SOstring)

        if not impactdir is None: variant.addFlag('IMPACT', IMPACTstring)

        if self.options.args['givealt']:
            variant.addFlag('ALTANN', ALTANNstring)
            if self.options.args['ontology'].upper() in ['CLASS', 'BOTH']: variant.addFlag('ALTCLASS', ALTCLASSstring)
            if self.options.args['ontology'].upper() in ['SO', 'BOTH']: variant.addFlag('ALTSO', ALTSOstring)

        if (not self.options.args['givealt']) or self.options.args['givealtflag']:
            variant.addFlag('ALTFLAG', ALTFLAGstring)

        return variant


#######################################################################################################################


# Class representing the dbSNP dataset
class dbSNP(object):
    # Constructor
    def __init__(self, options):
        # Openning tabix file representing the dbSNP database
        self.tabixfile = pysam.Tabixfile(options.args['dbsnp'])

    # Annotating a variant based on dbSNP data
    def annotate(self, variant):
        # Checking if variant is a SNP at all
        if variant.is_substitution:
            # Fetching data from dbSNP database
            goodchrom = core.convert_chrom(variant.chrom,self.tabixfile.contigs)
            if goodchrom is None:
                variant.addFlag('DBSNP', '')
                return variant
            reg = goodchrom + ':' + str(variant.pos) + '-' + str(variant.pos)
            lines = self.tabixfile.fetch(region=reg)
            for line in lines:
                cols = line.split('\t')
                # Adding DBSNP annotation to the variant
                alts = cols[3].split(',')
                if variant.alt in alts:
                    variant.addFlag('DBSNP', cols[0])
            if not 'DBSNP' in variant.flags: variant.addFlag('DBSNP', '')
        else:
            variant.addFlag('DBSNP', '')
        return variant


#######################################################################################################################

# Class representing the reference genome dataset
class Reference(object):
    # Constructor
    def __init__(self, options):
        # Openning tabix file representing the reference genome
        try:
            self.fastafile = pysam.FastaFile(options.args['reference'])
        except:
            try:  # old API, version prior to 0.8.1
                self.fastafile = pysam.Fastafile(options.args['reference'])
            except:
                sys.stderr.write("CAVA:ERROR, Error, pysam API invalid\n")

        lengths = self.fastafile.lengths
        references = self.fastafile.references
        self.reflens = dict()
        for i in range(0,len(lengths)):
            chrom = str(references[i])
            self.reflens[chrom] = lengths[i]
            self.reflens[chrom.upper()] = lengths[i]
            self.reflens[chrom.lower()] = lengths[i]
            if chrom.startswith("chr"):
                self.reflens[chrom[3:]] = lengths[i]
                self.reflens[chrom[3:].upper()] = lengths[i]
                self.reflens[chrom[3:].lower()] = lengths[i]
            else:
                self.reflens["chr"+chrom] = lengths[i]
                self.reflens["CHR"+chrom.upper()] = lengths[i]
                self.reflens["chr" + chrom.upper()] = lengths[i]
                self.reflens["chr"+chrom.lower()] = lengths[i]

            if chrom in ["chrM","chrMT","M","MT","chrM"]:
                self.reflens["chrM"] = lengths[i]
                self.reflens["chrMT"] = lengths[i]
                self.reflens["CHRM"] = lengths[i]
                self.reflens["CHRMT"] = lengths[i]
                self.reflens["M"] = lengths[i]
                self.reflens["MT"] = lengths[i]
                self.reflens["chrm"] = lengths[i]
                self.reflens["chrmt"] = lengths[i]
                self.reflens["m"] = lengths[i]
                self.reflens["mt"] = lengths[i]
        self.cache = ""
        self.CACHESIZE = 5000000
        self.PAD = 1000000  # Left Padding for normalization of large indels(5000bp) or  exons for transcripts from - strand
        self.chrom = ""
        self.start0 = 0
        self.endpos = 0

    # Private method .. do not call because assume that chrom is in sequence index
    #  also assumes end is less than chrom length .. and that start >=1

    def __getseq_from_cache_or_file(self, chrom, start, endpos):
        # XXX-HS To make faster, could retrieve large blocks of sequence
        # For normalization.
        # There are about 322,000 indels per 3XE9 /person, so 0.3/1000bp
        # so for a cache to be useful, it needs to include variants, so 30KB min.
        #  for exomes. Most genes are < 1MB .. and 2000 bp cds.. so would only include 1-2 indels/sample
        # Assuming the tabix index is in RAM, reading from disk is seek time (5-10ms) + read time (1MB/ms)
        #    so reading 10-20MB takes the same amount of time as a single disk seek.
        # For exome, there is one indel/gene and each gene is 1MB (average), so a 5MB cache insures 5X cache reuse(hits)
        # if dealing with SNPs or whole genome data, other elements, a smaller cache will reap benefits too .. up to about 5MB.
        if chrom != self.chrom or start < self.start0+1 or endpos > self.endpos:
            fetch_start = start - self.PAD
            if fetch_start<1:
                fetch_start =1
            fetch_end = start + self.CACHESIZE
            if fetch_end < endpos:
                fetch_end = endpos + self.PAD
            if fetch_end > self.reflens[chrom]:
                fetch_end = self.reflens[chrom]
            self.start0 = fetch_start-1
            self.endpos = fetch_end
            self.chrom = chrom
 #           st_time = time.perf_counter_ns()
 #           self.cache = self.fastafile.fetch(self.chrom, self.start0-110000, self.start0-109997)
 #           end_time = time.perf_counter_ns()
 #           sys.stdout.write("Time for small retrieval=" + str(end_time - st_time) + "\n")  # 15 ms
 #           st_time = time.perf_counter_ns()
            self.cache = self.fastafile.fetch(chrom, self.start0, self.endpos)
 #           end_time = time.perf_counter_ns()
 #           sys.stdout.write("Time for 6MB retrieval="+str(end_time-st_time)+"\n") # 73 milliseconds on mforge for 6MB

        return self.cache[(start-1-self.start0):(endpos-self.start0)]




        # Retrieving the sequence of a genomic region
    def getReference(self, chrom, start, end):
        # Checking if chromosome name exists
        # XXX-HS To make faster, could retrieve large blocks of sequence (2K)
        # .. and cache them .. then next retrieval would be against the cache.
        #
        goodchrom = core.convert_chrom(chrom,self.fastafile.references)
        if goodchrom is None:
            return None
            # Fetching data from reference genome
        if end < start:
            return core.Sequence('')
        if start < 1:
            #start = 1
            # Changed this code .. because calling function needs to be aware of the coordinate violations
            raise Exception("CAVA: getReference: Position requested before first base\n")
        try:
            last = self.fastafile.get_reference_length(goodchrom)
        except:
            try:
                last = self.fastafile.getReferenceLength(goodchrom)
            except:
                ichrom = self.fastafile.references.index(goodchrom)
                if ichrom>=0:
                    last = self.fastafile.lengths()[ichrom]
                else:
                    sys.stderr.write("CAVA:ERROR: Invalid chromosome "+goodchrom)  # should not happen
                    return None


        if end > last:
            #end = last
            # Changed this code .. because calling function needs to be aware of the coordinate violations
            raise Exception("CAVA: getReference: Position requested after last base\n")
#        seq = self.fastafile.fetch(goodchrom, start - 1, end)  # 0-based position

        seq = self.__getseq_from_cache_or_file(goodchrom, start, end)
        return core.Sequence(seq.upper())

#######################################################################################################################
