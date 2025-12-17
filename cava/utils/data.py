#!/usr/bin/env python3


# Classes providing interfaces with annotation databases and the reference genome
#######################################################################################################################

import os
import sys
import importlib.resources

from . import conseq
from . import core
from . import csn
#from CAVA.cava import ensembldb
import re


# import time
import pathlib
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + '/pysamdir')
import pysam

#######################################################################################################################
class Seldata(object):
    #
    # The Seldqta object is meant to be hung off the Transcript object.
    # The transcript object will provide the strand and coordinates
    #
    def __new__(cls, line):
        if line[0].startswith("#") or line[0] == "gene":
            return None
        return super(Seldata,cls).__new__(cls)
    def __init__(self, line):
        self.gene = line[0]
        self.accn = line[3]
        accns = self.accn.split('.')
        self.accn0 = accns[0]
        try:
            self.SECIS_maxstart = int(line[1])
            self.SECIS_maxend = int(line[2])
            self.last_sec_pos = int(line[4])
            self.cds_start = int(line[5])
            self.cds_end = int(line[6])
            self.RNA_len = int(line[7])
            self.last_sel_fromATG = self.last_sec_pos - self.cds_start + 1
            self.cds_len = self.cds_end - self.cds_start + 1
            self.spacing = 111
            self.minspacing = 51

        except:
            sys.stderr.write("ERROR: Invalid format for SECIS config file: line="+"\t".join(line)+"\n")

    def update_maxpos(self,maxpos_in_protein):
        if maxpos_in_protein > self.last_sel_fromATG:
            self.last_sel_fromATG = maxpos_in_protein
    def secis_active(self,cDNA_pos):
        #
        # If tramscript is not an exact match for the version, the only difference is the UTR.. so using
        # CDS-relative position workds.
        # Because Seldata is hung on a transcript, the chromosome shouhld match.. besides, we won't have the chromosome
        # in the SECIS annotation file.
        if cDNA_pos <= self.last_sel_fromATG: # Before known Selenocystering codon .. guaranteed to be translated into Sel
            return 1
        elif cDNA_pos <= self.SECIS_maxstart - self.cds_start +1 - self.spacing: # SECIS most likely active, TAG recoded as Sel(U)
            return 2
        elif cDNA_pos <= self.SECIS_maxstart - self.cds_start +1 - self.minspacing: # SECIS may be active, TAG recoded as Sel(U)
            return 3
        else:
            return 0



# Class representing the Ensembl (transcript) dataset ( can be any database)
class Ensembl(object):
    # Constructor
    def __init__(self, options, genelist, transcriptlist, codon_usage, reference):
        self.options = options
        # Openning tabix file representing the Ensembl database
        self.contigs = dict()
        try:
            self.tabixfile = pysam.Tabixfile(options.args['ensembl'])
        except Exception as e:
            try:
                self.tabixfile = pysam.TabixFile(options.args['ensembl'])
            except Exception as e:
                sys.stderr.write("CAVA: ERROR: error trying to open Tabix file for " + options.args['ensembl'] + "\n")
        for chrom in self.tabixfile.contigs:
            if chrom in reference.reflens:
                self.contigs[chrom] = reference.reflens[chrom]

        self.proteinSeqs = dict()
        self.exonSeqs = dict()
        self.exoncache_hit = dict()
        self.cds_ref = dict()
        self.utr5_ref = dict()
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
        # current gene symbols and old synonyms.
        self.selenogenes =[ 'DIO1', '1733', 'TXDI1', 'THMA2',
                            'DIO2', '1734', 'TXDI2', 'SELENOY', 'SELY', 'DIOII',
                            'DIO3', '1735', 'TXDI3', '5DIII', 'DIOIII',
                            'GPX1', '2876', 'GPXD', 'GSHPX1',
                            'GPX2', '2877', 'GPRP', 'GSHPX-GI', 'GPX-2', 'GI-GPX', 'GPRP-2', 'GPX-GI', 'GSHPX-2',
                            'GPX3', '2878', 'GPX-P', 'GSHPX-3', 'GSHPX-P',
                            'GPX4', '2879', 'MCSP', 'PHGPX', 'SMDS',  'GPX-4', 'GSHPX-4', 'SNGPX', 'SNPHGPX',
                            'GPX6', '257202', 'GPXP3', 'GPX5P', 'GPX-6', 'GSHPX-6', 'DJ1186N24', 'DJ1186N24.1',
                            'SELENOF', '9403', 'SEP15',
                            'SELENOH', '280636', 'C11orf31', 'C17orf10', 'SELH',
                            'SELENOI', '85465', 'SELI', 'EPT1', 'SEPI', 'SPG81',
                            'SELENOK', '58515', 'HSPC030', 'HSPC297', 'SELK',
                            'SELENOM', '140606', 'SELM', 'SEPM',
                            'SELENON', '57190', 'CFTD', 'CMYO3', 'CMYP3', 'MDRS1', 'RSMD1', 'RSS', 'SELN', 'SEPN1',
                            'SELENOO', '83642', 'SELO',
                            'SELENOP', '6414', 'SELP', 'SEPP', 'SEPP1', 'SEP',
                            'SELENOS' '55829', 'AD-015', 'ADO15', 'SBBI8', 'SELS', 'SEPS1', 'VIMP',
                            'SELENOT', '51714', 'SELT',
                            'SELENOU',
                            'SELENOV', '348303', 'SELV',
                            'SELENOW',  '6415', 'SEPW1', 'SELW',
                            'SELENOP1',  'SELENOPZ', 'SEPP1',
                            'SELENOP2', 'SELPB', 'SEPP1L', 'SEPP2',
                            'MSRB1',  '51734', 'SELENOX', 'SELENOR', 'SELR', 'SELX', 'SEPX1', 'SEPR', 'HSOC270',
                            'SEPHS2', '22928',  'SPS2', 'SPS2B',
                            'TXNRD1',  '7296', 'TXNR', 'GRIM-12', 'TRXR1', 'TXNR1', 'TR', 'TR1', 'TRXR1',
                            'TXNRD2', '10587', 'SELZ', 'TR',  'TRXR2',  'TR3', 'TXNR2', 'GCCD5',  'TR',  'TR-BETA',
                            'TXNRD3', '114112',  'TXNRD3NB', 'TXNRD3IT1', 'TR2', 'TRXR3', 'TGR', 'TR2IT1', 'TXNRD3NT1', 'TXNR3', 'TGR']

# Note, SELENOU and SELENOP1 and SELENOP2 are not in humans


        self.selenodata = dict()
        if 'selenogenes' in options.args:
            if len(str(options.args['selenogenes'])) > 1:
                extraselenogenes = str(options.args['selenogenes']).split(",")
                for gs in extraselenogenes:
                    gs = gs.upper()
                    if not gs in self.selenogenes:
                        self.selenogenes.append(gs)

        fid = None
        if 'selenofile' in options.args:
            selenofile = options.args['selenofile']
            if selenofile.startswith("/"):
                try:
                    fid = open(selenofile,'r')
                except IOError:
                    sys.stderr.write("ERROR: Error opening CESIS File="+selenofile+"\n")
            else:
                sys.path.append(str(pathlib.Path().resolve().parents[0]))
                try:
                    fid = importlib.resources.open_text('ensembldb', selenofile)
                except IOError:
                    sys.stderr.write("ERROR: Error opening CESIS File=ensembldb/" + selenofile + "\n")
        else:
            sys.path.append(str(pathlib.Path().resolve().parents[0]))
            fid = importlib.resources.open_text('ensembldb', "SECIS_in_refseq_pos.txt")

        if fid is not None:
            secis_lines = fid.readlines()
            self.load_CESIS(secis_lines)
            fid.close()


    def load_CESIS(self,secis_lines):
        for l in secis_lines:
            line = l.rstrip().split("\t")
            sel = Seldata(line)
            if sel is not None:
                self.selenodata[sel.accn0] = sel
                self.selenodata[sel.accn] = sel  # accession.version

    # Get the list of transcript lines overlapping
    # returns either an iterator over a tabix file .. or a list (that can be iterated over)
    def fetch_overlapping_transcripts(self, chrom, startpos0, endpos1):  # Give 0-base coordinate for start and 1-base for stop
        # If current chromosome is not loaded, then
        #      get tabix iterator figure out length, figu.. and load all chromosomes lines
        # Check self.tabixfile.l
        if chrom not in self.contigs:
            return list()
        if self.loadalltranscripts is False:
            # This fetch consumes most of the runtime for CAVA (78%) as long
            # caching should be faster.
            return self.tabixfile.fetch(reference=chrom, start=startpos0, end=endpos1)
        if self.chrom is None or chrom != self.chrom:
            # Flush cache and load all transcripts.
            self.chrom = chrom
            self.transcript_bins = [None] * (1 + int((1 + self.contigs[chrom]) / self.binsize))
            hits = self.tabixfile.fetch(reference=chrom)
            for line in hits:
                linedat = line.split("\t", 8)
                transcriptid = linedat[0]
                transcriptStart = int(linedat[6])  # lowest coordinate - base 0
                transcriptEnd = int(linedat[7])  # highest coordiate - base 1
                binstart = int(transcriptStart / self.binsize)
                binend = int(transcriptEnd / self.binsize)
                for ebin in range(binstart, binend + 1):
                    if self.transcript_bins[ebin] is None:
                        self.transcript_bins[ebin] = list([[transcriptid, transcriptStart, transcriptEnd, line]])
                    else:
                        (self.transcript_bins[ebin]).append([transcriptid, transcriptStart, transcriptEnd, line])
        lines = list()
        got_transcript = dict()
        binstart = int(startpos0 / self.binsize)
        binend = int(endpos1 / self.binsize)
        for ebin in range(binstart, binend + 1):
            bin_list = self.transcript_bins[ebin]
            if bin_list is not None:
                for bin_content in bin_list:  # [transcriptid, transcriptStart, transcriptEnd, line])
                    if ((startpos0 >= bin_content[1] and endpos1 <= bin_content[2]) or  # variant completely inside
                            (startpos0 >= bin_content[1] and startpos0 + 1 <= bin_content[2]) or  # partial overlap
                            (endpos1 - 1 >= bin_content[1] and endpos1 <= bin_content[2]) or  # partial overlap
                            (startpos0 < bin_content[1] and endpos1 > bin_content[2])):  # Wholly encompassing gene
                        transcriptid = bin_content[0]
                        if transcriptid not in got_transcript:
                            got_transcript[transcriptid] = 1
                            lines.append(bin_content[3])
        return lines

    #
    # Loading transcripts and the exons is very costly (need to read all exons), so caching them save a lot of disk access
    #

    def find_transcript_in_cache_or_in_file(self, line):
        transcriptid = line.split("\t")[0]
        self.nvar += 1
        self.transcript_nvar[transcriptid] = self.nvar
        if transcriptid in self.transcript_cache:
            transcript = self.transcript_cache[transcriptid]
        else:
            # with the reference sequence being cached, fetching a whole transcripts and exons takes 0.03-0.09 ms
            #        ... rather than 150-200 ms if the sequence was not cached.
            transcript = core.Transcript(line)
            if transcript.geneSymbol in self.selenogenes:
                transcript.is_selenocysteine = True
                if transcript.SECIS_data is None:
                    if transcript.TRANSCRIPT in self.selenodata:
                            transcript.SECIS_data = self.selenodata[transcript.TRANSCRIPT]
                    else: #this will work mostly for refseq .. and ensembl if there are ".version" transcripts.
                        transcriptid_split = transcriptid.split('.')
                        if len(transcriptid_split) == 2:
                            if transcriptid_split[0] in self.selenodata:
                                transcript.CESIS_data = self.selenodata[transcriptid_split[0]]

            self.transcript_cache[transcriptid] = transcript
            if len(self.transcript_cache) > self.CACHESIZE:
                vals = list(self.transcript_nvar.values())
                minval = min(vals)
                which_minval = vals.index(minval)
                rm_tr = '' + list(self.transcript_nvar.keys())[which_minval]
                self.transcript_cache.pop(rm_tr)
                self.transcript_nvar.pop(rm_tr)

        return transcript

    # Find transcripts overlapping with a variant

    def findTranscripts(self, variant):
        ret = dict()
        retOUT = dict()

        # Checking chromosome name
        goodchrom = core.convert_chrom(variant.chrom, self.tabixfile.contigs)
        if goodchrom is None:
            return ret, retOUT

        # Defining variant end points.
        # HS: Tabix uses 0-based indexing for start/pos
        if variant.is_insertion is False:
            start = variant.pos - 1
            end = variant.pos + len(variant.ref) - 1
        else:  # for insertion, position got shifted to be after the insertion point .. shift back
            start = variant.pos - 2
            end = variant.pos - 1
        if start < 0:
            start = 0
        if end <= start:
            end = start + 1

        if not variant.is_substitution:  # Insertion/deletion / MNP
            # HS notes:
            # tabix index is loaded in memory, and the data is buffered in 37K (2^16) blocks (from reading the pysam and htslib)
            #   the 'fetch' take about 0.004 - 0.049 ms .. MUCH faster than a single disk seek & read (10-15ms)
            # the iteration take 0.8 ms for the two .. must be cached .. but slow.
            # but uses 70% of the time of CAVA.
            # Loading transcript tables in memory and a pyranges implmementaion would help.
            #  ... so so it is clearly cached .. and no need for further optimization.
            # st_time = time.perf_counter_ns()
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
            hits1 = self.fetch_overlapping_transcripts(goodchrom, start, start + 1)
            hits2 = self.fetch_overlapping_transcripts(goodchrom, end - 1, end)
            # end_time0 = time.perf_counter_ns()
            # sys.stdout.write("Time for two tabix transcripts fetch=" + str(end_time - st_time) + "\n")

            hitdict1 = dict()
            hitdict2 = dict()
            for line in hits1:
                transcript = self.find_transcript_in_cache_or_in_file(line)
                # transcriptStart is 0-based, lowest coordinate .. transcriptEnd is 1-bases highest coordinate, start is 0-based
                if not (transcript.transcriptStart <= start < transcript.transcriptEnd): continue
                # if not strand == transcript.strand: continue
                hitdict1[transcript.TRANSCRIPT] = transcript
            for line in hits2:
                transcript = self.find_transcript_in_cache_or_in_file(line)
                if not (transcript.transcriptStart < end <= transcript.transcriptEnd): continue
                #  if not strand == transcript.strand: continue
                hitdict2[transcript.TRANSCRIPT] = transcript
            # end_time = time.perf_counter_ns()
            # sys.stdout.write("Two tabix fetch & iterate=" + str(end_time0 - st_time) + " "+ str(end_time - st_time) + "\n")

            for key, transcript in hitdict1.items():
                if len(self.genelist) > 0 and transcript.geneSymbol not in self.genelist: continue
                if len(self.transcriptlist) > 0 and transcript.TRANSCRIPT not in self.transcriptlist: continue
                if key in list(hitdict2.keys()):  # e.g. both ends of the variant are in transcript.
                    ret[key] = transcript
                else:  # either paartial overlap or insertion flush with transcript (right end)
                    #                    if variant.is_insertion: # insertions at the edges are TRULY outside, not even partially overlapping
                    retOUT[key] = transcript  # partial overlap downstream of transcript

            #            if not variant.is_insertion:
            for key, transcript in hitdict2.items():  # check for partial overlap upstream of transcript
                if len(self.genelist) > 0 and transcript.geneSymbol not in self.genelist: continue
                if len(self.transcriptlist) > 0 and transcript.TRANSCRIPT not in self.transcriptlist: continue
                if not key in list(hitdict1.keys()):
                    retOUT[key] = transcript

        else:  # Variant is single base Substitution
            hits1 = self.fetch_overlapping_transcripts(goodchrom, start, end)  # self.tabixfile.fetch(region=reg2)
            for line in hits1:
                transcript = self.find_transcript_in_cache_or_in_file(line)

                if len(self.genelist) > 0 and transcript.geneSymbol not in self.genelist: continue
                if len(self.transcriptlist) > 0 and transcript.TRANSCRIPT not in self.transcriptlist: continue

                if not (transcript.transcriptStart + 1 <= end <= transcript.transcriptEnd): continue
                ret[transcript.TRANSCRIPT] = transcript

        return ret, retOUT  # retOUT not populated for substitution

    # Find transcripts overlapping with a variant. This new version will also find variants not fully inside

    def findTranscriptsWide(self, variant):
        ret = dict()
        retOUT = dict()

        # Checking chromosome name
        goodchrom = core.convert_chrom(variant.chrom, self.tabixfile.contigs)
        if goodchrom is None:
            return ret, retOUT

        # Defining variant end points.
        # HS: Tabix uses 0-based indexing for start/pos
        if variant.is_insertion is False:
            start = variant.pos - 1
            end = variant.pos + len(variant.ref) - 1
        else:  # for insertion, position got shifted to be after the insertion point .. shift back
            start = variant.pos - 2
            end = variant.pos - 1
        if start < 0:
            start = 0
        if end <= start:
            end = start + 1

        if not variant.is_substitution:  # Insertion/deletion / MNP
            # HS notes:
            # tabix index is loaded in memory, and the data is buffered in 37K (2^16) blocks (from reading the pysam and htslib)
            #   the 'fetch' take about 0.004 - 0.049 ms .. MUCH faster than a single disk seek & read (10-15ms)
            # the iteration take 0.8 ms for the two .. must be cached .. but slow.
            # but uses 70% of the time of CAVA.
            # Loading transcript tables in memory and a pyranges implmementaion would help.
            #  ... so so it is clearly cached .. and no need for further optimization.
            # st_time = time.perf_counter_ns()
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
            hits1 = self.fetch_overlapping_transcripts(goodchrom, start, start + 1)
            hits2 = self.fetch_overlapping_transcripts(goodchrom, end - 1, end)
            # end_time0 = time.perf_counter_ns()
            # sys.stdout.write("Time for two tabix transcripts fetch=" + str(end_time - st_time) + "\n")

            hitdict1 = dict()
            hitdict2 = dict()
            for line in hits1:
                transcript = self.find_transcript_in_cache_or_in_file(line)
                # transcriptStart is 0-based, lowest coordinate .. transcriptEnd is 1-bases highest coordinate, start is 0-based
                if not (transcript.transcriptStart <= start < transcript.transcriptEnd): continue
                # if not strand == transcript.strand: continue
                hitdict1[transcript.TRANSCRIPT] = transcript
            for line in hits2:
                transcript = self.find_transcript_in_cache_or_in_file(line)
                if not (transcript.transcriptStart < end <= transcript.transcriptEnd): continue
                #  if not strand == transcript.strand: continue
                hitdict2[transcript.TRANSCRIPT] = transcript
            # end_time = time.perf_counter_ns()
            # sys.stdout.write("Two tabix fetch & iterate=" + str(end_time0 - st_time) + " "+ str(end_time - st_time) + "\n")

            for key, transcript in hitdict1.items():
                if len(self.genelist) > 0 and transcript.geneSymbol not in self.genelist: continue
                if len(self.transcriptlist) > 0 and transcript.TRANSCRIPT not in self.transcriptlist: continue

                if key in list(hitdict2.keys()):  # e.g. both ends of the variant are in transcript.
                    ret[key] = transcript
                else:
                    if variant.is_insertion:  # insertions at the edges are TRULY outside, not even partially overlapping
                        retOUT[key] = transcript  # partial overlap downstream of transcript

            if not variant.is_insertion:
                for key, transcript in hitdict2.items():  # check for partial overlap upstream of transcript
                    if len(self.genelist) > 0 and transcript.geneSymbol not in self.genelist: continue
                    if len(self.transcriptlist) > 0 and transcript.TRANSCRIPT not in self.transcriptlist: continue
                    if not key in list(hitdict1.keys()):
                        retOUT[key] = transcript

        else:  # Variant is Substitution
            hits1 = self.fetch_overlapping_transcripts(goodchrom, start, end)  # self.tabixfile.fetch(region=reg2)
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

    # Parse intron CSN coordinates .. only works well from the part before the "dup" string
    #return None id these are not intronic coordinates.
    def getIntronBases(self, x):
        idx = x.find('-')
        if idx < 1:
            idx = x.find('+')
            if idx < 1: return None
            return int(x[idx:])
        return int(x[idx:])
    # Parse Exon CSN coordinates .. only works well from the part before the "dup" string
    # return none, if these are not intronic coordinates
    # e.g. 123-12 or 122+12
    #  New for 2025
    #    edge cases: 5' UTR splice -134-11  or -135+12
    #                3' UTR splice *14+3 *15-5
    # starloc code is repetitive to optimize
    def getExonBases(self, x):
        # remove leading '*' (3') or '-' (5') .. and support remove incorrect leading '+'
        if x[0] == '*' or x[0]=='-' or x[0]=='+':
            x=x[1:]
        idx = x.find('-',1)
        if idx < 1: 
            idx = x.find('+',1)
            if idx < 1: return None
            return int(x[0:idx])
        return int(x[0:idx])



    def getExtractPosOrPosRange(self,x):
        if x is None:
            return None
        if x.startswith('c.'):
            x=x[2:]
        x = x.split('%3B')[0]
        if x.startswith('['):
            x=x[1:]
        if x.find('[')>=0: # repeat inside string
            x=re.sub(r'\[[0-9]+\]','',x)
            if x.endswith(']'): # has to occur after the regex
                x=x[0:len(x)-1]
        if x.endswith(']'): # bracket around a non-repeat
            x=x[0:len(x)-1]
        m=re.match(r'^([0-9\-\+_\*]+)ins.*inv',x)
        if m:
            x=m.group(1)

        x = re.sub(r'[A-Za-z]', '', x)  # remove bases and keywords like ins del or dup
        return x


    # argument is a partially parsed csn string
    # x as in c.x_p.proteinHGVS
    def parseRep(self,csnpart):
        # This is for a repeat that is not an insertion.
        # [position_positionGC[4]]%3G[position_positionGC[6]]
        # [position_-11position-4GC[4]]%3G[position-11_position-4GC[6]]

        if csnpart is None or csnpart.find('%3B')==-1:
            return None
        [x,y] = csnpart.split('%3B')
        if x.startswith('['):
            x=x[1:]
        if x.find('[')>=0: # repeat inside string
            m = re.search(r'^([0-9_+\-]+)([A-Z,a-z]+)\[([0-9]+)\]',x)
            if m:
                [crange,repeat_seq,n_ref] = [m.group(1),m.group(2),m.group(3)]
            else:
                return [None,None,None,None]
        else:
            return [None, None, None, None]
        if y.startswith('['):
            y=y[1:]
        if y.find('[')>=0: # repeat inside string
            m = re.search(r'^([0-9_+\-]+)([A-Z,a-z]+)\[([0-9]+)\]',y)
            if m:
                [crange,repeat_seq,n_alt ] = [m.group(1),m.group(2),m.group(3)]
                return [crange,repeat_seq,int(n_ref),int(n_alt)]
            else:
                return [None,None,None,None]
        else:
            return [None, None, None, None]
        return [x,None,None,None]


    # Check if variant is duplication overlapping SS boundary
    # ssrange is the size of the intron region (default 8, possible values should be >6)
    # This function needs to return True only if the duplication intervaal
    # touches the edge of the ss region and is inside.. so that class can be corrected
    # from SS to INT.
    # if the Duplication interval does not meet the creteria for shifting out, return False.
    #
    # XXX-HS This function is not 100 percent. If a repeat is annotated as a Dup, it is possible
    # that the repeat could extend beyond the dup region. This yields a conservative result
    #

    def isDupOverlappingSSBoundary(self, csnval, ssrange=8):
        if ssrange == 0 or csnval is None or len(csnval)==0:
            return False  # Nothing to correct, since nothing will overlap the splice region.
        if '_p' in csnval:
            [cpart, _] = csnval.split('_p')
        else:
            cpart = csnval
        if cpart.startswith("c."):
            cpart=cpart[2:]
        if cpart.find("del")>=0:
            return False
        isDup = False
        isIns = False
        if cpart.find('ins') !=-1:
            isIns = True
        elif cpart.find('dup') != -1:
            isDup = True
        elif cpart.find('[') != -1:
            isitoverlappingSSboundary = self.isRepOverlappingSSBoundary(cpart,ssrange)
            return isitoverlappingSSboundary
        else:
            return False

        cpart = self.getExtractPosOrPosRange(cpart)

        if '_' in cpart:
            xsys = cpart.split('_')
            if len(xsys)>2:
                sys.stderr.write("CAVA error: Too many underscores in interval :"+cpart+" from CSN_range="+csnval+"\n")
            [xs, ys] = xsys[0:2]
        else: # Single-base duplication
            x = self.getIntronBases(cpart)
            if x is None: return False
            return x == ssrange or x == -ssrange # dup at the edge of the region, splicosome will not see it as being inside

        x = self.getIntronBases(xs)
        y = self.getIntronBases(ys)
#        if x is None or y is None: return False
        # pre-2025, this was not correct two ways:
        #      Duplication can only shift out if one edge of the dup is
        #      either right at the SS boundary or overlapping it .. (and the other end does not have to be in splice regioin)

        if x is None and y is None: # No possibility of shifting
            return False
        xe= self.getExonBases(xs)
        ye=self.getExonBases(ys)
        if isIns is True:
            if x is not None:
                if x==(-ssrange-1):
                    return True
                if x == ssrange:
                    return True
            return False # Nothing to rescue.
        if xe is None or ye is None or xe!=ye: # Duplication across different exons or might cover whole ss area
            return False
        elif isDup is True:
            # Next few lines don't do the 'proper' thing for a very short intron ..
            #   they will not rescue events spanning the middle of the intron.(coordinate change)
            # .. but we'll assume that very small introns will get disrupted by any sequence insertioin, so we don't want to rescur anything
            if xe is not None and ye is not None:
                if xe != ye: # two ends of dup in different exons, do not rescue.
                    return False
                elif x>ssrange or y< -ssrange: # Don't change anything. outside intron SS regions
                    return False
                elif x==ssrange or y == -ssrange:
                    return True
            if x is None or y is None or xe is None or ye is None: # One of the ends of the dup, crosses the intron/exon boundary, so if it were to be
                return False             # to be rescued, it would create an entire splice site, likely has some functional impact .. so do not rescue
            if y> -ssrange: # Dup will insert before acceptor site, doesn't matter if y is in intron.
                if x<= -ssrange:
                    return True
                return False
            if x < ssrange:  # Dup will insert before acceptor site, doesn't matter if y is in intron.
                if y >= ssrange:
                    return True
                return False
            return False
        return False

        # Now deal with isRep

#pre102025        if x is None or y is None: return False
#pre102025        return self.inrange(x, y, ssrange) or self.inrange(x, y, -ssrange)


    # input is partially parsed csn, e.g. c.cpart_p.proteinHGVS
    # Returns true if could shift an event targetting the splice range to outside
    # False, means to leave class as-is (either in intron or in SS)
    #
    def isRepOverlappingSSBoundary(self, cpart, ssrange=8):
        [cpart, repeat_seq,n_ref,n_alt] = self.parseRep(cpart)
        if repeat_seq is None:
            return False
        repeat_len = len(repeat_seq)
        if '_' in cpart:
            [xs, ys] = cpart.split('_')
        else:  # Single base repeat with one copy
               # Cannot be a deletion .. has to be a repeat expansion..
               # can we insert the repeated bases outside the splice region.
            x = self.getIntronBases(cpart)
            if x is None: return False
            return x == ssrange or x == -ssrange  # 3 or more repeats at the edge of the region, splicosome will not see it as being inside

        x = self.getIntronBases(xs)
        y = self.getIntronBases(ys)
        #        if x is None or y is None: return False
        # pre-2025, this was not correct two ways:
        #      Duplication can only shift out if one edge of the dup is
        #      either right at the SS boundary or overlapping it .. (and the other end does not have to be in splice regioin)

        if x is None and y is None:  # No possibility of shifting solely within ssregion.
            return False
        xe = self.getExonBases(xs)
        ye = self.getExonBases(ys)
        if xe is None or ye is None: # One of them in exon, no rescue.
            return False
        if n_ref>n_alt: # Deletion. need at least (n_ref-n_alt)*repeat_len past edge of ssboundary.
            keepbases = (n_ref-n_alt)*repeat_len
            if y<0 and x<0 and y>=-ssrange and  x+(keepbases-1) < -ssrange:
                return True
            if x>0 and y>0 and x<=ssrange and y-(keepbases-1) > ssrange:
                return True
        elif n_ref<n_alt: # Insertion, can insert anywhere as long as art of repeat is past ssrange.
            if x<= -ssrange or y>= ssrange:
                return True
        return False


    # Correct CLASS annotations for duplications overlapping SS boundary because the splicing biology
    # will no lomger use the 3' shifting rule ..and the dup will be read as being in the intron
    # .. though it could still affect the ESS

    def correctClasses(self, mycsn, class_plus, class_minus):
        if self.isDupOverlappingSSBoundary(mycsn, ssrange=int(self.options.args['ssrange'])):
            if class_plus == 'SS' and class_minus == 'INT': return 'INT', 'INT'
            if class_plus == 'INT' and class_minus == 'SS': return 'INT', 'INT'
        return class_plus, class_minus

    # Correct SO annotations for duplications overlapping SS boundary
    def correctSOs(self, mycsn, so_plus, so_minus):
        if self.isDupOverlappingSSBoundary(mycsn):
            if so_plus == 'intron_variant|splice_region_variant' and so_minus == 'intron_variant': return 'intron_variant', 'intron_variant'
            if so_plus == 'intron_variant' and so_minus == 'intron_variant|splice_region_variant': return 'intron_variant', 'intron_variant'
        return so_plus, so_minus

        # Annotating a variant based on Ensembl data

    def annotate(self, variant, reference, impactdir):
        # Create left-aligned and right-aligned versions of the variant
        if variant.is_deletion or variant.is_insertion:  # optimization, MNP or substitutions cannot be shifted
            variant_plus = variant.alignOnPlusStrand(reference)
            variant_minus = variant.alignOnMinusStrand(reference)
            # Checking if variant alignment makes any difference
            if variant_plus.pos == variant_minus.pos:
                difference = False
            else:
                difference = True
        else:
            variant_plus = variant
            variant_minus = variant
            difference = False
        # If givealt is True, then will evaluate csn for the variant in the 'wrong/opposite' left-right shifting
        givealt=False
        if 'givealt' in self.options.args and self.options.args['givealt'] is True:
            givealt=True

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
        CSNALTstring = ''

        # Collecting transcripts that overlap with the variant
        # transcript overlapping left-end and right end of variant ("OUT" is empty for single-bases)
        transcripts_plus, transcriptsOUT_plus = self.findTranscripts(variant_plus)
        transcripts_minus, transcriptsOUT_minus = self.findTranscripts(variant_minus)

        if variant.is_deletion or variant.is_complex:
            # If variant is Deletion(or Del+repl),
            #         being partial at 3' end has predictable functional impact
            #         Deletion 5' end should be annotated because CAP site is not always correct)
            #         we want the CSN produced
            transcripts_allplus = set(list(transcripts_plus.keys()) + list(transcriptsOUT_plus.keys()))
            transcripts_allminus = set(list(transcripts_minus.keys()) + list(transcriptsOUT_minus.keys()))
            transcripts = set(list(transcripts_plus.keys()) + list(transcripts_minus.keys()) + list(
                transcriptsOUT_plus.keys()) + list(transcriptsOUT_minus.keys()))
            transcriptsOUT = set()
        else:  # This will only lead to an "OUT" for insertions right at the edge.
            transcripts_allplus = set(list(transcripts_plus.keys()))
            transcripts_allminus = set(list(transcripts_minus.keys()))
            transcripts = set(list(transcripts_plus.keys()) + list(transcripts_minus.keys()))
            transcriptsOUT = set(list(transcriptsOUT_plus.keys()) + list(transcriptsOUT_minus.keys()))  # SNP has no out

        transcripts = sorted(list(transcripts))
        transcriptsOUT = sorted(list(transcriptsOUT))
        transcripts_allplus = sorted(transcripts_allplus)
        transcripts_allminus = sorted(transcripts_allminus)

        # Annotating with transcripts that only partially overlap with the variant
        for TRANSCRIPT in transcriptsOUT:  # either plus or minus shifted variant is partially outside
            if TRANSCRIPT in transcripts: continue  # skip, if want the CSN annotated
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
                CSNALTstring += ':'

            TRANSCRIPTstring += TRANSCRIPT
            GENEstring += transcript.geneSymbol
            GENEIDstring += transcript.geneID
            TRINFOstring += transcript.TRINFO

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
            CSNALTstring += '.'

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
                CSNALTstring += ':'

            # Creating the TRANSCRIPT, GENE and TRINFO annotations
            TRANSCRIPTstring += TRANSCRIPT
            GENEstring += transcript.geneSymbol
            GENEIDstring += transcript.geneID
            TRINFOstring += transcript.TRINFO

            # Creating the LOC annotation
            if TRANSCRIPT in list(transcripts_plus.keys()):
                loc_plus = transcript.whereIsThisVariant(variant_plus)
            elif TRANSCRIPT in list(transcriptsOUT_plus.keys()):
                # loc_plus = 'OUT'
                loc_plus = transcript.whereIsThisVariant(variant_plus)
            else:
                loc_plus = '.'

            if difference is True:
                if TRANSCRIPT in list(transcripts_minus.keys()):
                    loc_minus = transcript.whereIsThisVariant(variant_minus)
                elif TRANSCRIPT in list(transcriptsOUT_minus.keys()):
                    loc_minus = 'OUT'
                    loc_minus = transcript.whereIsThisVariant(variant_minus)
                else:
                    loc_minus = '.'
            else:
                loc_minus = loc_plus

            # Creating reference and mutated protein sequence
            #  not_exonic_plus= means outside transcript .. which can (now) happen for deletion or MNP
            # Cannot call protein when variant spans intron/exon boundary, because they affect splicing
            #   so any predicted protein would be very uncertain .. bettter to trigger a p.?
            #
            # XXX-HS need to add Start Gain in 5'UTR.. from SNP/insertion/deletions
            #
            # notexonic_plus = (
            #         ( loc_plus in ['5UTR','3UTR, '<-5UTR','3UTR->','-','OUT','.']) or ('In' in loc_plus) )
            # if loc_plus.startswith("5UTR-) and variant.is_insertion is False: # Variants spanning start site.
            #    notexonic_plus = False
            notexonic_plus = transcript.isOutsideTranslatedRegion(variant_plus)
            in_utr5_plus = False
            in_utr5_minus = False
            start_aa = "M"
            is_methionine = True
            if variant.is_insertion:
                if transcript.isPositionOutsideCDS_5prime(variant_plus.pos) or \
                        (transcript.strand == 1 and variant_plus.pos == transcript.transcriptStart + 1):
                    in_utr5_plus = True
            elif transcript.isPositionOutsideCDS_5prime(variant_plus.pos) and transcript.isPositionOutsideCDS_5prime(
                    variant_plus.pos + len(variant_plus.alt) - 1):
                in_utr5_plus = True

            if difference is True:
                # notexonic_minus = (
                #        ('5UTR' == loc_minus) or ('3UTR' == loc_minus) or ('-' in loc_minus) or (
                #        'In' in loc_minus) or (loc_minus == 'OUT') or (loc_minus == '.'))  # Variants overlapping the ends of the protein
                # if loc_minus.startswith("5UTR-Ex") and variant.is_insertion is False: # variants Spanning Start site.
                #    notexonic_minus = False
                notexonic_minus = transcript.isOutsideTranslatedRegion(variant_minus)
                if variant.is_insertion:
                    if transcript.isPositionOutsideCDS_5prime(variant_minus.pos) or \
                            (transcript.strand == -1 and variant_minus.pos == transcript.transcriptEnd + 1):
                        in_utr5_minus = True
                elif transcript.isPositionOutsideCDS_5prime(
                        variant_minus.pos) and transcript.isPositionOutsideCDS_5prime(
                    variant_minus.pos + len(variant_plus.alt) - 1):
                    in_utr5_minus = True
            else:
                notexonic_minus = notexonic_plus
                in_utr5_minus = in_utr5_plus
            exonseqs = None
            utr5_ref = ''
            utr5_mut = ''

            if in_utr5_plus is True or in_utr5_minus is True or notexonic_plus is False or notexonic_minus is False:
                self.cache_num += 1
                self.exoncache_hit[transcript.TRANSCRIPT] = self.cache_num
                if not transcript.TRANSCRIPT in list(self.proteinSeqs.keys()):
                    protein, exonseqs, cds_ref, utr5_ref = transcript.getProteinSequence(reference, None, None,
                                                                                         self.codon_usage)
                    self.proteinSeqs[transcript.TRANSCRIPT] = protein
                    self.exonSeqs[transcript.TRANSCRIPT] = exonseqs
                    self.cds_ref[transcript.TRANSCRIPT] = cds_ref
                    transcript.exonseqs = exonseqs
                    self.utr5_ref[transcript.TRANSCRIPT] = utr5_ref
                    if len(self.proteinSeqs) > self.CACHESIZE:  # Cache of proteins and exons data
                        vals = list(self.exoncache_hit.values())
                        minval = min(vals)
                        which_minval = vals.index(minval)
                        rm_tr = '' + list(self.exoncache_hit.keys())[which_minval]
                        self.exonSeqs.pop(rm_tr)
                        self.exoncache_hit.pop(rm_tr)
                        self.proteinSeqs.pop(rm_tr)
                        self.cds_ref.pop(rm_tr)
                        self.utr5_ref.pop(rm_tr)
                else:
                    protein = self.proteinSeqs[transcript.TRANSCRIPT]
                    exonseqs = self.exonSeqs[transcript.TRANSCRIPT]
                    cds_ref = self.cds_ref[transcript.TRANSCRIPT]
                    utr5_ref = self.utr5_ref[transcript.TRANSCRIPT]
                if len(protein) > 0 and protein[0] != 'M':  # Support alternate start codon at the price of not supporting partial as well.
                    # Support the following alt start codons AUA,AUU,CUG,GUG,ACG,UUG,AUC,AAG,AGG
                    codon = cds_ref[0] + cds_ref[1] + cds_ref[2]
                    if codon in ['ATA', 'ATT', 'CTG', 'GTG', 'ACG', 'TTG', 'ATC', 'AAG', 'AGG']:
                        start_aa = protein[0]
                        is_methionine = False
                        protein = "M" + protein[1:]
            else:
                protein = ''
                cds_ref = ''
                utr5_ref = ''

            if (transcript.strand == -1 and  givealt is False)  or (notexonic_plus is True and in_utr5_plus is False): # No chance to have a protein or even an alt-start
                mutprotein_plus = None
                utr5_plus = None
            else:
                mutprotein_plus, exonseqsalt_plus, cds_mut_plus, utr5_plus = transcript.getProteinSequence(reference,
                                                                                                           variant_plus,
                                                                                                           exonseqs,
                                                                                                           self.codon_usage)
                # support alternate start codon .. but don't support one alternate mutating to another alternate start codon, just
                # because these kinds of change lead to a large change in expression and these would be invible with a synonymous Methionine->Methionine mutation
                # We support alternate start codon, because we respect non-partial annotation.
                if mutprotein_plus is not None and is_methionine is False and len(mutprotein_plus)>0:
                    if mutprotein_plus[0] == start_aa and cds_mut_plus[0]==cds_ref[0] and cds_mut_plus[1]==cds_ref[1] and cds_mut_plus[2]==cds_ref[2]:
                        mutprotein_plus='M'+mutprotein_plus[1:]
                if mutprotein_plus is not None and (transcript.geneSymbol in self.selenogenes) and cds_ref is not None:
                    mutprotein_plus = transcript.trimSelenoCysteine(cds_ref, cds_mut_plus, protein, mutprotein_plus,
                                                                    variant_plus,transcript)

            if (transcript.strand == 1 and  givealt is False) or (notexonic_minus is True and in_utr5_minus is False):
                mutprotein_minus = None
                utr5_minus = None
            else:
                mutprotein_minus, exonseqsalt_minus, cds_mut_minus, utr5_minus = transcript.getProteinSequence(
                    reference, variant_minus, exonseqs, self.codon_usage)
            # support alternate start codon .. but don't support one alternate mutating to another alternate start codon
                if mutprotein_minus is not None and is_methionine is False and len(mutprotein_minus)>0:
                    if mutprotein_minus[0] == start_aa and cds_mut_minus[0]==cds_ref[0] and cds_mut_minus[1]==cds_ref[1] and cds_mut_minus[2]==cds_ref[2]:
                        mutprotein_minus='M'+mutprotein_minus[1:]
                if mutprotein_minus is not None and transcript.geneSymbol in self.selenogenes and cds_ref is not None:
                    mutprotein_minus = transcript.trimSelenoCysteine(cds_ref, cds_mut_minus, protein,
                                                                     mutprotein_plus, variant_minus,transcript)

            # Creating the CSN annotations both for left and right aligned variant
            if TRANSCRIPT in transcripts_allplus and (givealt is True or transcript.strand == 1):
                csn_plus, protchange_plus,csn_plus_alt = csn.getAnnotation(variant_plus, transcript, reference, protein,
                                                              mutprotein_plus)
                csn_plus_str = csn_plus.getAsString()
                if csn_plus_alt is not None:
                    csn_plus_alt_str = csn_plus_alt.getAsString()
                else:
                    csn_plus_alt_str = '.'

            else:
                csn_plus_str, protchange_plus,csn_plus_alt_atr = '.', ('.', '.', '.'),'.'
            if TRANSCRIPT in transcripts_allminus and (givealt is True or transcript.strand == -1):
                csn_minus, protchange_minus,csn_minus_alt = csn.getAnnotation(variant_minus, transcript, reference, protein,
                                                                    mutprotein_minus)
                csn_minus_str = csn_minus.getAsString()
                if csn_minus_alt is not None:
                    csn_minus_alt_str = csn_minus_alt.getAsString()
                else:
                    csn_minus_alt_str = '.'
            else:
                csn_minus_str, protchange_minus,csn_minus_alt_str = '.', ('.', '.', '.'), '.'

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
                                                           int(self.options.args['ssrange']),
                                                           reference, exonseqs, utr5_ref, utr5_plus)
                else:
                    class_plus = '.'

                if TRANSCRIPT in transcripts_allminus:
                        class_minus = conseq.getClassAnnotation(variant_minus, transcript, protein, mutprotein_minus,
                                                                loc_minus, int(self.options.args['ssrange']),
                                                                reference, exonseqs, utr5_ref, utr5_minus)
                else:
                        class_minus = '.'

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

                if TRANSCRIPT in transcripts_allminus:
                    so_minus = conseq.getSequenceOntologyAnnotation(variant_minus, transcript, protein,
                                                                    mutprotein_minus, loc_minus)
                else:
                    so_minus = '.'

            # Deciding which is the correct CSN and CLASS annotation
            if transcript.strand == 1:
                class_plus, class_minus = self.correctClasses(csn_plus_str, class_plus, class_minus)
                so_plus, so_minus = self.correctSOs(csn_plus_str, so_plus, so_minus)
                if class_plus in ["ESS", 'SS5', 'SS'] or so_plus.startswith("splice"):
                    csn_plus_arr = csn_plus_str.split("_p.")
                    if len(csn_plus_arr) == 1 or csn_plus_arr[1] in ["=", "(=)", "="]:
                        csn_plus_str = csn_plus_arr[0] + "_p.?"
                if class_minus in ["ESS", 'SS5', 'SS'] or so_minus.startswith("splice"):
                    csn_minus_arr = csn_minus_str.split("_p.")
                    if len(csn_minus_arr) == 1 or csn_minus_arr[1] in ["=", "(=)", "="]:
                        csn_minus_str = csn_minus_arr[0] + "_p.?"
                if class_plus == "IG":
                    so_plus = so_plus + '|' + 'initiation_gain_in_5utr'
                if class_minus == "IG":
                    so_minus = so_minus + '|' + 'initiation_gain_in_5utr'
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
                CSNALTstring += csn_plus_alt_str
            else:
                class_plus, class_minus = self.correctClasses(csn_minus_str, class_plus, class_minus)
                so_plus, so_minus = self.correctSOs(csn_minus_str, so_plus, so_minus)
                if class_plus in ["ESS", 'SS5', 'SS'] or so_plus.startswith("splice"):
                    csn_plus_arr = csn_plus_str.split("_p.")
                    if len(csn_plus_arr) == 1 or csn_plus_arr[1] in ["=", "(=)"]:
                        csn_plus_str = csn_plus_arr[0] + "_p.?"
                if class_minus in ["ESS", 'SS5', 'SS'] or so_minus.startswith("splice"):
                    csn_minus_arr = csn_minus_str.split("_p.")
                    if len(csn_minus_arr) == 1 or csn_minus_arr[1] in ["=", "(=)"]:
                        csn_minus_str = csn_minus_arr[0] + "_p.?"
                if class_plus == "IG":
                    so_plus = so_plus + '|' + 'initiation_gain_in_5utr'
                if class_minus == "IG":
                    so_minus = so_minus + '|' + 'initiation_gain_in_5utr'
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
                CSNALTstring += csn_minus_alt_str

            if  givealt is True:
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

            if givealt is True:
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
        variant.addFlag('CSNALT', CSNALTstring)

        if self.options.args['ontology'].upper() in ['CLASS', 'BOTH']: variant.addFlag('CLASS', CLASSstring)
        if self.options.args['ontology'].upper() in ['SO', 'BOTH']: variant.addFlag('SO', SOstring)

        if not impactdir is None: variant.addFlag('IMPACT', IMPACTstring)

        if givealt is True:
            variant.addFlag('ALTANN', ALTANNstring)
            if self.options.args['ontology'].upper() in ['CLASS', 'BOTH']: variant.addFlag('ALTCLASS', ALTCLASSstring)
            if self.options.args['ontology'].upper() in ['SO', 'BOTH']: variant.addFlag('ALTSO', ALTSOstring)

        if (givealt is True and ('givealtflag' in self.options.args and self.options.args['givealtflag'] is True)):
            variant.addFlag('ALTFLAG', ALTFLAGstring)

        return variant


#######################################################################################################################


# Class representing the dbSNP dataset
class DbSnp(object):
    # Constructor
    def __init__(self, options):
        # Opening tabix file representing the dbSNP database
        self.tabixfile = pysam.Tabixfile(options.args['dbsnp'])

    # Annotating a variant based on dbSNP data
    def annotate(self, variant):
        # Checking if variant is a SNP at all
        if variant.is_substitution:
            # Fetching data from dbSNP database
            goodchrom = core.convert_chrom(variant.chrom, self.tabixfile.contigs)
            if goodchrom is None:
                variant.addFlag('DBSNP', '')
                return
            reg = goodchrom + ':' + str(variant.pos) + '-' + str(variant.pos)
            lines = self.tabixfile.fetch(region=reg)
            for line in lines:
                cols = line.split('\t')
                # Adding DBSNP annotation to the variant
                alts = cols[3].split(',')
                if variant.alt in alts:
                    variant.addFlag('DBSNP', cols[0])
            if not 'DBSNP' in variant.flags:
                variant.addFlag('DBSNP', '')
        else:
            variant.addFlag('DBSNP', '')
        return


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
        for i in range(0, len(lengths)):
            chrom = str(references[i])
            self.reflens[chrom] = lengths[i]
            self.reflens[chrom.upper()] = lengths[i]
            self.reflens[chrom.lower()] = lengths[i]
            if chrom.startswith("chr") and chrom.find('_') == -1:
                self.reflens[chrom[3:]] = lengths[i]
                self.reflens[chrom[3:].upper()] = lengths[i]
                self.reflens[chrom[3:].lower()] = lengths[i]
            else:
                self.reflens["chr" + chrom] = lengths[i]
                self.reflens["CHR" + chrom.upper()] = lengths[i]
                self.reflens["chr" + chrom.upper()] = lengths[i]
                self.reflens["chr" + chrom.lower()] = lengths[i]

            if chrom in ["chrM", "chrMT", "M", "MT", "chrM"]:
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
        if chrom != self.chrom or start < self.start0 + 1 or endpos > self.endpos:
            fetch_start = start - self.PAD
            if fetch_start < 1:
                fetch_start = 1
            fetch_end = start + self.CACHESIZE
            if fetch_end < endpos:
                fetch_end = endpos + self.PAD
            if fetch_end > self.reflens[chrom]:
                fetch_end = self.reflens[chrom]
            self.start0 = fetch_start - 1
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

        return self.cache[(start - 1 - self.start0):(endpos - self.start0)]

        # Retrieving the sequence of a genomic region

    def getReference(self, chrom, start, end):
        # Checking if chromosome name exists
        # XXX-HS To make faster, could retrieve large blocks of sequence (2K)
        # .. and cache them .. then next retrieval would be against the cache.
        #
        goodchrom = core.convert_chrom(chrom, self.fastafile.references)
        if goodchrom is None:
            return None
            # Fetching data from reference genome
        if end < start:
            return core.Sequence('')
        if start < 1:
            # start = 1
            # Changed this code .. because calling function needs to be aware of the coordinate violations
            raise Exception("CAVA: getReference: Position requested before first base\n")
        try:
            last = self.fastafile.get_reference_length(goodchrom)
        except:
            try:
                last = self.fastafile.getReferenceLength(goodchrom)
            except:
                ichrom = self.fastafile.references.index(goodchrom)
                if ichrom >= 0:
                    last = self.fastafile.lengths[ichrom]
                else:
                    sys.stderr.write("CAVA:ERROR: Invalid chromosome " + goodchrom)  # should not happen
                    return None

        if end > last:
            # end = last
            # Changed this code .. because calling function needs to be aware of the coordinate violations
            raise Exception("CAVA: getReference: Position requested after last base\n")
        #        seq = self.fastafile.fetch(goodchrom, start - 1, end)  # 0-based position

        seq = self.__getseq_from_cache_or_file(goodchrom, start, end)
        return core.Sequence(seq.upper())

#######################################################################################################################
