#!/usr/bin/env python3


# Collection of basic classes and functions
#######################################################################################################################


import gzip
import logging
import os
import sys
import time
import re
from . import csn


#######################################################################################################################

# Class representing a single variant call
class Variant(object):
    # Constructor
    def __init__(self, chrom, pos, vcf_ref, vcf_alt):
        vcf_ref = vcf_ref.upper()
        vcf_alt = vcf_alt.upper()
        self.chrom = chrom
        self.id = chrom + ":" + str(pos) + ":" + vcf_ref + "/" + vcf_alt
        self.pos = pos
        self.vcf_padded_base = ''
        shift_pos, a, b = self.trimCommonStart(vcf_ref, vcf_alt)
        self.pos = self.pos + shift_pos
        shiftpos2, x, y = self.trimCommonEnd(a, b)
        if (len(x) == 0 or len(y) == 0) and shift_pos > 0:
            self.vcf_padded_base = vcf_ref[shift_pos - 1]

        self.ref = Sequence(x)
        self.alt = Sequence(y)
        self.flags = []
        self.flagvalues = []
        self.alignonplus = None  # Cache the result of alignment because it is expensive and can be done many times.
        self.alignonminus = None
        lx = len(x)
        ly = len(y)
        if lx==1 and ly ==1:
            self.is_substitution = True
            self.is_deletion = False
            self.is_insertion = False
            self.is_complex = False
            self.is_in_frame = True
        elif lx==0 and ly>0:
            self.is_substitution = False
            self.is_deletion = False
            self.is_insertion = True
            self.is_complex = False
            self.is_in_frame = ((lx -ly) % 3 == 0)
        elif lx >0 and ly == 0:
            self.is_substitution = False
            self.is_deletion = True
            self.is_insertion = False
            self.is_complex = False
            self.is_in_frame = ((lx - ly) % 3 == 0)
        else:
            self.is_substitution = False
            self.is_deletion = False
            self.is_insertion = False
            self.is_complex = True
            self.is_in_frame = ((lx - ly) % 3 == 0)
        self.left_result = None
        self.right_result = None
        self.full_reult = None
        self.seleno_early_stop = sys.maxsize


    # Getting basic information about variant
    def info(self):
        return self.getFullPosition() + '_' + self.ref + '>' + self.alt

    # Getting chromosome name and position combined
    def getFullPosition(self):
        return self.chrom + ':' + str(self.pos)

    # Checking if the variant is a base substitution
    def isSubstitution(self):
        return len(self.ref) == 1 and len(self.alt) == 1

    # Checking if the variant is an insertion
    #def isInsertion(self):
    #    return len(self.ref) == 0 and len(self.alt) > 0

    # Checking if the variant is a deletion
    #def isDeletion(self):
    #    return len(self.ref) > 0 and len(self.alt) == 0

    # Checking if the variant is a complex indel (pre-computed.. not used)
    def isComplex(self):
        return len(self.ref) > 0 and len(self.alt) > 0 and not self.is_substitution


    # Checking if the variant is in-frame (pre-computed .. not used)
    def isInFrame(self):
        return (len(self.alt) - len(self.ref)) % 3 == 0

    # Checking if the variant overlaps with a genomic region
    def overlap(self, start, end):
        if self.is_insertion:
            return (self.pos - 1 >= start) and (self.pos <= end)
        if self.pos > start:
            return end >= self.pos
        if self.pos == start:
            return True
        if self.pos + len(self.ref) - 1 >= start:
            return self.pos + len(self.ref) - 1 <=end
        if self.pos<start and self.pos+len(self.ref)-1>end:
            return True
        return False

    # Getting the value of a given annotation flag
    def getFlag(self, flag):
        return self.flagvalues[self.flags.index(flag)]

    # Adding annotation flag to variant
    def addFlag(self, flag, value):
        self.flags.append(flag)
        self.flagvalues.append(value)

    # Annotating variant
    def annotate(self, ensembl, dbsnp, reference, impactdir):
        self.annotateWithType()
        if ensembl is not None:
            self = ensembl.annotate(self, reference, impactdir)
        if dbsnp is not None:
            dbsnp.annotate(self)

    # Annotating variant with its type
    def annotateWithType(self):
        if self.is_substitution:
            self.addFlag('TYPE', 'Substitution')
        elif self.is_insertion:
            self.addFlag('TYPE', 'Insertion')
        elif self.is_deletion:
            self.addFlag('TYPE', 'Deletion')
        elif self.is_complex:
            self.addFlag('TYPE', 'Complex')

    def remove_same_bases_ref_alt(self):  # Only apply for SNV or MNP with padded bases to remove.
        if self.is_substitution or self.is_insertion or self.is_deletion:
            return [self.pos, self.ref, self.alt]
        if self.ref == self.alt:
            return [self.pos, self.ref[0], self.alt[0]]
        left, seq1, seq2 = self.trimCommonStart(self.ref, self.alt)
        right, newref, newalt = self.trimCommonEnd(seq1, seq2)
        return [self.pos + left, newref, newalt]

    # Aligning variant on the GENOMIC plus strand, shifting 3' (right) on DNA (NOT cDNA)
    def alignOnPlusStrand(self, reference):  # right shift, so appropriate for transcripts on plus strand.
        if self.alignonplus is not None:
            return self.alignonplus
        if self.chrom not in reference.reflens or not (self.is_insertion or self.is_deletion):
            [self.pos, self.ref, self.alt] = self.remove_same_bases_ref_alt()  # Trime redundant bases from SNP/MNP
            self.alignononplus = self
            return self
        reflen = reference.reflens[self.chrom]
        maxreplen = abs(len(self.ref)-len(self.alt))
        PADLEN = max(100,1+5*maxreplen) # the bigger the repeat, the less likely there are multiple exact copies, so 5 is OKfSequence
#        seq1_0 = reference.getReference(self.chrom, self.pos, self.pos + len(self.ref) - 1 + PADLEN)
#        seq2_0 = self.alt + reference.getReference(self.chrom, self.pos + len(self.ref),
#                                                 self.pos + len(self.ref) + PADLEN - 1)
        maxpos = self.pos + len(self.ref) - 1 + PADLEN
        if maxpos > reflen:
            maxpos = reflen
        EXTRA_PREBASE = 1
        if self.pos == 1:
            EXTRA_PREBASE = 0
        seq12_0 = reference.getReference(self.chrom, self.pos-EXTRA_PREBASE, maxpos)
        seq1_0 = seq12_0[EXTRA_PREBASE:]
        seq2_0 = self.alt+seq12_0[len(self.ref)+EXTRA_PREBASE:]
        left, seq1, seq2 = self.rightAlign(seq1_0, seq2_0)
        while left >= PADLEN - (maxreplen -  1) and maxpos!=reflen:  # Shifted almost all the way to the end, make sure there is not more to shift
            PADLEN += PADLEN
            maxpos = self.pos + len(self.ref) - 1 + PADLEN
            if maxpos > reflen:
                maxpos = reflen
#            seq1_0 = reference.getReference(self.chrom, self.pos, maxpos)
#            seq2_0 = self.alt + reference.getReference(self.chrom, self.pos + len(self.ref),
#                                                       self.pos + len(self.ref) + PADLEN - 1)
            seq12_0 = reference.getReference(self.chrom, self.pos-EXTRA_PREBASE, maxpos)
            seq1_0 = seq12_0[EXTRA_PREBASE:]
            seq2_0 = seq12_0[len(self.ref)+EXTRA_PREBASE:]
            left, seq1, seq2 = self.rightAlign(seq1_0, seq2_0)

        base = ''
        if len(seq1) == 0 or len(seq2) == 0:  # e.g. if insertion/deletion .. should always be the case
            if (left >= 1 and EXTRA_PREBASE == 1):
                # base = reference.getReference(self.chrom, self.pos + left, self.pos + left)
                base = seq12_0[left]
            elif left == 0 and EXTRA_PREBASE == 1:
                base = seq12_0[0]
            else:  # Variant on First base of chrom (currently only Mitochorndria has first base not an 'N')
                base = ""
            seq1, seq2 = base + seq1, base + seq2
        ret = Variant(self.chrom, self.pos + left-len(base), seq1, seq2)  # have to use VCF padded coordinates
        ret.flags = self.flags
        self_flagvalues = self.flagvalues
        ret.flagvalues = self_flagvalues
        self.alignonplus = ret
        return ret

    # Aligning variant on the minus strand (left shift), so appropriate for variants on mimis strand.
    def alignOnMinusStrand(self, reference):
        if self.alignonminus is not None:
            return self.alignonminus
        if self.pos == 1 or (self.chrom not in reference.reflens or not (self.is_insertion or self.is_deletion)):
            [self.pos, self.ref, self.alt] = self.remove_same_bases_ref_alt()
            self.alignonminus = self
            return self
        else:
            if self.is_insertion:
                if not len(self.vcf_padded_base) == 0:
                    if self.vcf_padded_base[0] != self.alt[-1]:  # See if there is the potential for another copy of the alt allele
                        self.alignonminus = self
                        return self
            elif self.is_deletion:
                if not len(self.vcf_padded_base) == 0:  # Skip cases where initial VCF was invalid and did not have pre-base as can happen after splitting multiple alleles
                    if self.vcf_padded_base[0] != self.ref[-1]:
                        self.alignonminus = self
                        return self

        reflen = reference.reflens[self.chrom]

        maxreplen = abs(len(self.ref)-len(self.alt))
        PADLEN = max(100,1+5*maxreplen)
        if self.pos - PADLEN<1:
            PADLEN = self.pos -1
        seq1_0 = reference.getReference(self.chrom, self.pos - PADLEN, self.pos + len(self.ref) - 1)
        #s = reference.getReference(self.chrom, self.pos - PADLEN, self.pos - 1)
        s = seq1_0[0:PADLEN]
        seq2_0 = s + self.alt
        N = len(s)
        left, seq1, seq2 = self.leftAlign(seq1_0, seq2_0)
        while left <= maxreplen and self.pos-PADLEN > 0:  # Shifted  all the way to the beginning, make sure there is not more to shift
            PADLEN += PADLEN
            if self.pos - PADLEN < 1:  # This will only occur for Mitochondria
                PADLEN = self.pos-1
            seq1_0 = reference.getReference(self.chrom, self.pos - PADLEN, self.pos + len(self.ref) - 1)
            #s = reference.getReference(self.chrom, self.pos - PADLEN, self.pos - 1)
            s = seq1_0[0:PADLEN]
            seq2_0 = s + self.alt
            left, seq1, seq2 = self.leftAlign(seq1_0, seq2_0)
            N = len(s)

        if len(seq1) == 0 or len(seq2) == 0:
            # Unless Variant got shifted all the way to position 1(left==0), there should be enough bases for padding.
            if left > 0:
                left = left - 1
                #            base = reference.getReference(self.chrom, self.pos + left - PADLEN, self.pos + left - PADLEN)
                base = seq1_0[left]  # if eliminate 1st common base, take seq1_0[0] ..etc
            else:  # left == 0 means no shared bases in chunk, so it means variant fills the chunk .. all the way to the first base
                base = ''
            seq1, seq2 = base + seq1, base + seq2
        else:
            sys.stderr.write(
                "WARNING: Attempted to left shift non-indel variant .. bug? for variant:" + self.id + "\n")


        ret = Variant(self.chrom, self.pos + left - N, seq1, seq2)
        ret.flags = self.flags
        ret.flagvalues = self.flagvalues
        self.alignonminus = ret
        return ret

    # Right-aligning two sequences
    def rightAlign(self, seq1, seq2):
        left, seq1, seq2 = self.trimCommonStart(seq1, seq2)
        right, seq1, seq2 = self.trimCommonEnd(seq1, seq2)
        return left, seq1, seq2

    # Left-aligning two sequences
    def leftAlign(self, seq1, seq2):
        right, seq1, seq2 = self.trimCommonEnd(seq1, seq2)
        left, seq1, seq2 = self.trimCommonStart(seq1, seq2)
        return left, seq1, seq2

    # Trimming common starting subsequence of two sequences
    def trimCommonStart(self, s1, s2):
        counter = 0
        while True:
            if len(s1)<= counter or len(s2) <=counter or s1[counter] != s2[counter]:
                return counter, s1[counter:], s2[counter:]
            #s1, s2 = s1[1:], s2[1:]
            counter += 1

    # Trimming common ending subsequence of two sequences
    def trimCommonEnd(self, s1, s2):
        counter = 1
        if len(s1) == 0 or len(s2) ==0 or (s1[-1]!=s2[-1]):
            return 0, s1, s2
        while True:
            if counter > len(s1) or counter > len(s2):  # implicitely counter != 1
                return counter-1, s1[:-(counter-1)], s2[:-(counter-1)]
            if s1[-counter] != s2[-counter]:
                    return counter-1, s1[:-(counter-1)], s2[:-(counter-1)]
            counter += 1



#######################################################################################################################
def trim_all_alleles(ref,alts): # This is a VCF-style trimming, must leave 1 base
    ntrim = 0
    if len(ref)==1:
        return 0
    for ir in range(0, len(ref)-1):
        refbase = ref[ir]
        for ia in range(0, len(alts)): # loop over alleles
            if ntrim >= len(alts[ia])-1:  # Don't trim to the point of leaving empty alleles
                return ntrim
            if refbase != alts[ia][ir]:
                return ntrim
        ntrim = ntrim +1
    return ntrim



# Class representing a single VCF record     
class Record(object):
    # Constructor
    def __init__(self, line, options, targetBED, reference):

        self.targetBED = targetBED
        # Reference genome pysam/tabix object
        self.reference = reference
        # Parsing VCF format
        if 'build' in options.args:
            self.build = options.args['build']
        else:
            self.build = 'GRCh38'
        if options.args['inputformat'].upper() == 'VCF':
            cols = line.strip().split("\t")
            self.chrom = cols[0]

            if self.chrom.startswith('chr'):
                self.chrom_chr_prefix = True
                self.chrom = self.chrom[3:]
            else:
                self.chrom_chr_prefix = False
            try:
                self.pos = int(cols[1])
            except:
                sys.stderr.write("Init VCF: invalid position to line:"+line+"\n")
                sys.exit(1)
            self.id = cols[2]
            self.ref = cols[3]
            alts = cols[4].split(",")
            self.qual = cols[5]
            self.filter = cols[6]
            if len(cols) > 7:
                self.info = cols[7]
            else:
                self.info = ''
            if len(cols) > 8:
                self.rest = cols[8:]
            else:
                self.rest = []
        elif options.args['inputformat'].upper() == 'TXT':  # Parsing TXT format
            cols = line.strip().split("\t")
            self.id = cols[0]
            self.chrom = cols[1]
            self.pos = int(cols[2])
            self.ref = cols[3]
            alts = cols[4].split(",")

            self.qual = ''
            self.filter = 'PASS'
            self.info = ''
            self.rest = []
        else:
            raise Exception("Invalid Input format type")

        altsclean = []
        for alt in alts:
            if not alt.startswith('<'):
                altsclean.append(alt)
        # Creating a Variant object for each variant/alt-allele  in the record
        # as long as they pass filtering.
        #
        # Trim only if all alleles have an extra base.. but this really should be done by a normalizing code
        # And one should really remove all ale alleles.
        self.variants = []
        self.alts = []
        in_target_region = False

        ntrim = trim_all_alleles(self.ref,altsclean)
        if ntrim > 0:
            self.ref = self.ref[ntrim:]
            self.pos += ntrim
            for ia in range(0,len(alts)):
                alt = alts[ia]
                if not alt.startswith('<'):
                    alts[ia] = alts[ia][ntrim:]

        for alt in alts:
            # Keep original alt alleles because each alt-allele may normalize differently and have different positions
            self.alts.append(alt)
            # Initializing each Variant object with different alt allele
            if alt.startswith('<'):
                var = None
            else:
                var = Variant(self.chrom, self.pos, self.ref, alt)

                # DO  NOT FILTER
                """
                if 'N' in self.ref or 'N' in alt:
                    logging.info('Variant ignored as allele contains unknown base (\'N\'): ' + self.chrom + ':' + str(
                        self.pos) + ' ' + self.ref + '>' + alt)
                    continue
               
                if alt == '.':
                    logging.info("Variant ignored because it is monomorphic reference: " + self.chrom + ':' + str(
                        self.pos) + ' ' + self.ref + '>' + alt)
                    continue
    
                if not alt.strip('ACGT') is '' or not self.ref.strip('ACGT') is '':
                    logging.info("Variant ignored as format of alt allele is not supported: " + self.chrom + ':' + str(
                        self.pos) + ' ' + self.ref + '>' + alt)
                    continue
                """
                # Filtering by variant type (i.e. substitution, indel, insertion, deletion, complex indel), if required
                # Filtering via those options is NOT recommended, especially when dealing with multiple alt-alleles.

                if options.args['type'].upper() == 'SUBSTITUTION' and not var.is_substitution:
                    continue
                if options.args['type'].upper() == 'INDEL' and not var.is_insertion and not var.is_deletion \
                        and not var.is_complex:
                    continue
                if options.args['type'].upper() == 'INSERTION' and not var.is_insertion:
                    continue
                if options.args['type'].upper() == 'DELETION' and not var.is_deletion:
                    continue
                if options.args['type'].upper() == 'COMPLEX' and not var.is_complex:
                    continue
                # Adding Variant object to this record
                # 10/2021-New benavior, adding all alt-alleles variants if any of the alleles is in the target region .. and is not explicitely filtered out by type
                # The logic is that the user should be able to decide which variant has the best support in the sequencing data (read count)
                # .. and not be biased by the most deleterious option.
                # At least 1 must pass the BED filtering - if requested (otherwise all alt-alleles (whole variant) are excluded)
            self.variants.append(var)

            # Filtering by BED file, if required
            if targetBED is not None and var is not None:
                goodchrom = convert_chrom(var.chrom+"", targetBED.contigs)
                if goodchrom is None:
                    continue
                if not var.is_insertion:
                    start = var.pos - 1  # correct for 0-base
                    end = var.pos + len(var.ref) - 1
                else:
                    start = var.pos - 2 # correct for insertion AND 0-base
                    end = var.pos
                foundstart = False
                for _ in self.targetBED.fetch(region=goodchrom + ':' + str(start) + '-' + str(start)):
                    foundstart = True
                foundend = False
                for _ in self.targetBED.fetch(region=goodchrom + ':' + str(end) + '-' + str(end)):
                    foundend = True
                if not (foundstart or foundend):
                    if not (var.is_deletion or var.is_insertion):
                        continue  # No shifting possible, so definitively outside bed region.
                    else:
                        var_plus = var.alignOnPlusStrand(self.reference)
                        var_minus = var.alignOnMinusStrand(self.reference)
                        if var_plus.pos == var_minus.pos:  # No change, so no shifting
                            continue
                        else:
                            if not var_plus.is_insertion:
                                start = var_plus.pos -1 # correct for -base
                                end = var_plus.pos + len(var_plus.ref) - 1
                            else:
                                start = var_plus.pos - 2 # correct for insertion and 0-base
                                end = var_plus.pos
                            foundstart = False
                            for _ in self.targetBED.fetch(region=goodchrom + ':' + str(start) + '-' + str(start)):
                                foundstart = True
                            foundend = False
                            for _ in self.targetBED.fetch(region=goodchrom + ':' + str(end) + '-' + str(end)):
                                foundend = True
                            if not (foundstart or foundend):
                                if not var_minus.is_insertion:
                                    start = var_minus.pos -1  # correct for 0-base
                                    end = var_minus.pos + len(var_minus.ref) - 1
                                else:
                                    start = var_minus.pos - 2 # correct for 0-base and insertion
                                    end = var_minus.pos
                                foundstart = False
                                for _ in self.targetBED.fetch(region=goodchrom + ':' + str(start) + '-' + str(start)): foundstart = True
                                foundend = False
                                for _ in self.targetBED.fetch(region=goodchrom + ':' + str(end) + '-' + str(end)): foundend = True
                                if not (foundstart or foundend):
                                    continue
                in_target_region = True
            else:
                in_target_region = True

        if in_target_region is False:
            self.variants = []

    # Annotating record
    def annotate(self, ensembl, dbsnp, reference, impactdir):
        # Annotating each variant in the record
        for variant in self.variants:
            if variant is not None:
                variant.annotate(ensembl, dbsnp, reference, impactdir)
                if variant.alignonplus is not None:
                    variant.addFlag('HGVSg', csn.get_genomic_Annotation(variant.alignonplus, self.build ,reference))
                else:
                    variant.addFlag('HGVSg', csn.get_genomic_Annotation(variant, self.build, reference))

    # Writing record (a ref, multiple alts) to output file
    def output(self, outformat, outfile, options, genelist, transcriptlist, snplist, stdout):
        outvariants = []
        outalts = []
        build = "GRCh38"
        if "build" in options.args:
            build = options.args["build"]

        # Values of Entries for Separate Alleles are separated by ','
        # Entries for Multiple transcripts (for 1 allele) are separated by ':'
        #  This allows each allele to have different transcripts
        #  (for multiple allele, shifting may lead to different alleles in shifted vs non-shifted alt-alleles)

        for i in range(len(self.variants)):
            variant = self.variants[i]
            if variant is None:
                outalts.append(self.alts[i])
            else:
                isTRANSCRIPT = ('TRANSCRIPT' in variant.flags)
                if isTRANSCRIPT:
                    annTRANSCRIPT = variant.flagvalues[variant.flags.index('TRANSCRIPT')]
                else:
                    annTRANSCRIPT = ''
                isDBSNP = ('DBSNP' in variant.flags)
                if isDBSNP:
                    annDBSNP = variant.flagvalues[variant.flags.index('DBSNP')]
                else:
                    annDBSNP = ''

                if len(genelist) > 0 or len(transcriptlist) > 0:
                    if annTRANSCRIPT == '':
                        continue

                # Removing non-annotated variants/alleles, if required
                if 'nonannot' in options.args and options.args['nonannot'] is False and \
                        not ((isTRANSCRIPT and not annTRANSCRIPT == '') or (isDBSNP and not annDBSNP == '')):
                    continue

                # Filtering by gene, transcript or snp list, if required
                if isDBSNP:
                    if len(snplist) > 0 and annDBSNP not in snplist: continue


                outvariants.append(variant)
                outalts.append(self.alts[i])  # Original Alts.

        # Skipping record if all alt-allelic variants have been removed
        if len(outvariants) == 0 or len(outalts)==1 and outalts[0].startswith("<"):
            return

        # Writing output in VCF format
        if outformat.upper() == 'VCF':
            infos = self.info.split(';')
            newinfos = []
            rmfields = ['HGVSg','CAVA_HGVSg','HGVSp','CAVA_HGVSp','HGVSc','CAVA_HGVSc',
                'TYPE', 'TRANSCRIPT', 'GENE', 'GENEID', 'TRINFO', 'LOC', 'CSN',
                        'PROTPOS', 'PROTREF', 'PROTALT', 'CLASS', 'SO', 'ALTFLAG',
                        'CAVA_TYPE', 'CAVA_TRANSCRIPT', 'CAVA_GENE', 'CAVA_GENEID', 'CAVA_TRINFO',
                        'CAVA_LOC', 'CAVA_CSN',
                        'CAVA_PROTPOS', 'CAVA_PROTREF', 'CAVA_PROTALT', 'CAVA_CLASS', 'CAVA_SO', 'CAVA_ALTFLAG']
            for item in infos:
                items = item.split("=")
                if not ( items[0] in rmfields):
                    newinfos.append(item)
            self.info = ';'.join(newinfos)


            # Creating first part of the VCF record (up to FILTER field)

            chromstr = 'chr' + self.chrom if self.chrom_chr_prefix else self.chrom
            record = [chromstr, str(self.pos), self.id, self.ref, ",".join(outalts), self.qual, self.filter]

            # Preparing components of the String to be added to the INFO field
            flags = []
            flagvalues = []
            for variant in outvariants:
                for i in range(len(variant.flags)):
                    key = variant.flags[i]
                    value = variant.flagvalues[i]
                    if value == '':
                        value = '.'
                    if key in flags:
                        flagvalues[flags.index(key)].append(value)
                    else:
                        flags.append(key)
                        flagvalues.append([value])

            # Creating String added to the INFO field, separating multiple alleles by ',' (within alt alleles, multiple
            # transcripts separated by ':')
            added = ''
            for i in range(len(flags)):
                key = flags[i]
                value = ','.join(flagvalues[i])
                if len(added) > 0:
                    added += ';'
                if 'prefix' in options.args and options.args['prefix']:
                    added += 'CAVA_' + key + '=' + value
                else:
                    added += key + '=' + value

            ########################################################################
            # Add HGVS
            # HS: 10/21/2020 Added HGVSP and multi-transcript support to HGVSc
            #                added prefix support

            if 'prefix' in options.args and options.args['prefix'] is True:
                HGVSC_key = 'CAVA_HGVSc='
                HGVSP_key = 'CAVA_HGVSp='
 #               HGVSG_key = 'CAVA_HGVSg='

            else:
                HGVSC_key = 'HGVSc='
                HGVSP_key = 'HGVSp='
 #               HGVSG_key = 'HGVSg='

            HGVSC = ''
            HGVSP = ''
            # hgvsgs = flagvalues[flags.index('HGVSg')]
            # HGVSG = ''
            # if hgvsgs is not None:
            #     HGVSG = hgvsgs[0]
            contig = csn.get_contig_from_build(outvariants[0].chrom,build)
            if ('TRANSCRIPT' in flags) and ('GENE' in flags) and ('CSN' in flags):
                # List of 1 entry per variant (alt-allele) .. and within each multiple transcripts are separated by ','
                transcripts_list = flagvalues[flags.index('TRANSCRIPT')]
                gene_list = flagvalues[flags.index('GENE')]
                csn_list = flagvalues[flags.index('CSN')]
                #repeats_list = flagvalues[flags.index('REPEATS')]
                csn_hgvs_list = csn_list
                #csn_hgvs_list = flagvalues[flags.index('CSNHGVS')]
                if not (len(transcripts_list) == 0 or len(gene_list) != len(transcripts_list) or len(gene_list) != len(
                        csn_list)):
                    for ihg in range(0, len(transcripts_list)):  # Loop over alt-alleles
                        hgtranscripts = transcripts_list[ihg].split(":")  # split over alt-alleles
                        hggenes = gene_list[ihg].split(":")
                        hgcsns = csn_list[ihg].split(":")
                        hgcsns_hgvs = csn_hgvs_list[ihg].split(":")
                        #hgrepeats = repeats_list[ihg].split(":")

                        if not (len(hgtranscripts) == len(hggenes) and len(hgtranscripts) == len(hgcsns)):
                            print("ERROR transcripts GENE and CSN not same length\n")
                            print("--------------------------------------------------------------------\n")
                            if options.args['logfile']:
                                logging.error("Bug: TRANSCRIPT GENE and CSN not same length\n")
                            sys.exit(1)
                        else:
                            trHGVSC = ''
                            trHGVSP = ''
                            for itr in range(0, len(hgtranscripts)):  # Loop over transcripts within each alt-alleles
                                hgtranscript = hgtranscripts[itr]
                                hggene = hggenes[itr]
                                hgcsn_hgvs = hgcsns_hgvs[itr]
                                #hgrepeat = hgrepeats[itr]
                                if len(hgtranscript)==0 or hgtranscript == '.':
                                    tHGVSC = '.(.):'
                                else:
                                    tHGVSC = contig + '(' + hgtranscript + '):'
                                try:
                                    cdna, prot = hgcsn_hgvs.split('_p.')
                                except ValueError:  # Example c.802-51_802-14del38, splice
                                    cdna = hgcsn_hgvs
                                    prot = '.'
                                # HGVS does not allow nucs after del c.21_22delAA
                                del_nuc = re.match(r"^(.*del[ACGTacgtnN]+)$", cdna)
                                if del_nuc:
                                    cdna = del_nuc.group(1)
                                #if len(hgrepeat)>1:
                                #    tHGVSC += hgrepeat
                                #else:
                                #   tHGVSC += cdna
                                tHGVSC += cdna
                                if cdna == '.' or  len(cdna)==0 or tHGVSC == '.(.):.' or tHGVSC == '.(.):':
                                    tHGVSC = '.'
                                if prot == '.' or prot == '':
                                    tHGVSP = '.'
                                else:
                                    # Get Protein matching Transcript from options
                                    if hgtranscript in options.transcript2protein:
                                        protid = options.transcript2protein[hgtranscript]
                                        tHGVSP = protid + ':p.(' + prot + ')'
                                    else:
                                        # tHGVSP=hgtranscript+':p.('+prot+')'
                                        # If user wants HGVSP, they can always look in the CSN... but don't want to
                                        # provide incorrect HGVSP
                                        tHGVSP = '.'
                                        # Only report warning if transcript2protein mapping file was provided.
                                        if len(options.transcript2protein) > 0:
                                            print(
                                                "WARNING: transcript " + hgtranscript + " not in transcript2protein "
                                                                                        "file, HGVSp will be invalid\n")
                                            print(
                                                "--------------------------------------------------------------------"
                                                "\n")
                                            if options.args['logfile']:
                                                logging.info(
                                                    "WARNING: transcript " + hgtranscript + " not in transcript2protein"
                                                                                            " file\n")
                                if trHGVSC == '':
                                    trHGVSC = tHGVSC
                                    trHGVSP = tHGVSP
                                else:
                                    trHGVSC = trHGVSC + ':' + tHGVSC
                                    trHGVSP = trHGVSP + ':' + tHGVSP
                            if HGVSC == '':
                                HGVSC = trHGVSC
                                HGVSP = trHGVSP
                            else:
                                HGVSC = HGVSC + ',' + trHGVSC
                                HGVSP = HGVSP + ',' + trHGVSP

            if HGVSC == '':
                HGVSC = HGVSC_key + 'None'
                HGVSP = HGVSP_key + 'None'
            else:
                HGVSC = HGVSC_key + HGVSC
                HGVSP = HGVSP_key + HGVSP

            # Add multi-transcripts/multi-allele HGVS to output record
            if added == '':
                added = HGVSC + ';' + HGVSP
            else:
                added += ';' + HGVSC + ';' + HGVSP

            # Adding second part of the VCF record (starting from the INFO field)
            if self.info == '.' or self.info == '':
                record += [added]
                record += self.rest
            else:
                record += [self.info + ';' + added]
                record += self.rest

            # Writing record to the output file
            if stdout:
                print('\t'.join(record))
            else:
                outfile.write('\t'.join(record) + '\n')

        # Writing output in TSV format
        # 10/22/2020 added HGVS and HGVSp between core record and extras.

        if outformat.upper() == 'TSV':

            # Iterating through variants 
            #   Print one line for each alt-allele*Variant Combination
            c = 0
            flags = []
            flagvalues = []
            for variant in outvariants:
                for i in range(len(variant.flags)):
                    key = variant.flags[i]
                    value = variant.flagvalues[i]
                    if value == '': value = '.'
                    if key in flags:
                        flagvalues[flags.index(key)].append(value)
                    else:
                        flags.append(key)
                        flagvalues.append([value])

                # Standardize chromosome M notation if options specify
                if options.args['normalized_mitochondrial_chrom'] == 'MT':
                    if self.chrom == 'chrM':
                        self.chrom = 'chrMT'
                    if self.chrom == 'M':
                        self.chrom = 'MT'
                if options.args['normalized_mitochondrial_chrom'] == 'M':
                    if self.chrom == 'chrMT':
                        self.chrom = 'chrM'
                    if self.chrom == 'MT':
                        self.chrom = 'M'

                # Creating first part of the TSV record (up to FILTER field)
                # Common to all transcripts for that variant
                record = self.id + '\t' + self.chrom + '\t' + str(self.pos) + '\t' + self.ref + '\t' + outalts[
                    c] + '\t' + self.qual + '\t' + self.filter
                contig = csn.get_contig_from_build(self.chrom, build)
                # Number of transcripts overlapping with the variant
                if 'TRANSCRIPT' in variant.flags:
                    transcripts_list = variant.flagvalues[variant.flags.index('TRANSCRIPT')].split(":")
                    N = len(transcripts_list)
                    if N > 0 and ('GENE' in variant.flags) and ('CSN' in variant.flags):
                        gene_list = variant.flagvalues[flags.index('GENE')].split(":")
                        csn_list = variant.flagvalues[flags.index('CSN')].split(":")
                    else:
                        if len(transcripts_list) > 0:
                            print("ERROR: transcript " + transcripts_list[0] + " has no matching GENE and/or CSN\n")
                            print("--------------------------------------------------------------------\n")
                            if options.args['logfile']:
                                logging.info(
                                    "ERROR: transcript " + transcripts_list[0] + " has no matching GENE and/or CSN\n")
                            quit()
                        gene_list = []
                        csn_list = []
                        transcripts_list = []
                else:  # Allow having No matching transcripts (for variant outside genes/transcripts)
                    N = 1
                    transcripts_list = []
                    gene_list = []
                    csn_list = []

                # Iterating through the transcripts
                for i in range(N):
                    # Creating second part of the TSV record - Transcript-specific per line.
                    rest = ''
                    for j in range(len(variant.flags)):
                        # First, do not split annotation that is not transcript specific
                        if not variant.flags[j] in ['TRANSCRIPT', 'GENE', 'GENEID', 'TRINFO', 'LOC', 'CSN', 'CLASS',
                                                    'SO',
                                                    'IMPACT', 'ALTANN', 'ALTCLASS', 'ALTSO', 'ALTFLAG', 'PROTPOS',
                                                    'PROTREF', 'PROTALT']:
                            value = variant.flagvalues[j]
                            if value == '':
                                value = '.'
                            rest += '\t' + value
                            continue

                        values = variant.flagvalues[j].split(':')
                        value = values[i]
                        if value == '': value = '.'
                        rest += '\t' + value
                    if len(transcripts_list) == 0 or len(gene_list) != len(transcripts_list) or len(gene_list) != len(
                            csn_list) or len(transcripts_list) != N:
                        HGVSC = 'None'
                        HGVSP = 'None'
                    else:
                        hgtranscript = transcripts_list[i]
                        if len(hgtranscript) == 0 or hgtranscript == '.':
                            HGVSC = '.(.):'
                        else:
                            HGVSC = contig + '(' + hgtranscript + '):'
                        try:
                            cdna, prot = csn_list[i].split('_p.')
                            prot = prot.replace('X', "Ter")
                        except ValueError:  # Example c.802-51_802-14del38, splice
                            cdna = csn_list[i]
                            prot = '.'
                        HGVSC += cdna
                        if len(cdna)==0 or cdna == '.' or HGVSC == '.(.):.' or HGVSC == '():' or HGVSC == '.(.):':
                            HGVSC = '.'
                        if prot == '.' or prot == '':
                            HGVSP = '.'
                        else:
                            # Get Protein matching Transcript from options
                            if hgtranscript in options.transcript2protein:
                                protid = options.transcript2protein[hgtranscript]
                                HGVSP = protid + ':p.(' + prot + ')'
                            else:
                                #                                HGVSP=hgtranscript+':p.('+prot+')'
                                HGVSP = '.'  # If user wants HGVSP, they can always look in the CSN... but don't want to provide incorrect HGVSP
                                if len(
                                        options.transcript2protein) > 0:  # Only report warning if transcript2protein mapping file was provided.
                                    print(
                                        "WARNING: transcript " + hgtranscript + " not in transcript2protein file, HGVSp will be invalid\n")
                                    print("--------------------------------------------------------------------\n")
                                    if options.args['logfile']:
                                        logging.info(
                                            "WARNING: transcript " + hgtranscript + " not in transcript2protein file\n")
                    # Writing record to the output file
                    if stdout:
                        print(record + rest + "\t" + HGVSC + "\t" + HGVSP )
                    else:
                        outfile.write(record + rest + "\t" + HGVSC + "\t" + HGVSP + '\n')

                c += 1


#######################################################################################################################
class Tr_store:
    lasttr = dict
    
# Class representing a single Ensembl transcript
# noinspection PyUnresolvedReferences
class Transcript(object):
    # Constructor
    def __init__(self, line):
        self.exons = []
        cols = line.split('\t')
        self.TRANSCRIPT = cols[0]
        self.geneSymbol = cols[1]
        self.geneID = cols[2]
        self.TRINFO = cols[3]
        self.chrom = cols[4]
        self.strand = int(cols[5])
        self.transcriptStart = int(cols[6])  # 0-based, lowest coordinate if first exon
        self.transcriptEnd = int(cols[7])   # 1-bases upper coordinate
        self.codingStart = int(cols[8])  # in cDNA coordinated (1st base is 1)
        self.codingStartGenomic = int(cols[9])# coding start genomic 1-based
        self.codingEndGenomic = int(cols[10])# coding end genomic 1-based (includes stop codon)
        self.is_selenocysteine = False
        # Initializing and adding exons
        for i in range(1, len(cols) - 11, 2):
            self.exons.append(Exon(int((i + 1) / 2), int(cols[10 + i]), int(cols[11 + i])))   # Start is 0-based, ENd is 1-based -- like a bed file
        self.three_prime_len = 0
        self.exonseqs = None  # list of the exons for the reference sequence
        self.cds_len = 0
        foundStart = False
        for ex in self.exons:
            if self.strand == 1:
                if self.codingStartGenomic > ex.start and self.codingStartGenomic<=ex.end: # CDS starts within exon
                    if self.codingEndGenomic<=ex.end: # Both start and end within exon
                        self.cds_len = self.codingEndGenomic+1-self.codingStartGenomic
                        if self.codingEndGenomic<ex.end:  # is there anything left for 3' utr in this exon
                            self.three_prime_len += ex.end-self.codingEndGenomic  #ex.end+1-(self.codingEndGenomic+1)
                    else: # Start in this one, but ends in another exon
                        self.cds_len = ex.end +1 - self.codingStartGenomic
                elif self.codingStartGenomic <=ex.start: # Already Past start codon
                    if self.codingEndGenomic<=ex.end and self.codingEndGenomic> ex.start: # cds ends inside this exon
                        self.cds_len += (self.codingEndGenomic-ex.start)
                        if self.codingEndGenomic<ex.end:  # is there anything left for 3' utr in this exon
                            self.three_prime_len += ex.end-self.codingEndGenomic  #ex.end+1-(self.codingEndGenomic+1)
                    elif ex.end <= self.codingEndGenomic: # End in downstream exon
                        self.cds_len += ex.end-ex.start
                    else: # Already past cds
                        self.three_prime_len += ex.end-ex.start
            else:  # Minus strand (exons are in transcription order)
                if self.codingStartGenomic > ex.start and self.codingStartGenomic<=ex.end: # start within this exon
                    if self.codingEndGenomic>ex.start: # Both start and end within exon
                        self.cds_len = self.codingStartGenomic +1 - self.codingEndGenomic
                        if ex.start+1<self.codingEndGenomic:  # is there anything left for 3' utr in this exon
                            self.three_prime_len += (self.codingEndGenomic-1 - ex.start)
                    else: # end in another exon
                        self.cds_len = self.codingStartGenomic - ex.start
                elif self.codingStartGenomic > ex.end: # Already passed start codon
                    if self.codingEndGenomic < ex.start:  # End in downstream exon or includes whole exon
                        self.cds_len += ex.end - ex.start
                    elif self.codingEndGenomic<=ex.end and self.codingEndGenomic> ex.start: # end inside this exon
                        self.cds_len += (ex.end +1 - self.codingEndGenomic ) # already  covered case where start and stop both in same exon.
                        if ex.start + 1 < self.codingEndGenomic:  # is there anything left for 3' utr in this exon
                            self.three_prime_len += (self.codingEndGenomic - 1 - ex.start)
                    else: # Already past cds
                        self.three_prime_len += ex.end-ex.start

    def trimSelenoCysteine(self, ref_cds, mut_cds, protein, mutprotein, variant):
        # The protein and murprotein is translated (without truncation) to UGT.
        #     The protein is truncated according to CDS length.
        # Assume this was only called for a Selenocysteine gene.
        # Selenocystein genes have an element in the 3'UTR that directs translation of UGA (TGA) into selenocystein.
        # The location of this Selenocysteine Insertion Sequence (SECIS) is annotated on genbank records,
        #  Too close to the SECIS and it stops working, but there is no universal distance cutoff
        # The difficulty is to know when when (along transcript) this mechanisms stops to work
        #  .. the only thing we know for sure is that it stops working at the site of the stop codon (if there is
        #      a TGA as the stop codon.
        # We know that the location of the first to the last UGT.
        # We will assume that UGT translation works from the first base of the CDS>
        # a) Track the location of the original Stop Codon along the position of the mutated protein,
        #    translate the mutated protein with a Selenocysteine genetic code.
        if protein is None:
            return mutprotein
        allUindexes = [_x.start() for _x in re.finditer("U",protein)]
        if len(allUindexes) == 0:  # No U in translated protein==>not a selenocysteine gene
            sys.stderr.write(
                "WARNING: gene " + self.geneSymbol + " configured as selenocysteine gene, but no U found in translation of reference CDS\n")
            return mutprotein
        if mutprotein is None:
            return mutprotein
        allmutUindexes = [_x.start() for _x in re.finditer("U",mutprotein)]
        if len(allmutUindexes) == 0:  # No issue with U being possibly a stop codon
            return mutprotein
        # The reference protein should always be Stop terminated
        # For Selenocysteine, this could be a "U" ..
        ref_stop_AAindex = len(protein)-1
        ref_stop_codon = ref_cds[(3*ref_stop_AAindex):(3*ref_stop_AAindex+3)]
        AA = Sequence(ref_stop_codon).translate('7')
    # Compute the length change of the protein from an insertion or deletion
    # Even if it's a 1 bp frameshift, count it as 1 AA
        if len(variant.ref)==len(variant.alt):
            variant_AA_size = 0
        elif len(variant.ref)> len(variant.alt):  # deletion
            variant_AA_size = int(abs((len(variant.alt)-len(variant.ref)) +2)/3)
        else: #  insertion
            variant_AA_size = - int(abs((len(variant.ref)-len(variant.alt)) +2)/3)
        # Assume that variant is all inside CDS to compute the shifted_mut_stop_pos
        #  therefore Variant cannot be past the reference stop codon

        mut_stop_AAindex = ref_stop_AAindex + variant_AA_size
        LASTU = allUindexes[len(allUindexes) - 1]
        mut_stopX = mutprotein.find("X")
        if mut_stopX == -1:
            mut_stopX = len(mutprotein) -1
        # Find first AA that is different.
        isame = -1
        while len(protein) > isame + 1 and len(mutprotein) > isame + 1:
            if protein[isame + 1] == mutprotein[isame + 1]:
                isame = isame + 1
            else:
                break
        if isame == len(protein) -1: # Same up to and including Stop codon
            return mutprotein[0:len(protein)]
        different_index = isame + 1  # first AA position after AA are the same (if isame==len(protein
# if variant is before LASTU, shift LASTU position to adjust codinates
        if different_index < LASTU: # AA position (in mut coordinates) at which we know for sure CESIS works
            LASTU += variant_AA_size
            if LASTU<different_index: # Cannot shift lastU before Variant
                LASTU=different_index
        if ref_stop_codon == 'TGA': # Selenocysteine codon is not effective past this point.
            # If the Stop codon is TGA, there is a region
            # betwem the last U and the stop codon, where we cannot predict what a TGA will do
            for mutindex in allmutUindexes:
                if mut_stopX<mutindex:
                    return mutprotein[0:(mut_stopX+1)]
                if LASTU < mutindex < mut_stop_AAindex: # Region where we do not knowwether TGA acts
                    variant.seleno_early_stop = mutindex+1 # THis is a candidate for an early TGA stop codon, in AA positions
                    return mutprotein
            for mutindex in allmutUindexes:
                if mutindex > mut_stop_AAindex: # TGA Should be working as a Stop codon now after this position
                    return mutprotein[0:(mutindex+1)]
            return mutprotein[0:(mut_stopX+1)]
        else: # Reference Stop Codon is regular stop codon, so we don't know how far past position of last wildtype U TGA gets recoded Ter to vs U
            for mutindex in allmutUindexes:
                if mut_stopX<mutindex: # Early regular Stop wins
                    return mutprotein[0:(mut_stopX+1)]
                if mutindex>LASTU: # past point where we know what happens with Selenocystein stops
                    # Must be new TGA stop codon.
                    variant.seleno_early_stop = mutindex+1 # THis is a candidate for an early stop codon, in AA positions
                    return mutprotein
            return mutprotein[0:(mut_stopX + 1)]


# Get coding sequence of the Alternate Transcript.
# Don't try to detect effect of variants affecting the splice site, user us expected to check if overlap EE
    # raise exception if variant would change ACROSS intron/exon junction.

    def getCodingSequence(self, reference, variant, exonseqs):
        ret = ''  # transcript sequence including alt-allele
        ret_ref = ''  # transcript sequence with ref allele
        ret_exonseqs = []
        coding_start_alt = self.codingStart
        trans_pos = 0
        if exonseqs is None and self.exonseqs is not None:
            exonseqs = self.exonseqs

        for i in range(len(self.exons)):
            # Exons are ordered 5' to 3; (different order for strand = 1 or -1)
            exon = self.exons[i]
            if exonseqs is None:
                # bef = reference.getReference(self.chrom, exon.start + 1, variant.pos - 1)  # 1-base positions
                # aft = reference.getReference(self.chrom, variant.pos + len(variant.ref), exon.end)  # 1-base positions
                exonseqF = reference.getReference(self.chrom, exon.start + 1, exon.end)
                if self.strand == 1:
                    exonseq = exonseqF
                else:
                    exonseq = exonseqF.reverseComplement()
            else:
                exonseq = exonseqs[i]

            if variant is not None:
                # start of the exon is "0"-based (end is 1-based)
                # The position is not the purelyVCF position. The MNP/SNP indel position is 1-based and point to the affected based
                # but the insertion position, points 1-bp AFTER the variant (variant.ref="")
                if self.strand == 1:
                    if variant.is_insertion is False: # Del, MNP, SNP .. "regular" meaning of position
                        if exon.start < variant.pos <= exon.end:  #  same as exon.start+1 <= variant.pos <= exon.end
                            if  variant.pos + len(variant.ref) -1 > exon.end:  # variant starts in exon, but ends past exon
                                if i != len(self.exons) -1:
                                    raise Exception("Variant crossing intron/exon boundary.. Catch this exception and set mutprot=None")

                            trans_pos += (variant.pos - exon.start) # Position of variant in transcript coordinates
                            bef_end_index = variant.pos - exon.start -1 # Points to variant pos .. this index will not be included in sequence
                            upstream = exonseq[0:bef_end_index]  # note [0:0] == '', so if variant.pos is first base of exon, bef ==''
                            aft_start_index = variant.pos+len(variant.ref)- (exon.start+1)  # index points after variant, or other side of insertion.
                            downstream = exonseq[aft_start_index:]   # if variant.ref extend outside exon, will only grab sequence inside exon .. but realistically
                                                                # the splicing won't occur.. this is caught earlier in CAVA by checking type of variant
                                                                # and not protein annotating variants that span Extron-intron or Intron-Exon
                            temp_alt = upstream + variant.alt + downstream
                            if variant.pos+len(variant.ref)-1>exon.end: # Variant extending past last exon
                                temp_ref = upstream + variant.ref[0:exon.end + 1 - variant.pos]
                            else:
                                temp_ref = upstream + variant.ref+ downstream
                            ret_ref += temp_ref
                            ret += temp_alt
                            if exonseqs is None:
                                ret_exonseqs.append(temp_ref)
                            if trans_pos<coding_start_alt:
                                coding_start_alt += (len(variant.alt) - len(variant.ref))
                            continue
                        elif exon.start < variant.pos+ len(variant.ref)<= exon.end:  # variant starts before exon and ends in exon
                            raise Exception("Variant crossing intron/exon boundary.. Catch this exception and set mutprot=None")
                        elif variant.pos <exon.start and variant.pos+len(variant.ref)> exon.end: # variant encompasses whole exon
                            raise Exception("Variant changing whole exon .. Catch this exception and set mutprot=None")
                        else: #
                            pass  # Variant not affecting this exon

                    else: # pure insertion (ref len ==0). pos is base AFTER insertion site
                        if exon.start+1 <= variant.pos <= exon.end+1:  # Insertion allowed if touching the edges of the exon.II^EEEEII or OOEEEE^III OK
                            trans_pos += (variant.pos-1 - exon.start) # Position of variant in transcript coordinates
                            bef_end_index = variant.pos - exon.start -1 # 0-based index, Points to variant pos .. this index will not be included in sequence
                            upstream = exonseq[0:bef_end_index]  # note [0:0] == '', so if variant.pos is first base of exon, bef ==''
                            aft_start_index = variant.pos -(exon.start+1)  # 0-based index points after variant, or other side of insertion.
                            downstream = exonseq[aft_start_index:]   # if variant.ref extend outside exon, will only grab sequence inside exon .. but realistically
                                                                # the splicing won't occur.. this is caught earlier in CAVA by checking type of variant
                                                                # and not protein annotating variants that span Extron-intron or Intron-Exon
                                                                # if aft_start_index is past exon end .. downstream = ''
                            temp_alt = upstream + variant.alt + downstream
                            temp_ref = upstream + downstream
                            ret_ref += temp_ref
                            ret += temp_alt
                            if exonseqs is None:
                                ret_exonseqs.append(temp_ref)
                            if trans_pos<coding_start_alt:  # If inserted before START codon, must shift offset of coding codon.
                                coding_start_alt += len(variant.alt)
                            continue
                        else: #
                            pass  # Variant not affecting this exon
                else:  # Minus strand (exons are listed in coding order (reverse of genomic order)
                    if variant.is_insertion is False: # MNP, SNP, DEL .. "normal" ccoordinate meaning
                        if exon.start < variant.pos <= exon.end:  #  same as exon.start+1 <= variant.pos <= exon.end
                            if variant.pos + len(variant.ref) > exon.end:  # variant starts in exon, but ends past exon
                                if i!=0:
                                    raise Exception(
                                    "Variant crossing intron/exon boundary (not edge ones).. Catch this exception and set mutprot=None")
                            trans_pos += (exon.end + 1 - variant.pos)  # position of "end" (most 5', lowest genomic coord.. no need to adjust for variant overlapping end.
                            bef_end_index = exon.length - (variant.pos - exon.start -1)
                            downstream = exonseq[bef_end_index:]  ## if bef_end_index == exon.length ==> ==''
                            aft_start_index = exon.end -variant.pos-len(variant.ref) +1
                            if aft_start_index>0:
                                upstream = exonseq[:aft_start_index]  # if SNP is at end of exon exonseq[:0] == '', correctly
                            else:
                                upstream = Sequence('')


                            if variant.pos+len(variant.ref)-1>exon.end: # This might not even be transcribed
                                sys.stderr.write("CAVA: WARNING: Variants spanning first exon \n")
                                raise Exception("CAVA: WARNING: Variants spanning first (minus strand) exon \n")
                            else:
                                temp_alt = upstream + Sequence(variant.alt).reverseComplement() + downstream
                                temp_ref = upstream + Sequence(variant.ref).reverseComplement() + downstream

                            ret_ref += temp_ref
                            ret += temp_alt
                            if exonseqs is None:
                                ret_exonseqs.append(temp_ref)
                            if trans_pos<coding_start_alt: #Fully "past" variant
                                coding_start_alt += (len(variant.alt) - len(variant.ref))
                            continue
                        elif exon.start < variant.pos + len(variant.ref) <= exon.end:  # variant starts before exon and ends in exon
                            raise Exception(
                                "Variant crossing intron/exon boundary.. Catch this exception and set mutprot=None")
                        elif variant.pos < exon.start and variant.pos + len(
                                variant.ref) > exon.end:  # variant encompasses whole exon
                            raise Exception("Variant changing whole exon .. Catch this exception and set mutprot=None")
                        else:  #
                            pass  # Variant not affecting this exon
                    else: # insertion on minus strand
                        if exon.start+1 <= variant.pos <= exon.end+1:  # same as exon.start+1 <= variant.pos <= exon.end
                            trans_pos += (exon.end + 1 - variant.pos)
                            bef_end_index = exon.length - (variant.pos - exon.start - 1)
                            downstream = exonseq[bef_end_index:]  ## if bef_end_index == exon.length ==> ==''
                            aft_start_index = exon.end - variant.pos + 1
                            if aft_start_index > 0:
                                upstream = exonseq[:aft_start_index]  # if SNP is at end of exon exonseq[:0] == '', correctly
                            else:
                                upstream = Sequence('')
                            temp_alt = upstream + Sequence(variant.alt).reverseComplement() + downstream
                            temp_ref = upstream  + downstream

                            ret_ref += temp_ref
                            ret += temp_alt
                            if exonseqs is None:
                                ret_exonseqs.append(temp_ref)
                            if trans_pos < coding_start_alt:
                                coding_start_alt += len(variant.alt)
                            continue
                        else:  #
                            pass  # Variant not affecting this exon


            trans_pos += abs(exon.end-exon.start)

            if exonseqs is None:
                ret_exonseqs.append(exonseq)

            ret += exonseq
            ret_ref += exonseq

        ret_transcript_alt = ret + ""
        # HS: fixed bug that if variant was before coding start, "ret" would not point to ATG
        ret = ret[coding_start_alt  - 1:]
        if variant is None: # Respect the transcript annotation of the position of the Stop Codon.
            ret = ret[0:self.cds_len]

        self.transcriptseq = ret_ref
        if self.exonseqs is None:
            self.exonseqs = ret_exonseqs
        if exonseqs is None:
            return ret, ret_exonseqs, ret_transcript_alt, ret_ref
        else:
            return ret, exonseqs,     ret_transcript_alt, ret_ref

    def getTranscriptSeq(self,reference):
        if self.transcriptseq is None:
            codingsequencealt, exonseqalt, transcript_alt_seq, transcript_seq = self.getCodingSequence(reference,None,None)
            # self.transcriptseq = transcript_seq  .. this is done inside self.getCodingSequence
        return self.transcriptseq

    # Getting the translated protein sequence of the transcript
    # with the variant included as well as a list of the (reference) exon sequence.
    # The exonseqs should NOT include the alt allele
    # Ideally, call this function first with variant = None to get exonseq .. to use in 2nd call of this function.

    def getProteinSequence(self, reference, variant, exonseqs, codon_usage):
        # Translating coding sequence
        codon_usage = codon_usage
        # returns exonsseq of the alternate allele unless variant is None
        try:
            codingsequencealt, exonseqalt, transcript_alt_seq, transcript_seq = self.getCodingSequence(reference, variant, exonseqs)
        except:
            if variant is None:
                sys.stderr.write("ERROR: exception in getting reference sequence from getCodingSequence")
            return None, None, None # Mutated protein is None when variant crosses intron/exon boundary. 2nd term not used by calling function for variant

        if self.is_selenocysteine: # some plants have selenocysteine genes on their mitochondrion, so selenocysteine genes have to ve first.
            # the 'CDS' annotation of that transcript was trusted.. any TAG remaining will be recoded
            # as Sec/U = Selenocysteine .. and not as 'X'
            # Some 25+ human genes genes recode  the UGA stop codons as selenocysteines,
            # most have name SELENO* so cannot remove the "X" .. must trust the annotation
            #
            # DIO1,DIO2,DIO3,GPX1,GPX2,GPX3,GPX4,GPX6,SELENO[F,H,I,K,M,N,O,P,S,T,U,V,W],MSRB1,SEPHS2,TXNRD1,TXNRD2,TXNRD3
            ret = Sequence(codingsequencealt).translate('7')
            if variant is None and ret[-1] == 'X':
                ret = ret[0:(len(ret)-1)] + 'X'
        elif self.chrom in ['chrM','chrMT' , 'M','MT']: # No selenocysteine genes in mitochondrion genome
            ret = Sequence(codingsequencealt).translate('4')
        else:
            ret = Sequence(codingsequencealt).translate(codon_usage)
        ret = trim_prot_after_stop(ret)  # trim end, keeping from start to and including the first X (if any)

        return ret, exonseqalt, codingsequencealt[0:3*len(ret)]

    # Checking if a given position is outside the region between the start and stop codon
    # don't take account of insertion .. should be dealth with by calling function
    def isPositionOutsideCDS(self, pos):
        if self.strand == 1:
            return (pos < self.codingStartGenomic) or (pos > self.codingEndGenomic)
        else:
            return (pos > self.codingStartGenomic) or (pos < self.codingEndGenomic)

    # Checking if a given position is upstream the region between the start and stop codon
    # don't take account of insertion .. should be dealth with by calling function
    # 9/14/2022 Make sure position is not outside the transcript either.
    def isPositionOutsideCDS_5prime(self, pos):
        return (self.strand == 1 and pos < self.codingStartGenomic and pos>self.transcriptStart) or (
                self.strand == -1 and pos > self.codingStartGenomic and pos<=self.transcriptEnd)

    # Checking if a given position is downstream the region between the start and stop codon
    # don't take account of insertion .. should be dealth with by calling function
    def isPositionOutsideCDS_3prime(self, pos):
        return (self.strand == 1 and pos > self.codingEndGenomic and pos<=self.transcriptEnd) or \
               (self.strand == -1 and pos < self.codingEndGenomic and pos > self.transcriptStart)

    # Checking if the variant is completely (both ends) outside the coding region
    # if FALSE, then part of this inside coding region.
    # if variant is outside transcript, should return true
    def isOutsideTranslatedRegion(self, variant):  # in either UTR's'
        if self.strand == 1:
            if variant.is_insertion: # position points to after variant.
                # reject insertion is it's before start
                if variant.pos <= self.codingStartGenomic: return True
                # if insertion is right next to end of transcript .. it will be cleaved off .. so don't include it.
                if variant.pos - 1 >= self.codingEndGenomic: return True
                return False
            else:
                if variant.pos + len(variant.ref) - 1 < self.codingStartGenomic: return True
                if variant.pos > self.codingEndGenomic : return True
                return False
        else:
            if variant.is_insertion:
                if variant.pos <= self.codingEndGenomic: return True
                if variant.pos - 1 >= self.codingStartGenomic: return True
                return False
            else:
                if variant.pos + len(variant.ref) - 1 < self.codingEndGenomic: return True
                if variant.pos > self.codingStartGenomic: return True
                return False

    # Checking if the given variant is outside of the translated region of the transcript, +/- the first and last 3 bases of the coding sequence
    def isOutsideTranslatedRegionPlus3(self, variant):
        if self.strand == 1:
            if variant.is_insertion:
                if variant.pos <= self.codingStartGenomic:# + 3:
                    return True
                if variant.pos - 1 >= self.codingEndGenomic:# - 3:
                    return True
                return False
            else:
                if variant.pos + len(variant.ref) - 1 < self.codingStartGenomic:# + 3:
                    return True
                if variant.pos > self.codingEndGenomic:# - 3:
                    return True
                return False
        else:
            if variant.is_insertion:
                if variant.pos <= self.codingEndGenomic:# + 3:
                    return True
                if variant.pos - 1 >= self.codingStartGenomic:# - 3:
                    return True
                return False
            else:
                if variant.pos + len(variant.ref) - 1 < self.codingEndGenomic:# + 3:
                    return True
                if variant.pos > self.codingStartGenomic:# - 3:
                    return True
                return False

    # Checking if the given variant overlaps with splicing region, only matters if it affects protein
    def isInSplicingRegion(self, variant, ssrange):
        if not self.isOutsideTranslatedRegion(variant):
            for exon in self.exons:
                if variant.overlap(exon.end + 1, exon.end + ssrange): return True
                if variant.overlap(exon.start - (ssrange - 1), exon.start): return True
            return False
        else:
            return False

    # Checking if the given variant affects an essential splice site, only matters if it affects protein
    def isInEssentialSpliceSite(self, variant):
        if not self.isOutsideTranslatedRegion(variant):
            for exon in self.exons:
                if variant.overlap(exon.end + 1, exon.end + 2): return True
                if variant.overlap(exon.start - 1, exon.start): return True
            return False
        else:
            return False

    # Checking if the given variant affects a +5 essential splice site, only matters if it affects protein
    def isIn_SS5_Site(self, variant):
        if not self.isOutsideTranslatedRegion(variant):
            if self.strand == 1:
                for exon in self.exons:
                    if variant.overlap(exon.end + 1, exon.end + 5):
                        if not variant.is_substitution:
                            if not (variant.pos == exon.end + 3 and len(variant.ref) == 2 and len(
                                    variant.alt) == 2): return True
                        else:
                            if variant.pos == exon.end + 5: return True
                return False
            else:
                for exon in self.exons:
                    if variant.overlap(exon.start - 4, exon.start):
                        if not variant.is_substitution:
                            if not (variant.pos == exon.start - 3 and len(variant.ref) == 2 and len(
                                    variant.alt) == 2): return True
                        else:
                            if variant.pos == exon.start - 4: return True
                return False
        else:
            return False

    # Checking if the given variant affects the first or last 3 bases of an exon
    def isInFirstOrLast3BaseOfExon(self, variant):
        if not self.isOutsideTranslatedRegionPlus3(variant):
            for exon in self.exons:
                if variant.overlap(exon.start + 1, exon.start + 3): return True
                if variant.overlap(exon.end - 2, exon.end): return True
            return False
        else:
            return False

    # Checking where a given genomic position is located in the transcript
    # already adjusted for introns .. because whereIsThisVariant adjusts position
    def whereIsThisPosition(self, pos): # pos is 1-based
        # Iterating through exons and introns and checking if genomic position is located within

        for exon in self.exons:
            if exon.index > 1 and ((self.strand == 1 and prevexonend < pos <= exon.start) or (
                    self.strand == -1 and exon.end < pos <= prevexonend)):
                if self.intronLength(exon.index) > 5 or self.intronLength(exon.index) == 3:
                    return 'In' + str(exon.index - 1) + '/' + str(exon.index)
                else:
                    return 'fsIn' + str(exon.index - 1) + '/' + str(exon.index)
            if exon.start < pos <= exon.end:  # exon.start is 0-based , codingStartGenomic is 1-based
                if (self.strand == 1 and pos < self.codingStartGenomic) or (
                        self.strand == -1 and pos > self.codingStartGenomic):
                    return '5UTR'
                if (self.strand == 1 and pos > self.codingEndGenomic) or (
                        self.strand == -1 and pos < self.codingEndGenomic):
                    return '3UTR'
                return 'Ex' + str(exon.index)
            prevexonend = exon.end if self.strand == 1 else exon.start

        return '.'

    # Checking where a given variant is located in the transcript
    def whereIsThisVariant(self, variant):
        # Getting the locations of both end points of the variant
        if variant.is_insertion:
            first = self.whereIsThisPosition(variant.pos - 1)
            second = self.whereIsThisPosition(variant.pos)
        else:
            first = self.whereIsThisPosition(variant.pos)
            second = self.whereIsThisPosition(variant.pos + len(variant.ref) - 1)
        if first == second:
            return first
        if self.strand == 1:
            if (len(first) == 0 or first == '.') and (len(second) == 0 or second == '.'):
                return '.'
            elif first == '.' or len(first) == 0:
                return second
            elif second == '.' or len(second) == 0:
                return first
            else:
                return first + '-' + second
        else:
            if (len(first) == 0 or first == '.') and (len(second) == 0 or second == '.'):
                return '.'
            elif first == '.' or len(first) ==0:
                return second
            elif second == '.' or len(second)==0:
                return first
            else:
                return second + '-' + first

    # Getting the length of an intron, where idx is the index of the succeeding exon
    def intronLength(self, idx):
        if idx <= 0:
            return 0
        for exon in self.exons:
            if exon.index == idx:
                if self.strand == 1:
                    return exon.start - prev
                else:
                    return prev - exon.end
            if self.strand == 1:
                prev = exon.end
            else:
                prev = exon.start
        return 0


#######################################################################################################################

# Class representing a single exon
class Exon(object):
    # Constructor, Start is 0-based, ENd is 1-based -- like a bed
    def __init__(self, index, start, end):
        self.index = index
        self.start = start
        self.end = end
        self.length = end - start

    def contains(self, pos):
        return self.start + 1 <= pos <= self.end


#######################################################################################################################

# Class representing a DNA sequence, inherits from str class
class Sequence(str):
    # Translating to amino acid sequence
    def translate(self, letter):
        gencode = {}
        if letter == '1':
            gencode = {
                'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
                'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
                'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
                'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
                'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
                'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
                'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
                'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
                'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
                'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
                'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
                'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
                'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
                'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
                'TAC': 'Y', 'TAT': 'Y', 'TAA': 'X', 'TAG': 'X',
                'TGC': 'C', 'TGT': 'C', 'TGA': 'X', 'TGG': 'W'}
        elif letter == '3': # 3-letter code version of 1 (Human)
            gencode = {
                'ATA': 'Ile', 'ATC': 'Ile', 'ATT': 'Ile', 'ATG': 'Met',
                'ACA': 'Thr', 'ACC': 'Thr', 'ACG': 'Thr', 'ACT': 'Thr',
                'AAC': 'Asn', 'AAT': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
                'AGC': 'Ser', 'AGT': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
                'CTA': 'Leu', 'CTC': 'Leu', 'CTG': 'Leu', 'CTT': 'Leu',
                'CCA': 'Pro', 'CCC': 'Pro', 'CCG': 'Pro', 'CCT': 'Pro',
                'CAC': 'His', 'CAT': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
                'CGA': 'Arg', 'CGC': 'Arg', 'CGG': 'Arg', 'CGT': 'Arg',
                'GTA': 'Val', 'GTC': 'Val', 'GTG': 'Val', 'GTT': 'Val',
                'GCA': 'Ala', 'GCC': 'Ala', 'GCG': 'Ala', 'GCT': 'Ala',
                'GAC': 'Asp', 'GAT': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
                'GGA': 'Gly', 'GGC': 'Gly', 'GGG': 'Gly', 'GGT': 'Gly',
                'TCA': 'Ser', 'TCC': 'Ser', 'TCG': 'Ser', 'TCT': 'Ser',
                'TTC': 'Phe', 'TTT': 'Phe', 'TTA': 'Leu', 'TTG': 'Leu',
                'TAC': 'Tyr', 'TAT': 'Tyr', 'TAA': 'X', 'TAG': 'X',
                'TGC': 'Cys', 'TGT': 'Cys', 'TGA': 'X', 'TGG': 'Trp'}
        elif letter == '4': # Mitochondrgenetic code for vertebrates
            gencode = {
                'ATA': 'M', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
                'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
                'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
                'AGC': 'S', 'AGT': 'S', 'AGA': 'X', 'AGG': 'X',
                'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
                'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
                'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
                'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
                'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
                'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
                'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
                'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
                'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
                'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
                'TAC': 'Y', 'TAT': 'Y', 'TAA': 'X', 'TAG': 'X',
                'TGC': 'C', 'TGT': 'C', 'TGA': 'W', 'TGG': 'W'}
        elif letter == '5': # 3-letter version of Mitochondria
            gencode = {
                'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAT': 'N',
                'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
                'AGA': 'R', 'AGC': 'S', 'AGG': 'R', 'AGT': 'S',
                'ATA': 'I', 'ATC': 'I', 'ATG': 'M', 'ATT': 'I',
                'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAT': 'H',
                'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
                'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
                'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
                'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAT': 'D',
                'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
                'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
                'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
                'TAA': 'X', 'TAC': 'Y', 'TAG': 'X', 'TAT': 'Y',
                'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
                'TGA': 'X', 'TGC': 'C', 'TGG': 'W', 'TGT': 'C',
                'TTA': 'L', 'TTC': 'F', 'TTG': 'L', 'TTT': 'F'}
        elif letter == '7': # Human genome for one of the 25 selenoproteins (TAG is Se-Cys (U/Sec)
            gencode = {
                'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
                'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
                'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
                'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
                'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
                'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
                'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
                'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
                'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
                'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
                'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
                'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
                'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
                'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
                'TAC': 'Y', 'TAT': 'Y', 'TAA': 'X', 'TAG': 'X',
                'TGC': 'C', 'TGT': 'C', 'TGA': 'U', 'TGG': 'W'}
        elif letter == '8': # 3-letter code version of 1 (Human) with SelenoCysteinr for ATG
            gencode = {
                'ATA': 'Ile', 'ATC': 'Ile', 'ATT': 'Ile', 'ATG': 'Met',
                'ACA': 'Thr', 'ACC': 'Thr', 'ACG': 'Thr', 'ACT': 'Thr',
                'AAC': 'Asn', 'AAT': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
                'AGC': 'Ser', 'AGT': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
                'CTA': 'Leu', 'CTC': 'Leu', 'CTG': 'Leu', 'CTT': 'Leu',
                'CCA': 'Pro', 'CCC': 'Pro', 'CCG': 'Pro', 'CCT': 'Pro',
                'CAC': 'His', 'CAT': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
                'CGA': 'Arg', 'CGC': 'Arg', 'CGG': 'Arg', 'CGT': 'Arg',
                'GTA': 'Val', 'GTC': 'Val', 'GTG': 'Val', 'GTT': 'Val',
                'GCA': 'Ala', 'GCC': 'Ala', 'GCG': 'Ala', 'GCT': 'Ala',
                'GAC': 'Asp', 'GAT': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
                'GGA': 'Gly', 'GGC': 'Gly', 'GGG': 'Gly', 'GGT': 'Gly',
                'TCA': 'Ser', 'TCC': 'Ser', 'TCG': 'Ser', 'TCT': 'Ser',
                'TTC': 'Phe', 'TTT': 'Phe', 'TTA': 'Leu', 'TTG': 'Leu',
                'TAC': 'Tyr', 'TAT': 'Tyr', 'TAA': 'X', 'TAG': 'X',
                'TGC': 'Cys', 'TGT': 'Cys', 'TGA': 'Sec', 'TGG': 'Trp'}
        ret_list = []
        index = 0
        while index + 3 <= len(self):
            #codon = self[index:index + 3].upper()   # Upper is not necessary .. it is done in Referenge.getReference
            codon = self[index:index + 3]
# SLOWEST line of code.           if codon.replace('A', '').replace('C', '').replace('G', '').replace('T', ''):
            if codon not in gencode:
                # if 'N' in codon:
                ret_list.append('?')
                index += 3
                continue
            ret_list.append(gencode[codon])
            index += 3
        return ''.join(ret_list)

    # Getting reverse complement sequence
    def reverseComplement(self):
        # 5/13/21 Added support for reverse complement of "."
        complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N", "a": "t", "t": "a", "c": "g", "g": "c",
                      "n": "n", "*": "*", ".": "."}
        ret = ''
        for base in self[::-1]:
            ret += complement[base]
        return ret


#######################################################################################################################

# Class representing the collection of options specified in the configuration file
class Options(object):
    # Constructor
    def __init__(self, configfn):

        self.args = dict()
        self.defs = dict()
        self.configfn = configfn
        self.transcript2protein = dict()

        # Defining option flags
        self.defs['reference'] = ('string', '.')
        self.defs['ensembl'] = ('string', '.')
        self.defs['dbsnp'] = ('string', '.')
        self.defs['filter'] = ('boolean', False)
        self.defs['inputformat'] = ('string', 'VCF')
        self.defs['outputformat'] = ('string', 'VCF')
        self.defs['type'] = ('string', 'ALL')
        self.defs['logfile'] = ('boolean', False)
        self.defs['target'] = ('string', '.')
        self.defs['genelist'] = ('string', '.')
        self.defs['transcriptlist'] = ('string', '.')
        self.defs['snplist'] = ('string', '.')
        self.defs['nonannot'] = ('boolean', True)
        self.defs['givealt'] = ('boolean', True)
        self.defs['givealtflag'] = ('boolean', True)
        self.defs['ssrange'] = ('string', '8')
        self.defs['ontology'] = ('string', 'both')
        self.defs['impactdef'] = ('string', 'SG,ESS,FS|SS5,IM,SL,EE,IF,NSY|SY,SS,INT,5PU,3PU')
        self.defs['prefix'] = ('boolean', False)
        self.defs['codon_usage'] = ('string', '1')
        self.defs['normalized_mitochondrial_chrom'] = ('string', 'not_normalized')

        self.defs['transcript2protein'] = ('string', '.')
        self.defs['loadalltranscripts'] = ('boolean', True)

        # Reading options from file
        self.read()

        if self.args['ensembl'] == '.' or self.args['ensembl'] == '':
            d = os.path.dirname(os.path.realpath(__file__))
            dirn = d[:d.rfind('env/lib')]
            self.args['ensembl'] = dirn + '/defaultdb/ensembl75s.gz'

    # Reading options from configuration file
    def read(self):
        for line in open(self.configfn):
            line = line.strip()
            if line.startswith('@'):
                key = line[1:line.index('=')].strip()
                if key in list(self.defs.keys()):
                    (typeofvar, default) = self.defs[key]
                    if typeofvar == 'string': self.args[key] = line[line.find('=') + 1:].strip()
                    if typeofvar == 'list': self.args[key] = line[line.find('=') + 1:].strip().split(',')
                    if typeofvar == 'boolean': self.args[key] = (line[line.find('=') + 1:].strip().upper() == 'TRUE')
        for key, (typeofvar, default) in self.defs.items():
            if not key in list(self.args.keys()): self.args[key] = default


#######################################################################################################################            

# Other basic utility functions 

#
# Convert various chromosome nomenclatures to that found in file.
#
def convert_chrom(chrom, contigs):
        # Checking if chromosome name exists
        if chrom in contigs:
            return chrom
        if len(chrom)>3 and chrom.startswith("chr"):
            goodchrom = chrom[3:]
            if goodchrom in contigs:
                return goodchrom

            if goodchrom == "MT":
                goodchrom = 'chrM'
                if goodchrom in contigs:
                    return goodchrom
                goodchrom = "M"
                if goodchrom in contigs:
                    return goodchrom
                else:
                    return None
        elif chrom == "MT": # Cannot be chrMT .. because those starting with chr dealth above
            goodchrom = 'chrMT'
            if goodchrom in contigs:
                return goodchrom
            goodchrom = 'chrM'
            if goodchrom in contigs:
                return goodchrom
            goodchrom = "M"
            if goodchrom in contigs:
                return goodchrom
        elif chrom == "M":  # Cannot be chrM .. because those starting with chr dealth above
                goodchrom = "MT"
                if goodchrom in contigs:
                    return goodchrom
                goodchrom = "chrMT"
                if goodchrom in contigs:
                    return goodchrom
                goodchrom = "chrM"
                if goodchrom in contigs:
                    return goodchrom
        else:
            goodchrom = "chr"+chrom
            if goodchrom in contigs:
                return goodchrom
        return None


# Reading gene, transcript or snp list from file
def readSet(options, tag):
    ret = set()
    if tag in list(options.args.keys()) and not (options.args[tag] == '' or options.args[tag] == '.'):
        for line in open(options.args[tag]):
            line = line.strip()
            if line == '' or line == '.':
                continue
            ret.add(line)
        if options.args['logfile']:
            txt = ''
            if tag == 'genelist':
                txt = 'Gene list'
            if tag == 'transcriptlist':
                txt = 'Transcript list'
            if tag == 'snplist':
                txt = 'SNP list'
            logging.info(txt + ' loaded.')
    return ret


# Reading dictionary (from transcript -> Protein)
def read_dict(options, tag):
    ret = dict()
    tx_to_prot_source_file = options.args['ensembl'].replace('db.gz', 'txt').replace('.gz', '.txt')
    print("Reading option file dictionary : " + tx_to_prot_source_file + "\n")
    with open(tx_to_prot_source_file, 'r') as f:
        for line in f:
            if line == '' or line == '.' or line.startswith("#"): continue
            linedat = line.strip().split("\t")
            if len(linedat) == 3:
                ret[linedat[2]] = '.'
            else:
                ret[linedat[2]] = linedat[3]

    if options.args['logfile']:
        txt = options.args[tag]
        logging.info(txt + ' dictionary loaded from file ' + tx_to_prot_source_file)
    return ret


# Writing header information to output file
def writeHeader(options, header, outfile, stdout, version):
    if options.args['prefix']:
        prefix = 'CAVA_'
    else:
        prefix = ''
    headerinfo = '##INFO=<ID=' + prefix + 'TYPE,Number=.,Type=String,Description=\"Variant type: Substitution, Insertion, Deletion or Complex\",Source=\"CAVA\",Version=\"' + version + '\">\n'
    headerinfo += '##INFO=<ID=' + prefix + 'GENE,Number=.,Type=String,Description=\"HGNC gene symbol\",Source=\"CAVA\",Version=\"' + version + '\">\n'
    headerinfo += '##INFO=<ID=' + prefix + 'TRANSCRIPT,Number=.,Type=String,Description=\"Transcript identifier\",Source=\"CAVA\",Version=\"' + version + '\">\n'
    headerinfo += '##INFO=<ID=' + prefix + 'GENEID,Number=.,Type=String,Description=\"Gene identifier\",Source=\"CAVA\",Version=\"1.2.4\">\n'
    headerinfo += '##INFO=<ID=' + prefix + 'TRINFO,Number=.,Type=String,Description=\"Transcript information: Strand/Length of transcript/Number of exons/Length of coding DNA + UTR/Protein length\",Source=\"CAVA\",Version=\"' + version + '\">\n'
    headerinfo += '##INFO=<ID=' + prefix + 'LOC,Number=.,Type=String,Description=\"Location of variant in transcript\",Source=\"CAVA\",Version=\"' + version + '\">\n'
    headerinfo += '##INFO=<ID=' + prefix + 'CSN,Number=.,Type=String,Description=\"CSN annotation\",Source=\"CAVA\",Version=\"' + version + '\">\n'
    headerinfo += '##INFO=<ID=' + prefix + 'PROTPOS,Number=.,Type=String,Description=\"Protein position\",Source=\"CAVA\",Version=\"' + version + '\">\n'
    headerinfo += '##INFO=<ID=' + prefix + 'PROTREF,Number=.,Type=String,Description=\"Reference amino acids\",Source=\"CAVA\",Version=\"' + version + '\">\n'
    headerinfo += '##INFO=<ID=' + prefix + 'PROTALT,Number=.,Type=String,Description=\"Alternate amino acids\",Source=\"CAVA\",Version=\"' + version + '\">\n'
    headerinfo += '##INFO=<ID=' + prefix + 'CLASS,Number=.,Type=String,Description=\"5PU: Variant in 5 prime untranslated region, 3PU: Variant in 3 prime untranslated region, INT: Intronic variant that does not alter splice site bases, SS: Intronic variant that alters a splice site base but not an ESS or SS5 base, ESS: Variant that alters essential splice site base (+1,+2,-1,-2), SS5: Variant that alters the +5 splice site base, but not an ESS base, SY: Synonymous change caused by a base substitution (i.e. does not alter amino acid), NSY: Nonsynonymous change (missense) caused by a base substitution (i.e. alters amino acid), IF: Inframe insertion and/or deletion (variant alters the length of coding sequence but not the frame), IM: Variant that alters the start codon, SG: Variant resulting in stop-gain (nonsense) mutation, SL: Variant resulting in stop-loss mutation, FS: Frameshifting insertion and/or deletion (variant alters the length and frame of coding sequence), EE: Inframe deletion, insertion or base substitution which affects the first or last three bases of the exon\",Source=\"CAVA\",Version=\"' + version + '\">\n'
    headerinfo += '##INFO=<ID=' + prefix + 'SO,Number=.,Type=String,Description=\"Sequence Ontology term\",Source=\"CAVA\",Version=\"' + version + '\">\n'
    headerinfo += '##INFO=<ID=' + prefix + 'ALTFLAG,Number=.,Type=String,Description=\"None: variant has the same CSN annotation regardless of its left or right-alignment, AnnNotClass/AnnNotSO/AnnNotClassNotSO: indel has an alternative CSN but the same CLASS and/or SO, AnnAndClass/AnnAndSO/AnnAndClassNotSO/AnnAndSONotClass/AnnAndClassAndSO: Multiple CSN with different CLASS and/or SO annotations\",Source=\"CAVA\",Version=\"' + version + '\">\n'
    headerinfo += '##INFO=<ID=' + prefix + 'ALTANN,Number=.,Type=String,Description=\"Alternate CSN annotation\",Source=\"CAVA\",Version=\"' + version + '\">\n'
    headerinfo += '##INFO=<ID=' + prefix + 'ALTCLASS,Number=.,Type=String,Description=\"Alternate CLASS annotation\",Source=\"CAVA\",Version=\"' + version + '\">\n'
    headerinfo += '##INFO=<ID=' + prefix + 'ALTSO,Number=.,Type=String,Description=\"Alternate SO annotation\",Source=\"CAVA\",Version=\"' + version + '\">\n'
    headerinfo += '##INFO=<ID=' + prefix + 'IMPACT,Number=.,Type=String,Description=\"Impact group the variant is stratified into\",Source=\"CAVA\",Version=\"' + version + '\">\n'
    headerinfo += '##INFO=<ID=' + prefix + 'DBSNP,Number=.,Type=String,Description=\"rsID from dbSNP\",Source=\"CAVA\",Version=\"' + version + '\">\n'
    headerinfo += '##INFO=<ID=' + prefix + 'HGVSc,Number=.,Type=String,Description=\"HGVS Nomenclature for cDNA changes\",Source=\"CAVA\",Version=\"' + version + '\">\n'
    headerinfo += '##INFO=<ID=' + prefix + 'HGVSp,Number=.,Type=String,Description=\"HGVS Nomenclature for protein changes\",Source=\"CAVA\",Version=\"' + version + '\">\n'
    headerinfo += '##INFO=<ID=' + prefix + 'HGVSg,Number=.,Type=String,Description=\"HGVS Nomenclature for genomic changes, right-shifted\",Source=\"CAVA\",Version=\"' + version + '\">\n'

    dateline = '##fileDate=' + time.strftime("%Y-%m-%d")

    if options.args['outputformat'] == 'VCF':
        if header == '':
            if stdout:
                print(
                    '##fileformat=VCFv4.1\n' + dateline + '\n' + headerinfo + '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO')
            else:
                outfile.write(
                    '##fileformat=VCFv4.1\n' + dateline + '\n' + headerinfo + '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        else:
            headerv = header.split('\n')
            headerm = []
            startline = '##fileformat=VCFv4.1'
            for x in headerv:
                if x.startswith('##fileDate') or x.startswith('##filedate'):
                    continue
                if x.startswith('##fileformat'):
                    startline = x
                    continue
                headerm.append(x)
            if stdout:
                print(startline + '\n' + dateline + '\n' + headerinfo + '\n'.join(headerm) + '\n')
            else:
                outfile.write(startline + '\n' + dateline + '\n' + headerinfo + '\n'.join(headerm) + '\n')

    if options.args['outputformat'] == 'TSV':
        hstr = 'ID\tCHROM\tPOS\tREF\tALT\tQUAL\tFILTER\tTYPE'
        if (not options.args['ensembl'] == '.') and (not options.args['ensembl'] == ''):
            if options.args['ontology'].upper() == 'CLASS':
                hstr += '\tTRANSCRIPT\tGENE\tGENEID\tTRINFO\tLOC\tCSN\tPROTPOS\tPROTREF\tPROTALT\tCLASS'
            if options.args['ontology'].upper() == 'SO':
                hstr += '\tTRANSCRIPT\tGENE\tGENEID\tTRINFO\tLOC\tCSN\tPROTPOS\tPROTREF\tPROTALT\tSO'
            if options.args['ontology'].upper() == 'BOTH':
                hstr += '\tTRANSCRIPT\tGENE\tGENEID\tTRINFO\tLOC\tCSN\tPROTPOS\tPROTREF\tPROTALT\tCLASS\tSO'

            if not (options.args['impactdef'] == '.' or options.args['impactdef'] == ''):
                hstr += '\tIMPACT'

            if options.args['givealt']:
                if options.args['ontology'].upper() == 'CLASS': hstr += '\tALTANN\tALTCLASS'
                if options.args['ontology'].upper() == 'SO': hstr += '\tALTANN\tALTSO'
                if options.args['ontology'].upper() == 'BOTH': hstr += '\tALTANN\tALTCLASS\tALTSO'

            if (not options.args['givealt']) or options.args['givealtflag']:
                hstr += '\tALTFLAG'

        if (not options.args['dbsnp'] == '.') and (not options.args['dbsnp'] == ''):
            hstr += '\tDBSNP'

        hstr += '\tHGVSG\tHGVSC\tHGVSP'

        if stdout:
            print(hstr)
        else:
            outfile.write(hstr + '\n')


# Counting number of records in a file
def countRecords(filename):
    ret = 0
    if filename.endswith('.gz') or filename.endswith('.bgz'):
        inputf = gzip.open(filename, 'rt',encoding='utf-8')
    else:
        inputf = open(filename,encoding="utf-8")
    for line in inputf:
        line = line.strip()
        if not (line.startswith("#") or line == ''): ret += 1
    return ret


# Checking if options are correct
def checkOptions(options):
    # Checking if @inputformat was given correct value
    optstr = options.args['inputformat'].upper()
    if not (optstr == 'VCF' or optstr == 'TXT'):
        print('ERROR: incorrect value of the tag @inputformat.')
        print('(Allowed values: \'VCF\' or \'TXT\')')
        print('\nNo output file written. CAVA quit.')
        print("--------------------------------------------------------------------\n")
        if options.args['logfile']:
            logging.error('Incorrect value of the tag @inputformat.')
            logging.info('No output file written. CAVA quit.')
        quit()

    # Checking if @outputformat was given correct value
    optstr = options.args['outputformat'].upper()
    if not (optstr == 'VCF' or optstr == 'TSV'):
        print('ERROR: incorrect value of the tag @outputformat.')
        print('(Allowed values: \'VCF\' or \'TSV\')')
        print('\nNo output file written. CAVA quit.')
        print("--------------------------------------------------------------------\n")
        if options.args['logfile']:
            logging.error('Incorrect value of the tag @outputformat.')
            logging.info('No output file written. CAVA quit.')
        quit()

    # Checking if @type was given correct value
    optstr = options.args['type'].upper()
    if not (
            optstr == 'ALL' or optstr == 'SUBSTITUTION' or optstr == 'INDEL' or optstr == 'INSERTION' or optstr == 'DELETION' or optstr == 'COMPLEX'):
        print('ERROR: incorrect value of the tag @type.')
        print('(Allowed values: \'all\', \'substitution\', \'indel\', \'insertion\', \'deletion\' or \'complex\')')
        print('\nNo output file written. CAVA quit.')
        print("--------------------------------------------------------------------\n")
        if options.args['logfile']:
            logging.error('Incorrect value of the tag @type.')
            logging.info('No output file written. CAVA quit.')
        quit()

    # Checking if @ssrange was given correct value
    ssrange = int(options.args['ssrange'])
    if not ssrange >= 6:
        print('ERROR: incorrect value of the tag @ssrange.')
        print('(Minimum value allowed is 6.)')
        print('\nNo output file written. CAVA quit.')
        print("--------------------------------------------------------------------\n")
        if options.args['logfile']:
            logging.error('Incorrect value of the tag @ssrange.')
            logging.info('No output file written. CAVA quit.')
        quit()

    # Checking if @ontology was given correct value
    optstr = options.args['ontology'].upper()
    if not (optstr == 'CLASS' or optstr == 'SO' or optstr == 'BOTH'):
        print('ERROR: incorrect value of the tag @ontology.')
        print('(Allowed values: \'CLASS\' or \'SO\' or \'both\')')
        print('\nNo output file written. CAVA quit.')
        print("--------------------------------------------------------------------\n")
        if options.args['logfile']:
            logging.error('Incorrect value of the tag @ontology.')
            logging.info('No output file written. CAVA quit.')
        quit()

    # Checking if @reference file exists
    if not os.path.isfile(options.args['reference']):
        print('ERROR: the file given as @reference does not exist.')
        print('\nNo output file written. CAVA quit.')
        print("--------------------------------------------------------------------\n")
        if options.args['logfile']:
            logging.error('The file given as @reference does not exist.')
            logging.info('No output file written. CAVA quit.')
        quit()

    # Checking if @reference index file exists
    if not os.path.isfile(options.args['reference'] + '.fai'):
        print('ERROR: the .fa.fai index file for @reference is not found.')
        print('\nNo output file written. CAVA quit.')
        print("--------------------------------------------------------------------\n")
        if options.args['logfile']:
            logging.error('The .fa.fai index file for @reference is not found.')
            logging.info('No output file written. CAVA quit.')
        quit()

    # Checking if @ensembl file exists
    if not (options.args['ensembl'] == '.' or options.args['ensembl'] == '') and not os.path.isfile(
            options.args['ensembl']):
        print('ERROR: the file given as @ensembl does not exist.')
        print('\nNo output file written. CAVA quit.')
        print("--------------------------------------------------------------------\n")
        if options.args['logfile']:
            logging.error('The file given as @ensembl does not exist.')
            logging.info('No output file written. CAVA quit.')
        quit()

    # Checking if @ensembl index file exists
    if not (options.args['ensembl'] == '.' or options.args['ensembl'] == '') and not os.path.isfile(
            options.args['ensembl'] + '.tbi'):
        print('ERROR: the .gz.tbi index file for @ensembl is not found.')
        print('\nNo output file written. CAVA quit.')
        print("--------------------------------------------------------------------\n")
        if options.args['logfile']:
            logging.error('The .gz.tbi index file for @ensembl is not found.')
            logging.info('No output file written. CAVA quit.')
        quit()

    # Checking if @dbsnp file exists
    if not (options.args['dbsnp'] == '.' or options.args['dbsnp'] == '') and not os.path.isfile(options.args['dbsnp']):
        print('ERROR: the file given as @dbsnp does not exist.')
        print('\nNo output file written. CAVA quit.')
        print("--------------------------------------------------------------------\n")
        if options.args['logfile']:
            logging.error('The file given as @dbsnp does not exist.')
            logging.info('No output file written. CAVA quit.')
        quit()

    # Checking if @dbsnp index file exists
    if not (options.args['dbsnp'] == '.' or options.args['dbsnp'] == '') and not os.path.isfile(
            options.args['dbsnp'] + '.tbi'):
        print('ERROR: the .gz.tbi index file for @dbsnp is not found.')
        print('\nNo output file written. CAVA quit.')
        print("--------------------------------------------------------------------\n")
        if options.args['logfile']:
            logging.error('The .gz.tbi index file for @dbsnp is not found.')
            logging.info('No output file written. CAVA quit.')
        quit()

    # Checking if @target file exists
    if not (options.args['target'] == '.' or options.args['target'] == '') and not os.path.isfile(
            options.args['target']):
        print('ERROR: the file given as @target does not exist.')
        print('\nNo output file written. CAVA quit.')
        print("--------------------------------------------------------------------\n")
        if options.args['logfile']:
            logging.error('The file given as @target does not exist.')
            logging.info('No output file written. CAVA quit.')
        quit()

    # Checking if @target index file exists
    if not (options.args['target'] == '.' or options.args['target'] == '') and not os.path.isfile(
            options.args['target'] + '.tbi'):
        print('ERROR: the .bed.tbi index file for @target is not found.')
        print('\nNo output file written. CAVA quit.')
        print("--------------------------------------------------------------------\n")
        if options.args['logfile']:
            logging.error('The .bed.tbi index file for @target is not found.')
            logging.info('No output file written. CAVA quit.')
        quit()

    # Checking if @genelist file exists
    if not (options.args['genelist'] == '.' or options.args['genelist'] == '') and not os.path.isfile(
            options.args['genelist']):
        print('ERROR: the file given as @genelist does not exist.')
        print('\nNo output file written. CAVA quit.')
        print("--------------------------------------------------------------------\n")
        if options.args['logfile']:
            logging.error('The file given as @genelist does not exist.')
            logging.info('No output file written. CAVA quit.')
        quit()

    # Checking if @transcriptlist file exists
    if not (options.args['transcriptlist'] == '.' or options.args['transcriptlist'] == '') and not os.path.isfile(
            options.args['transcriptlist']):
        print('ERROR: the file given as @transcriptlist does not exist.')
        print('\nNo output file written. CAVA quit.')
        print("--------------------------------------------------------------------\n")
        if options.args['logfile']:
            logging.error('The file given as @transcriptlist does not exist.')
            logging.info('No output file written. CAVA quit.')
        quit()

    # Checking if @snplist file exists
    if not (options.args['snplist'] == '.' or options.args['snplist'] == '') and not os.path.isfile(
            options.args['snplist']):
        print('ERROR: the file given as @snplist does not exist.')
        print('\nNo output file written. CAVA quit.')
        print("--------------------------------------------------------------------\n")
        if options.args['logfile']:
            logging.error('The file given as @snplist does not exist.')
            logging.info('No output file written. CAVA quit.')
        quit()

    # Checking if @transcript2protein file exists
    #  
    if not (('transcripts2protein' not in options.args) or options.args['transcript2protein'] == '.'
            or options.args['transcript2protein'] == '') and not os.path.isfile(options.args['transcript2protein']):
        print('ERROR: the file given as @transcript2protein does not exist.')
        print('\nNo output file written. CAVA quit.')
        print("--------------------------------------------------------------------\n")
        quit()

    # Checking if @normalized_mitochondrial_chrom was given correct value
    optstr = options.args['normalized_mitochondrial_chrom']
    if optstr != 'not_normalized' and not (optstr == 'M' or optstr == 'MT'):
        print('ERROR: incorrect value of the tag @normalized_mitochondrial_chrom.')
        print('(Allowed values: \'M\' or \'MT\')')
        print('\nNo output file written. CAVA quit.')
        print("--------------------------------------------------------------------\n")
        if options.args['logfile']:
            logging.error('Incorrect value of the tag @normalized_mitochondrial_chrom.')
            logging.info('No output file written. CAVA quit.')
        quit()

# trim protein, up to including the stop codon .. to fix annotation errors.

def trim_prot_after_stop(seq):
    if 'X' in seq:
        seq = seq[0:(seq.index('X') + 1)]
    elif 'x' in seq:
        seq = seq[0:(seq.index('x') + 1)]
    return seq

######################################################################################################################
