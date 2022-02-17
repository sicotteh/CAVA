import pysam
import sys


#
# Convert various chromosome nomenclatures to that found in file.
#
def convert_chrom(chrom, contigs):
        # Checking if chromosome name exists
        if chrom in contigs:
            return chrom
        if chrom.startswith("chr") and len(chrom)>3:
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



# Class representing the reference genome dataset
class Reference(object):
    # Constructor
    def __init__(self, fasta_reference):
        # Openning faidx reference fasta representing the reference genome
        #
        # Switch from hasattr to try except to be compatible with python2 version of library
        try:
            self.fastafile = pysam.FastaFile(fasta_reference)
        except:  # old API, version prior to 0.8.1
            try:
                self.fastafile = pysam.Fastafile(fasta_reference)
            except:
                raise Exception("pysam library API changes, cannot find Fastafile of FastaFile. Use older version\n")
        # Need pysam >= 0.7.7 for lengths and references to work.
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

            if chrom in ["chrM","chrMT","M","MT","chrM"]:
                self.reflens["chrM"] = lengths[i]
                self.reflens["chrMT"] = lengths[i]
                self.reflens["M"] = lengths[i]
                self.reflens["MT"] = lengths[i]
                self.reflens["chrm"] = lengths[i]
                self.reflens["chrmt"] = lengths[i]
                self.reflens["m"] = lengths[i]
                self.reflens["mt"] = lengths[i]

        self.cache = ""
# caching will only wins if get more than 4-5 cache hits/5MB. Will not win for indel-only normalization of exome without off-target elements.
#  in all other cases, it will win.
        self.CACHESIZE = 5000000
        self.PAD = 5000  # Default Left Padding for normalization of large indels(5000bp) .. so first cache hit will have padding.
        self.chrom = ""
        self.start0 = 0
        self.endpos = 0

    # Private method .. do not call because assume that chrom is in sequence index
    #  also assumes end is less than chrom length.

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
        if self.CACHESIZE == 0:
            return  self.fastafile.fetch(chrom, start-1, endpos)
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
            self.cache = self.fastafile.fetch(chrom, self.start0, self.endpos)

        return self.cache[(start - 1 - self.start0):(endpos  - self.startpos0)]

    # Retrieving the sequence of a genomic region
    def getReference(self, chrom, start, end):
        # XXX-HS To make faster, could retrieve large blocks of sequence
        # For normalization.
        # There are about 322,000 indels per 3XE9 /person, so 0.3/1000bp
        # so for a cache to be useful, it needs to include variants, so 30KB min.
        #  for exomes. Most genes are < 1MB .. and 2000 bp cds.. so would only include 1-2 indels/sample
        # Assuming the tabix index is in RAM, reading from disk is seek time (5-10ms) + read time (1MB/ms)
        #    so reading 10-20MB takes the same amount of time as a single disk seek.
        # For exome, there is one indel/gene and each gene is 1MB (average), so a 5MB cache insures 5X cache reuse(hits)
        # if dealing with SNPs or whole genome data, other elements, a smaller cache will reap benefits too .. up to about 5MB.

        # Checking if chromosome name exists
        goodchrom = convert_chrom(chrom,self.fastafile.references)
        if goodchrom is None:
            return ''
            # Fetching data from reference genome
        if end < start:
            return ''
        if start < 1:
            start = 1

        try:
            last = self.fastafile.get_reference_length(goodchrom)
        except:
            try:
                last = self.fastafile.getReferenceLength(goodchrom)  # if pysam.__version__ in ['0.7.7', '0.7.8', '0.8.0']:
            except:
                sys.stderr.write("API for pysam changed, may have issues at ends of genome")
                last = end

        if end > last:
            end = last

 #       seq = self.fastafile.fetch(goodchrom, start - 1, end)

        seq = self.__getseq_from_cache_or_file(goodchrom, start, end)
        return seq.upper()


# Class representing a single variant call
class Variant(object):
    # Constructor
    def __init__(self, chrom, pos, vcf_ref, vcf_alt):
        vcf_ref = vcf_ref.upper()
        vcf_alt = vcf_alt.upper()
        self.id = chrom + ":" + str(pos)+":" +vcf_ref + "/"+vcf_alt
        self.chrom = chrom
        self.pos = pos
        self.vcf_padded_base = ''
        shift_pos, a, b = self.trimCommonStart(vcf_ref, vcf_alt)
        self.pos = self.pos + shift_pos
        shiftpos2, x, y = self.trimCommonEnd(a, b)
        if (len(x) == 0 or len(y)==0)  and shift_pos>0:
            self.vcf_padded_base = vcf_ref[shift_pos-1]
        self.ref = x
        self.alt = y

    def get_id(self):
        return self.id
    def get_vcf_pos_ref_alt(self):
        if self.is_insertion():
            return [self.pos-1,self.vcf_padded_base+self.ref,self.vcf_padded_base+self.alt]
        elif self.is_deletion():
            return [self.pos-1,self.vcf_padded_base+self.ref,self.vcf_padded_base+self.alt]
        else:
            return [self.pos, self.ref, self.alt]  # for SNPs or fixed MNPs.



    # Checking if the variant is an insertion
    def is_insertion(self):
        return len(self.ref) == 0 and len(self.alt) > 0

    # Checking if the variant is a deletion
    def is_deletion(self):
        return len(self.ref) > 0 and len(self.alt) == 0

    # Checking if the variant overlaps with a genomic region
    def overlap(self, start, end):
        if self.is_insertion():
            return (self.pos - 1 >= start) and (self.pos <= end)
        int1_start = self.pos
        int1_end = self.pos + len(self.ref) - 1
        int2_start = start
        int2_end = end
        if int1_start == int2_start:
            return True
        if int1_start > int2_start:
            return int2_end >= int1_start
        else:
            return int1_end >= int2_start

    def remove_same_bases_ref_alt(self):  # Only apply for SNV or MNP with padded bases to remove.
        if len(self.ref) != len(self.alt) or len(self.ref)==0:
            return [self.pos,self.ref,self.alt]
        if self.ref == self.alt:
                return [self.pos, self.ref[0], self.alt[0]]
        left, seq1, seq2 = self.trimCommonStart(self.ref, self.alt)
        right, newref, newalt = self.trimCommonEnd(seq1, seq2)
        return [self.pos+left,newref,newalt]


    # Aligning variant on the plus strand
    def alignOnPlusStrand(self, reference):
        if self.chrom not in reference.reflens or not (self.is_insertion() or self.is_deletion()):
            [self.pos,self.ref,self.alt] = self.remove_same_bases_ref_alt()    # Trime redundant bases from SNP/MNP
            return self
        reflen = reference.reflens[self.chrom]
        maxreplen = abs(len(self.ref) - len(self.alt))
        PADLEN = max(100, 1 + 5 * maxreplen)  # the bigger the repeat, the less likely there are multiple exact copies, so 5 is OKfSequence
#        seq1_0 = reference.getReference(self.chrom, self.pos, self.pos + len(self.ref) - 1 + PADLEN)
#        seq2_0 = self.alt + reference.getReference(self.chrom, self.pos + len(self.ref),
#                                                   self.pos + len(self.ref) + PADLEN - 1)
        maxpos = self.pos + len(self.ref) - 1 + PADLEN
        if maxpos > reflen:
            maxpos = reflen
        EXTRA_PREBASE = 1
        if self.pos==1: # Cannot have extra VCF padding base at first base
            EXTRA_PREBASE = 0
        seq12_0 = reference.getReference(self.chrom, self.pos-EXTRA_PREBASE, maxpos)
        # seq12_0 now includes an extra padded base
        seq1_0 = seq12_0[EXTRA_PREBASE:]
        seq2_0 = self.alt + seq12_0[len(self.ref)+EXTRA_PREBASE:]
        left, seq1, seq2 = self.rightAlign(seq1_0, seq2_0)
        while left >= PADLEN - (maxreplen - 1):  # Shifted almost all the way to the end, make sure there is not more to shift
            PADLEN += PADLEN
            maxpos = self.pos + len(self.ref) - 1 + PADLEN
            if maxpos > reflen:
                maxpos = reflen
#            seq1_0 = reference.getReference(self.chrom, self.pos, self.pos + len(self.ref) - 1 + PADLEN)
#            seq2_0 = self.alt + reference.getReference(self.chrom, self.pos + len(self.ref),
#                                                       self.pos + len(self.ref) + PADLEN - 1)
            seq12_0 = reference.getReference(self.chrom, self.pos-EXTRA_PREBASE, maxpos)
            seq1_0 = seq12_0[EXTRA_PREBASE:]
            seq2_0 = self.alt + seq12_0[len(self.ref)+EXTRA_PREBASE:]
            left, seq1, seq2 = self.rightAlign(seq1_0, seq2_0)
        if len(seq1) == 0 or len(seq2) == 0: # e.g. if insertion/deletion .. should always be the case
            if(left >= 1 and EXTRA_PREBASE == 1):
                #base = reference.getReference(self.chrom, self.pos + left, self.pos + left)
                base = seq12_0[left]
            elif left == 0 and EXTRA_PREBASE ==1:
                base = seq12_0[0]
            else: # Variant on First base of chrom (currently only Mitochorndria has first base not an 'N')
                base = ""
            seq1, seq2 = base + seq1, base + seq2
        ret = Variant(self.chrom, self.pos + left, seq1, seq2)
        return ret

    # Aligning variant on the minus strand
    def alignOnMinusStrand(self, reference):
        if self.pos == 1 or (self.chrom not in reference.reflens or not (self.is_insertion() or self.is_deletion())):
            [self.pos, self.ref, self.alt] = self.remove_same_bases_ref_alt()
            return self
        else:
            if self.is_insertion():
                if not len(self.vcf_padded_base) == 0:
                    if self.vcf_padded_base[0] != self.alt[-1]:# See if there is the potential for another copy of the alt allele
                        return self
            elif self.is_deletion():
                if not len(self.vcf_padded_base) == 0: # Skip cases where initial VCF was invalid and did not have pre-base as can happen after splitting multiple alleles
                    if self.vcf_padded_base[0] != self.ref[-1]:
                        return self
        reflen = reference.reflens[self.chrom]

        maxreplen = abs(len(self.ref) - len(self.alt))
        PADLEN = max(100, 1 + 5 * maxreplen)
        if self.pos - PADLEN<1:
            PADLEN = self.pos -1
        # Assume that variant is real .. so that self.pos+len(self.ref) should be valid
        seq1_0 = reference.getReference(self.chrom, self.pos - PADLEN, self.pos + len(self.ref) - 1)
#        s = reference.getReference(self.chrom, self.pos - PADLEN, self.pos - 1)
        s = seq1_0[0:PADLEN]
        N = len(s)
        seq2_0 = s + self.alt
        left, seq1, seq2 = self.leftAlign(seq1_0, seq2_0)
        while left <= maxreplen  and self.pos-PADLEN>1:  # Shifted  all the way to the beginning, make sure there is not more to shift
            PADLEN += PADLEN
            if self.pos-PADLEN<1: # This will only occur for Mitochondria
                PADLEN = self.pos-1
            seq1_0 = reference.getReference(self.chrom, self.pos - PADLEN, self.pos + len(self.ref) - 1)
            #s = reference.getReference(self.chrom, self.pos - PADLEN, self.pos - 1)
            s = seq1_0[0:PADLEN]
            seq2_0 = s + self.alt
            left, seq1, seq2 = self.leftAlign(seq1_0, seq2_0)
            N = len(s)


        if len(seq1) == 0 or len(seq2) == 0:
            # Unless Variant got shifted all the way to position 1(left==0), there should be enough bases for padding.
            if left>0:
                left = left -1
#            base = reference.getReference(self.chrom, self.pos + left - PADLEN, self.pos + left - PADLEN)
                base = seq1_0[left]  # if eliminate 1st common base, take seq1_0[0] ..etc
            else: # left == 0 means no shared bases in chunk, so it means variant fills the chunk .. all the way to the first base
                base = ''
            seq1, seq2 = base + seq1, base + seq2
        else:
            sys.stderr.write("WARNING: Attempted to left shift non-indel variant .. bug? for variant:"+self.id+"\n")
        # Even if no shifting occurs  self.pos has be be shifted by -1
        ret = Variant(self.chrom, self.pos + left - N, seq1, seq2)

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
            if len(s1) == 0 or len(s2) == 0:
                return counter, s1, s2
            if s1[0] != s2[0]:
                return counter, s1, s2
            s1, s2 = s1[1:], s2[1:]
            counter += 1

    # Trimming common ending subsequence of two sequences
    def trimCommonEnd(self, s1, s2):
        counter = 0
        while True:
            if len(s1) == 0 or len(s2) == 0:
                return counter, s1, s2
            if s1[-1] != s2[-1]:
                return counter, s1, s2
            s1, s2 = s1[:-1], s2[:-1]
            counter += 1


# find if the alt-allele (or deletion) contains a repeated pattern (N full repeats )
#
# Can you shift a deletion not multiple of local unit (4 base del in a triplet repeat)? No! .. example
# ACGA[CGAC]GACG -> ACGAGACG
# ACG[ACGA]CGACG -> ACGCGACG  .. not same
# AC[GACG]ACGACG -> ACACGACG   .. not same
# A[CGAC]GACGACG -> AGACGACG   .. not same
# However, you can delete multiple copies of the repeat unit
# ACG[ACGACG]ACG -> ACGACG
# AC[GACGAC]GACG -> ACGACG
# A[CGACGA]CGACG -> ACGACG
# [ACGACG]ACGACG -> ACGACG

# This is why this function tries to split a deleted/inserted sequence into elementary repeat.
def find_repeat_unit(extraseq):
    if len(extraseq) == 0:
        return ["", 0]
    elif len(extraseq) == 1:
        return [extraseq, 1]
    else:
        for irlen in range(1, min(len(extraseq) - 1, int((len(extraseq) + 1) / 2)) + 1):  # max repeat leg
            nsegs = int((len(extraseq)) / irlen)
            seg0 = extraseq[0:irlen]
            allsegsmatch = True
            for iseg in range(1, nsegs):
                segi = extraseq[(iseg * irlen):((iseg + 1) * irlen)]
                if seg0 != segi:
                    allsegsmatch = False
                    break
            if allsegsmatch is True:
                roll_seq = extraseq[(nsegs * irlen):]
                if len(roll_seq) == 0:  # success, found full repeat,
                    return [seg0, nsegs]
        return [extraseq, 1]


# An indel from NGS could be a repeat expansion..
# Priority to the repeat uniq that leads to the most shifting
#   inserted or deleted sequence is mono or exact repeat (no partial)
#   .. if not, perfect repeat
#              shift all the way right, then all the way left.. and
#              scan for smallest possible repeat pattern that covers the interval
#              Once have that pattern, scan for possible expansion of range
def scan_for_repeat(variant, reference):
    if not( variant.is_deletion() or variant.is_insertion()):
        return [None,None,None]
    else:
        if variant.is_insertion():
            rep0 = variant.alt
            pos_left_of_variant = variant.pos - 1  # where do you look for previous repeat
            pos_right_of_variant = variant.pos  # Where do you look for next repeat
        else:
            rep0 = variant.ref
            pos_left_of_variant = variant.pos - 1
            pos_right_of_variant = variant.pos + len(variant.ref)
        [rep, nrep] = find_repeat_unit(rep0)  # The insertion or deletion could be more than 1 repeat unit long.

        lseq = ""
        lrep = len(rep)
        left_context = ""
        nrep_left = 0
        # get reference operations are very expensive, better to get a big chunk.
        left_end = pos_left_of_variant + 1  # Just adding +1 to get loop initialized
        next_left = pos_left_of_variant - len(rep) + 1
        match_rep = True
        while match_rep is True:  # Scan for repeats and load bigger chunks of data as scan toward the end.
            left_right_end = left_end - 1
            left_end = left_end - 40 * lrep
            left_context = reference.getReference(variant.chrom, left_end, left_right_end) + left_context
            while match_rep is True and next_left > left_end + lrep:  # scan for repeat in current chunk
                pseq = left_context[next_left - left_end:(next_left + lrep - left_end)]
                if pseq != rep:
                    match_rep = False
                    lseq = pseq
                else:
                    next_left = next_left - lrep
                    nrep_left += 1

        # Now scan all the way to the right.
        rseq = ""
        nrep_right = 0
        # get reference operations are very expensive, better to get a big chunk.
        right_left_end = pos_right_of_variant  # left end of fasta sequence block
        right_end = right_left_end - 1  # initialize so loop works (will need right_left_end to be variant.pos+1)
        next_right = right_left_end  # start position of next repeat
        match_rep = True
        right_context = ""
        while match_rep is True:
            right_left_end = right_end + 1
            right_end = right_end + 40 * lrep
            right_context = right_context + reference.getReference(variant.chrom, right_left_end,
                                                                   right_end)  # positions are inclusive
            while match_rep is True and next_right + lrep - 1 <= right_end - lrep:
                pseq = right_context[next_right - right_left_end:(next_right + lrep - right_left_end)]
                if pseq != rep:
                    match_rep = False
                    rseq = pseq
                else:
                    next_right = next_right + lrep
                    nrep_right += 1
        # We need to check if changing the repeat unit might  allow a little more shifting
        # 3 cases arise
        # 1) extend on both sides, lead to more repeat units
        # e.g. CTATAT[AT]ATAC -> CTATATATAC was described as an AT deletion, CT[AT](3;2)AC or as a TA deletion C[TA](4;3)C
        # 2)  extend on left only, leading to more left-shifting, but same number of repeat units
        #  CTATAT[AT]ATC -> CTATATATC can be described as AT deletion CT[AT](4;3)C or TA deletion C[TA](4;3)C, the latter being more left-shifted
        # 3) extend to right only, leading to more right-shifting, but same number of repeat units.
        #   CATAT[AT]ATAC -> CATATATAC was described as an AT deletion, C[AT](4;3)AC or as CA[TA](4;3)C
        right_pad = 0
        left_pad = 0
        if lrep > 1 and len(lseq)== lrep and len(rseq) == lrep:
            # see if "rolling" the repeat pattern could extend the repeated section
            # We already know that we can only find a partial pattern on left and/or right

            for npad in range(1, lrep):
                # Check left side extension
                if rep[lrep - npad:] == lseq[lrep - npad:]:
                    left_pad = npad
                # Check right side extension of pattern
                if rep[0:npad] == rseq[0:npad]:
                    right_pad = npad

        extra_rep = 0
        if left_pad + right_pad >= lrep:
            extra_rep = 1
        sys.stdout.write("left_pad="+str(left_pad)+"\n")
        sys.stdout.write("right_pad="+str(right_pad)+"\n")

        # return fully left_padded variant with ref_repeats, alt_repeats, rep_unit, start_pos
        if variant.is_insertion():
            nrep_ref = 0
            nrep_alt = nrep
        else: # Has to be the case that variant.is_deletion() because we excluded other options at beginning
            nrep_ref = nrep
            nrep_alt = 0

        newleft_rep = rep[lrep - left_pad:] + rep[0:lrep - left_pad]
        left_result = [next_left + lrep - left_pad,  # left shifted position of indel for variant obj
                       next_left + lrep - left_pad,  # left shifted position of repeat for HGVS
                       nrep_left + nrep_right + extra_rep + nrep_ref,  # Number of ref repeats
                       nrep_left + nrep_right + extra_rep + nrep_alt,  # number of alt repeats
                       newleft_rep,  # repeat_pattern
                       newleft_rep * nrep_ref,  # trimmed ref-allele
                       newleft_rep * nrep_alt]  # trimmed alt-allele
        # right shifting means two things.
        #    right-shifting insertion or deletion as far right as possible (for vcf-like reporting of variant)
        #    picking the repeat pattern that matches right_pad
        newright_rep = rep[right_pad:] + rep[0:right_pad]
        right_result = [next_right - lrep * nrep + right_pad,  # right shifted position of indel for variant obj
                        next_left + lrep + right_pad,  # right shifted position of repeat for HGVS
                        nrep_left + nrep_right + extra_rep + nrep_ref,  # Number of ref repeats
                        nrep_left + nrep_right + extra_rep + nrep_alt ,  # number of alt repeats
                        newright_rep,  # repeat_pattern
                        newright_rep * nrep_ref,  # trimmed ref-allele
                        newright_rep * nrep_alt]  # trimmed alt-allele
        full_result = [left_result[0],  # left-most position
                       next_right +lrep-1 + right_pad,  # rightmost position, including repeat sequencue
                       left_context[next_left - left_end + lrep - left_pad:] +
                       variant.ref +
                       right_context[:(next_right - right_left_end + right_pad)],
                       left_context[next_left - left_end + lrep - left_pad:] +
                       variant.alt +
                       right_context[:(next_right - right_left_end + right_pad)]
                       ]  # sequence
        return [left_result,right_result,full_result]


def testEqual(predicted, actual):
    if predicted != actual:
        sys.stderr.write("Error predicted:"+predicted+" is not "+actual)


if __name__ == "__main__":
    reference = Reference("/research/bsi/projects/staff_analysis/m037385/savi_data/hg38_nochr.fa")
    seq = "TTTTT"
    [rep, nrep1] = find_repeat_unit(seq)
    testEqual("T", rep)
    testEqual(5, nrep1)

    seq = "TCTCT"
    [rep, nrep2] = find_repeat_unit(seq)
    testEqual("TCTCT", rep)
    testEqual(1, nrep2)

    seq = "TCTCTC"
    [rep, nrep3] = find_repeat_unit(seq)
    testEqual("TC", rep)
    testEqual(3, nrep3)

    seq = "GTCGTCGTC"
    [rep, nrep4] = find_repeat_unit(seq)
    testEqual("GTC", rep)
    testEqual(3, nrep4)


    chrom="13"
    pos = 32329810
    vcf_ref = "ATA"
    vcf_alt = "A"
    sys.stdout.write("VCF:"+chrom+":"+str(pos)+":"+vcf_ref+"/"+vcf_alt+"\n")
    var = Variant(chrom,pos,vcf_ref,vcf_alt)
    sys.stdout.write("left_truth=32329816(vcf_pos=32329815), newl_vcf_ref=A,newl_vcf_alt=ATA\n")
    sys.stdout.write("full_range_of_repeat=32329807-32329818:[TA] repeat\n")

    [left_norm,right_norm, full_norm] = scan_for_repeat(var, reference)
    sys.stdout.write("left_norm =" + ",".join(list(map(str,left_norm)))+ "\n")
    sys.stdout.write("right_norm =" + ",".join(list(map(str,right_norm))) + "\n")
    sys.stdout.write("full_norm =" + ",".join(list(map(str,full_norm))) + "\n")

    sys.stdout.write("right_truth=32329817(vcf_pos=32329816), newr_vcf_ref=ATA,newr_vcf_alt=A\n")
    normr_var = var.alignOnPlusStrand(reference)
    [newposr,ref,alt] = normr_var.get_vcf_pos_ref_alt()
    sys.stdout.write("right_pos = " + str(newposr) + "\n")
    sys.stdout.write("right_ref = " + ref + "\n")
    sys.stdout.write("right_alt = " + alt + "\n")

    sys.stdout.write("left_truth=32329807(vcf_pos=32329806), newl_vcf_ref=TTA,newl_vcf_alt=T\n")
    norml_var = var.alignOnMinusStrand(reference)
    [newposl, ref, alt] = norml_var.get_vcf_pos_ref_alt()
    sys.stdout.write("left_vcf_pos = "+str(newposl)+"\n")
    sys.stdout.write("left_vcf_ref = "+ref+"\n")
    sys.stdout.write("left_vcf_alt = "+alt+"\n")


