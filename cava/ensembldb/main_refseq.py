# This is to import refseq databases.
import datetime
import gzip
import os
import pickle
import sys
import re
from operator import itemgetter
from pathlib import Path

import pybedtools
import requests
from urllib.parse import urlparse

import wget

requests.packages.urllib3.disable_warnings()

import pysam
from cmmodule.utils import read_chain_file
from cmmodule.mapgff import crossmap_gff_file

failed_conversions = dict()
failed_conversions['GENE'] = set()
failed_conversions['GENETYPE'] = set()
failed_conversions['TRANSTYPE'] = set()
failed_conversions['ENST'] = set()


def warn(transcript):
    global failed_conversions
    failed_conversions['GENE'].add(transcript.GENE)
    failed_conversions['GENETYPE'].add(transcript.GENETYPE)
    failed_conversions['TRANSTYPE'].add(transcript.TRANSTYPE)
    failed_conversions['ENST'].add(transcript.ENST)
    # print(f"Messed up: {transcript.ENST} {transcript.GENE}")


def replace_chrom_names(line):
    chrom = line.split('\t')
    if line.startswith('NC_0000'):
        base, v = chrom[0].split('.')
        base = int(base.replace('NC_0000', ''))
        if base == 23:
            base = 'X'
        if base == 24:
            base = 'Y'
        chrom[0] = base
        res = '\t'.join(str(x) for x in chrom)
        return res
    elif line.startswith('NC_012920'):
        _, _ = chrom[0].split('.')
        chrom[0] = 'MT'
        return '\t'.join(str(x) for x in chrom)
    elif line.startswith('chr'):
            base = chrom[0]
            base = base[3:]
            if base == 'M':
                base = 'MT'
            chrom[0] = base
            res = '\t'.join(str(x) for x in chrom)
            return res
    elif chrom[0] in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18',
                      '19', '20', '21', '22', '23', 'MT', 'X', 'Y', 'M']:
        if chrom[0] == 'M':
            chrom[0]= 'MT'
            res = '\t'.join(str(x) for x in chrom)
            return res
        return line
    else:
        return line



# Class representing a transcript
class Transcript(object):

    # Constructor
    def __init__(self):
        self.PROT = None
        self.ENST = None
        self.GENE = None
        self.ENSG = None
        self.CHROM = None
        self.STRAND = None
        self.POS = None
        self.POSEND = None
        self.GENETYPE = None
        self.TRANSTYPE = None
        self.CODING_START = -1
        self.CODING_END = -1
        self.CODING_START_RELATIVE = None
        self.CCDS = None
        self.EXONS = []
        self.PROTL = None
        self.CDNAL = None
        self.isComplete = None

    # Get summary information about the transcript
    def getInfoString(self):
        if self.STRAND == '1':
            ret = '+/'
        else:
            ret = '-/'
        cdna = self.getcDNALength()
        return ret + str(self.POSEND - self.POS + 1) + 'bp/' + str(len(self.EXONS)) + '/' + str(cdna) + 'bp/' + str(
            self.getProteinLength())

    # Get cDNA length of the transcript
    def getcDNALength(self):
        ret = 0
        for exon in self.EXONS:
            ret += exon.END - exon.START
        return ret

    # Get protein length of the transcript
    def getProteinLength(self):
        codingdna = 0
        if self.STRAND == '1':
            for exon in self.EXONS:
                if exon.END < self.CODING_START: continue
                if exon.START > self.CODING_END: continue
                if exon.START <= self.CODING_START <= exon.END:
                    start = self.CODING_START
                else:
                    start = exon.START + 1
                if exon.START <= self.CODING_END <= exon.END:
                    end = self.CODING_END
                else:
                    end = exon.END
                codingdna += end - start + 1
        else:
            for exon in self.EXONS:
                if exon.START > self.CODING_START: continue
                if exon.END < self.CODING_END: continue
                if exon.START <= self.CODING_START <= exon.END:
                    end = self.CODING_START
                else:
                    end = exon.END
                if exon.START <= self.CODING_END <= exon.END:
                    start = self.CODING_END
                else:
                    start = exon.START + 1
                codingdna += end - start + 1
        return int((codingdna - 3) / 3)

    # Check if it is a candidate transcript
    def isCandidate(self):
        return self.CODING_START > -1 and self.CODING_END > -1
        # if not (self.GENETYPE == 'protein_coding' and self.TRANSTYPE == 'protein_coding'):
        #    return False
        # return (self.CODING_START > -1 and self.CODING_END > -1) and self.isComplete

    # Output transcript
    def output(self, outfile, outfile_list):
        self.ENSG = self.GENE
        out = self.ENST + '\t' + self.GENE + '\t' + self.ENSG + '\t' + self.getInfoString() + '\t' + \
              self.CHROM + '\t' + self.STRAND + '\t' + str(self.POS)
        out += '\t' + str(self.POSEND) + '\t' + str(self.CODING_START_RELATIVE) + '\t' + str(self.CODING_START)
        out += '\t' + str(self.CODING_END)
        for exondata in self.EXONS: out += '\t' + str(exondata.START) + '\t' + str(exondata.END)
        outfile.write(out + '\n')

        outfile_list.write(self.ENSG + '\t' + self.GENE + '\t' + self.ENST + '\t' + self.PROT + '\n')

    # Finalize transcript
    def finalize(self):
        if self.STRAND == '1':
            self.POS = self.EXONS[0].START
            self.POSEND = self.EXONS[len(self.EXONS) - 1].END
            assert self.POS <= self.POSEND
            codingStartRelative = 0
            for exondata in self.EXONS:
                if exondata.START <= self.CODING_START <= exondata.END:
                    codingStartRelative += self.CODING_START - exondata.START
                    break
                else:
                    codingStartRelative += exondata.END - exondata.START
            self.CODING_START_RELATIVE = codingStartRelative
        else:
            self.POS = self.EXONS[len(self.EXONS) - 1].START
            self.POSEND = self.EXONS[0].END

            codingStartRelative = 0
            for exondata in self.EXONS:
                if exondata.START <= self.CODING_START <= exondata.END:
                    codingStartRelative += exondata.END - self.CODING_START + 1
                    break
                else:
                    codingStartRelative += exondata.END - exondata.START
            self.CODING_START_RELATIVE = codingStartRelative
        self.PROTL = self.getProteinLength()
        self.CDNAL = self.getcDNALength()


# Class representing an exon
class Exon(object):

    # Constructor
    def __init__(self, start, end):
        self.START = start
        self.END = end


# Class representing a gene
class Gene(object):

    # Constructor
    def __init__(self, symbol, ensg):
        self.SYMBOL = symbol
        self.ENSG = ensg
        self.TRANSCRIPTS = dict()

    # Select ICR transcript
    def selectTranscript(self):
        ccds_set = []
        nonccds_set = []
        for enst, transcript in self.TRANSCRIPTS.items():
            if transcript.CCDS:
                ccds_set.append(transcript)
            else:
                nonccds_set.append(transcript)

        if len(ccds_set) > 0:
            candidates = ccds_set
        else:
            candidates = nonccds_set

        # Sort candidates by ENST. In case there multiple selectable
        # transcripts, the selected one does not depend on the order they are
        # returned by .items()
        candidates.sort(key=lambda x: x.ENST)

        selected = Transcript()
        selected.PROTL = selected.CDNAL = -1
        for t in candidates:
            # Note that we return the *last* selected candidate since we
            # overwrite the variable `selected`.
            if t.PROTL > selected.PROTL:
                selected = t
            elif t.PROTL == selected.PROTL and t.CDNAL > selected.CDNAL:
                selected = t

        return selected

    # Output all or selected transcripts
    def output(self, outfile, outfile_list, select, target_transcripts):
        for t, transcript in self.TRANSCRIPTS.items():
            if select:
                if t in target_transcripts:
                    transcript.output(outfile, outfile_list)
            else:
                try:
                    transcript.output(outfile, outfile_list)
                except:
                    warn(transcript)


#######################################################################################################################
def write_temp(output_name, options, candidates, genesdata):
    outfile = open('temp.txt', 'w')

    # Initialize output list file if needed
    outfile_list = open(output_name, 'w')

#    outfile_list.write(
#        '# Created by CAVA refseq_db ' + options.version + ' based on refseq release ' + options.refseq + '\n')
    outfile_list.write(
        '# Created by CAVA refseq_db based on ' + options.url_gtf + '\n')
    outfile_list.write('#GENE1\tGENE2\tTRANSCRIPT\tPROTEIN\n')

    # Output transcripts of each gene
    for ensg, gene in genesdata.items():
        try:
            gene.output(outfile, outfile_list, options.select, candidates)
        except:
            print(f'Failed {gene.SYMBOL}, {gene.TRANSCRIPTS}')
    # Close temporary output files
    outfile.close()
    outfile_list.close()


def build_tx_to_prot_dict(opener, filename):
    print(f'\nBuilding transcript to protein mapping')
    tx_to_prot_dict = dict()
    for line in opener(filename, 'rt'):
        line = line.strip()
        if line.startswith('#'): continue
        cols = line.split('\t')
        tags = cols[8].split(';')
        enst_prot = getValue(tags, 'transcript_id')
        prot = getValue(tags, 'protein_id')
        if prot is not None:
            tx_to_prot_dict[enst_prot] = prot
    return tx_to_prot_dict


def parse_GTF(filename='', options=None, genesdata=None, transIDs=None):
    first = True
    prevenst = ''
    transcript = None

    if filename.endswith('gz'):
        opener = gzip.open
    else:
        opener = open

    tx_to_prot_dict = build_tx_to_prot_dict(opener, filename)

    print(f'Parsing {filename}', end="...")

    for line in opener(filename, 'rt'):
        line = line.strip()
        if line.startswith('#'): continue
        cols = line.split('\t')

        # Consider only certain types of lines
        if cols[2] not in ['exon', 'transcript', 'start_codon', 'stop_codon']:
            continue
        if cols[0] not in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16',
                           '17',
                           '18', '19', '20', '21', '22', '23', 'MT', 'X', 'Y', 'M'] and \
                cols[0] not in ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
                                'chr11',
                                'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
                                'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chr23', 'chrMT', 'chrX', 'chrY'] and \
                not cols[0].startswith('NC_0'):
            continue
        if cols[0].startswith("NC_0"):
            pmatch = re.match('NC_[0]+([1-9]+[0-9]*)\.', cols[0])
            if pmatch is None:
                continue
            id = pmatch.group(1)
            if id not in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16',
                          '17', '18',
                          '19', '20', '21', '22', '23', '24', '12920']:
                continue


        # Annotation tags
        tags = cols[8].split(';')
        # Retrieve transcript ID
        enst = getValue(tags, 'transcript_id')

        # Do not consider transcript if it is not on the custom transcript list
        if options.input is not None and enst not in transIDs: continue
        # Exclude non-NM transcripts
        if options.nm_only is not None and not enst.startswith('NM'):
            continue

        # Finalize and output transcript object
        if not enst == prevenst:

            # Finalize transcript and add to Gene object if candidate
            if not first:
                try:
                    transcript.finalize()
                except:
                    warn(transcript)

                if transcript.isCandidate():
                    if transcript.ENSG not in list(genesdata.keys()):
                        genesdata[transcript.ENSG] = Gene(transcript.GENE, transcript.ENSG)
                    genesdata[transcript.ENSG].TRANSCRIPTS[transcript.ENST] = transcript

            # Initialize new Transcript object
            transcript = Transcript()
            transcript.ENST = enst
            transcript.GENE = getValue(tags, 'gene')
            transcript.ENSG = getValue(tags, 'gene_id')
            try:
                transcript.PROT = tx_to_prot_dict[enst]
            except KeyError:
                print(f'refseq:{enst} not in protein database. No protein_id tag in gtf ')
                transcript.PROT = ''
                transcript.PROTL = 0

            transcript.CHROM = cols[0]
            if cols[6] == '+':
                transcript.STRAND = '1'
            else:
                transcript.STRAND = '-1'

        # If line represents an exon
        if cols[2] == 'exon':
            idx = 0
            for x in tags:
                x = x.strip()
                if x.startswith('exon_number'):
                    s = x[x.find('\"') + 1:]
                    idx = int(s[:s.find('\"')]) - 1
                    break
            start = int(cols[3]) - 1
            end = int(cols[4])
            if idx >= len(transcript.EXONS):
                for _ in range(len(transcript.EXONS), idx + 1): transcript.EXONS.append(None)
            transcript.EXONS[idx] = Exon(start, end)

        if cols[2] == 'start_codon':
            if transcript.STRAND == '1':
                if transcript.CODING_START < 0 or int(cols[3]) < transcript.CODING_START: transcript.CODING_START = int(
                    cols[3])
            else:
                if transcript.CODING_START < 0 or int(cols[4]) > transcript.CODING_START: transcript.CODING_START = int(
                    cols[4])

        if cols[2] == 'stop_codon':
            if transcript.STRAND == '1':
                if transcript.CODING_END < 0 or int(cols[4]) > transcript.CODING_END: transcript.CODING_END = int(
                    cols[4])
            else:
                if transcript.CODING_END < 0 or int(cols[3]) < transcript.CODING_END: transcript.CODING_END = int(
                    cols[3])

        # Check if transcript is complete and is a CCDS transcript
        if transcript.isComplete is None:
            transcript.isComplete = not (
                    getBooleanValue(tags, 'cds_start_NF') or getBooleanValue(tags, 'cds_end_NF'))
            if getValue(tags, 'ccds_id') is not None:
                transcript.CCDS = True
            else:
                transcript.CCDS = False

        prevenst = enst
        if first: first = False
    return transcript, prevenst, first, genesdata


def sort_tmpfile(f):
    # Sort temporary output file
    data = dict()
    counter = 0
    for line in open(f, 'r'):
# NM-only should already have been filtered out .. if specified as an option
# if not line.startswith('N'): continue
        counter += 1
        line.rstrip()
        record = line.split('\t')
        record[6] = int(record[6])
        chrom = record[4]

        if record[4] in list(data.keys()):
            data[record[4]].append(record)
        else:
            data[record[4]] = []
            data[record[4]].append(record)

    sys.stdout.write('OK\n')
    sys.stdout.write(f'Sorting {counter} transcripts... ')
    sys.stdout.flush()
    sortedRecords = sortRecords(data, 6, 7)
    return sortedRecords


# Retrieve tag value
def getValue(tags, tag):
    ret = None
    for x in tags:
        x = x.strip()
        if x.startswith(tag):
            s = x[x.find('\"') + 1:]
            ret = s[:s.find('\"')]
            break
    return ret


# Retrieve boolean tag value
def getBooleanValue(tags, tag):
    for x in tags:
        x = x.strip()
        if x.startswith('tag'):
            s = x[x.find('\"') + 1:]
            value = s[:s.find('\"')]
            if value == tag: return True
    return False


# Read transcript IDs from file
def readTranscriptIDs(inputfn):
    ret = set()
    for line in open(inputfn): ret.add(line.strip())
    return ret


# Sort records in file
def sortRecords(records, idx1, idx2):
    ret = []
    chroms = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19',
              '20', '21', '22', '23', 'X', 'Y','MT']
    allkeys = list(records.keys())
    ordered_chroms = chroms
    ordered_chroms.sort()
    for key in allkeys:
        if key not in ordered_chroms:
            ordered_chroms.append(key)
    for i in range(len(ordered_chroms)):
        chrom = ordered_chroms[i]
        if chrom in list(records.keys()):
            records[chrom] = sorted(records[chrom], key=itemgetter(idx1, idx2))
    for i in range(len(ordered_chroms)):
        chrom = ordered_chroms[i]
        if chrom in list(records.keys()):
            for record in records[chrom]: ret.append(record)
    return ret


# Write records to file
def writeToFile(sortedRecords, filename):
    outfile = open(filename, 'w')
    for record in sortedRecords:
        s = str(record[0]).rstrip()
        for i in range(1, len(record)): s += '\t' + str(record[i]).rstrip()
        outfile.write(s + '\n')
    outfile.close()


# Read records from file as a list
def readRecords(inputfn):
    ret = []
    for line in open(inputfn): ret.append(line.strip())
    return ret


# Process Ensembl data
def process_data(options):
    # Dictionary of Gene objects
    genesdata = dict()

    # Load custom transcript IDs
    transIDs = None
    if options.input is not None:
        transIDs = readTranscriptIDs(options.input)
        print('\nOnly ' + str(len(transIDs)) + ' transcripts read from ' + options.input + ' are considered\n')
    else:
        nm = 'All transcripts from the release are considered'
        if options.nm_only: nm = 'All NM transcripts from the release are considered'
        print(f'\n{nm}\n')

    # Load candidate and CCDS data for Ensembl <75
    dict()

    ######################################################################
    # Download RefSeq data if necessary
    url_compressed_gtf = options.url_gtf
    pu = urlparse(url_compressed_gtf)
    path = pu.path
    paths = path.split('/')
    source_compressed_gtf = paths[len(paths)-1]
    fname = source_compressed_gtf

    #source_compressed_gtf = options.refseq + '_genomic.gtf.gz'
    # https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/reference/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gtf.gz
    source_compressed_gtf = os.path.join(options.output_dir, source_compressed_gtf)

    if not os.path.exists(source_compressed_gtf):
        sys.stdout.write('Downloading RefSeq database... ')
        sys.stdout.flush()

        #url = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/reference/' + options.refseq + '/' + options.refseq + '_genomic.gtf.gz'
        url = options.url_gtf
        try:
            wget.download(url)
            sys.stdout.flush()
            # Convert chromosome names #Note we will lose unmapped transcripts here!
            # fname = options.refseq + "_genomic.gtf.gz"
            print(f'\nUnzipping {fname}')
            #cmd = 'bgzip -d ' + options.refseq + '_genomic.gtf.gz'
            cmd = 'bgzip -d ' + fname
            os.system(cmd)
            #gtf_fname = options.refseq + '_genomic.gtf'
            gtf_fname = Path(fname).stem
            out = open('temp.txt', 'w')
            print(f'Parsing {gtf_fname}')
            with open(gtf_fname, 'r') as g:
                for line in g:
                    if line.startswith('#'): continue
                    try:
                        new_line = replace_chrom_names(line)
                        line = new_line
                    except:
                        print(f'Failed: {line}')
                        exit()
                    if line is not None:
                        out.write(new_line)
            out.close()
            print(f'Compressing the GTF into: {fname}')
            source_compressed_gtf = os.path.join(options.output_dir,fname)
            cmd = 'bgzip -c temp.txt > ' + source_compressed_gtf
            os.system(cmd)
            os.remove('temp.txt')
            os.remove(gtf_fname)
        except Exception as e:
            print('\n\nCannot connect to RefSeq FTP site. No internet connection?\n')
            print(f'{e}\n{url}')
            quit()
    ################################################################
    # Use crossmap to get hg19 if desired
    #################################################################
    converted_gtf = None
    if options.build == 'GRCh38' and options.no_hg19 is not False:
        requests.packages.urllib3.util.ssl_.DEFAULT_CIPHERS += 'HIGH:!DH:!aNULL'  # Needed for UCSC
        # only download if necessary
        if not os.path.exists(os.path.join('data', 'hg38ToHg19.over.chain.gz')):
            sys.stdout.write('Downloading UCSC database... ')
            sys.stdout.flush()
            url = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz'
            try:
                p = requests.get(url, verify=False)
                with open(os.path.join('data', 'hg38ToHg19.over.chain.gz'), 'wb') as o:
                    o.write(p.content)

            except Exception as e:
                print('\n\nCannot connect to UCSC FTP site. No internet connection?\n')
                print(f'Exception: {e}')
                quit()
        if fname.endswith(".gtf.gz"):
            fname_path = os.path.splitext(fname)
            fname_path = os.path.splitext(fname_path[0])
            basename = fname_path[0]
        elif fname.endswith(".gtf"):
            fname_path = os.path.splitext(fname)
            basename = fname_path[0]
        else:
            basename = fname
#        converted_gtf = os.path.join(options.output_dir, 'Homo_sapiens.RefSeq.hg19_converted.' + options.refseq + '.gtf')
        converted_gtf = os.path.join(options.output_dir,
                                     'Homo_sapiens.RefSeq.hg19_converted.' + basename + '.gtf')

        if not os.path.exists(converted_gtf):
            sys.stdout.write('\nMaking a hg19-converted GTF file\n')
            mapTree, targetChromSizes, sourceChromSizes = read_chain_file(
                os.path.join('data', 'hg38ToHg19.over.chain.gz'))
            crossmap_gff_file(mapTree, source_compressed_gtf, converted_gtf)

            # Note this file is not sorted!
            a = pybedtools.BedTool(converted_gtf)
            a.sort().remove_invalid().saveas('tmp.txt')
            os.rename('tmp.txt', converted_gtf)

    ################################################################
    #
    #################################################################
    # Iterate through the lines in the refseq data file
    sys.stdout.write('Extracting transcript data from RefSeq...')

    transcript, prevenst, first, genesdata = parse_GTF(filename=source_compressed_gtf,
                                                       options=options,
                                                       genesdata=genesdata,
                                                       transIDs=transIDs)

    sys.stdout.write('Done\n')
    sys.stdout.flush()
    # Finalize last transcript and add to Gene object if candidate
    if transcript is not None:
        transcript.finalize()
        if transcript.isCandidate():
            if transcript.ENSG not in list(genesdata.keys()): genesdata[transcript.ENSG] = Gene(transcript.GENE,
                                                                                                transcript.ENSG)
            genesdata[transcript.ENSG].TRANSCRIPTS[transcript.ENST] = transcript

    # If no transcript ID from the input file was found in the Ensembl release
    if len(genesdata) == 0:
        print('\n\nNo transcripts found in this release.')
        print('\nNo transcript database created.')
        print("-----------------------------------------------------------------\n")
        quit()

    write_temp(os.path.join(options.output_dir, options.output + '.txt'), options, transIDs, genesdata)
    enst_records = sort_tmpfile('temp.txt')
    assert (len(enst_records) > 0)
    writeToFile(enst_records, os.path.join(options.output_dir, options.output))

    failed_conversions['GENE'] = set()
    failed_conversions['GENETYPE'] = set()
    failed_conversions['TRANSTYPE'] = set()
    failed_conversions['ENST'] = set()
    # ################################################################
    # Begin converted GTF conversion
    # ################################################################
    hg19_records = []
    if options.build == 'GRCh38' and options.no_hg19 is not False:
        sys.stdout.write('Extracting transcript data from hg19 version...')
        sys.stdout.flush()
        transcript, prevenst, first, genesdata = parse_GTF(filename=converted_gtf,
                                                           options=options,
                                                           genesdata=genesdata,
                                                           transIDs=transIDs)

        # Finalize last transcript and add to Gene object if candidate
        if transcript is not None:
            try:
                transcript.finalize()
            except:
                warn(transcript)

            if transcript.isCandidate():
                if transcript.ENSG not in list(genesdata.keys()): genesdata[transcript.ENSG] = Gene(transcript.GENE,
                                                                                                    transcript.ENSG)
                genesdata[transcript.ENSG].TRANSCRIPTS[transcript.ENST] = transcript

        # If no transcript ID from the input file was found in the release
        if len(genesdata) == 0:
            print('\n\nNo transcripts from ' + options.input + ' found in the release.')
            print('\nNo transcript database created.')
            print("-----------------------------------------------------------------\n")
            quit()

        write_temp(os.path.join(options.output_dir, options.output + '.hg19_converted.txt'), options, transIDs,
                   genesdata)
        sortedRecords = sort_tmpfile('temp.txt')
        writeToFile(sortedRecords, os.path.join(options.output_dir, options.output + '.hg19_converted'))
        sys.stdout.write('Completed hg19 version...')
        sys.stdout.flush()
        pickle.dump(failed_conversions,
                    open(os.path.join(options.output_dir, options.output + '_failed_conversions.pkl'), 'wb'))
        hg19_records = sortedRecords
    # ################################################################
    # END converted GTF conversion
    # ################################################################

    # Remove temporary files
    sys.stdout.write('OK\n')
    sys.stdout.write('Removing temporary files... ')
    sys.stdout.flush()
    os.remove('temp.txt')
    os.remove(source_compressed_gtf)

    print(f"Failed {failed_conversions['GENE'].__len__()} Genes and {failed_conversions['ENST'].__len__()} transcripts")

    # Return sorted records
    return len(enst_records), len(hg19_records)


# Use Tabix to index output file
def indexFile(f, options):
    sys.stdout.write(f'Compressing output file {f}... ')
    sys.stdout.flush()
    pysam.tabix_compress(os.path.join(options.output_dir, f), os.path.join(options.output_dir, f + '.gz'), force=True)
    sys.stdout.write('OK\n')
    sys.stdout.write(f'Indexing output file {f}... ')
    sys.stdout.flush()
    pysam.tabix_index(os.path.join(options.output_dir, f + '.gz'), seq_col=4, start_col=6, end_col=7, meta_char='#',
                      force=True)
    sys.stdout.write('OK\n')


# CHeck if string is a number (integer)
def is_number(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def run(options):
    # Checking if all required options specified

    if options.output is None:
        print('\nError: no output file name specified. Use option -h to get help!\n')
        quit()

    # Genome build
    # genome_build = options.genome
    genome_build = 'GRCh37' if 'GRCh37' in options.build else 'GRCh38'

    # Printing out version.py information
    print("\n---------------------------------------------------------------------------------------")
    print('CAVA ' + options.version + ' transcript database preparation tool (refseq_db) is now running')
    print('Started: ', datetime.datetime.now(), '\n')

    # Print info
    #print('RefSeq version:  ' + options.refseq)
    print('Reference genome: ' + genome_build)

    # Creating compressed output file
    enst_parsed, ens_lifted = process_data(options)
    print('\nA total of ' + str(enst_parsed) + ' transcripts have been retrieved\n')
    # Indexing output file with Tabix
    indexFile(options.output, options)
    # Printing out summary information
    print('')
    print('---------------------')
    print('Output files created:')
    print('---------------------')
    print(options.output + '.gz (transcript database)')
    print(options.output + '.gz.tbi (index file)')
    print(options.output + '.txt (list of transcripts)')

    if ens_lifted:
        print('\nA total of ' + str(ens_lifted) + ' transcripts have been lifted over\n')
        # Indexing output file with Tabix
        indexFile(options.output + '.hg19_converted', options)
        # Printing out summary information
        print('')
        print('---------------------')
        print('Output files created:')
        print('---------------------')
        print(options.output + '.hg19_converted' + '.gz (transcript database)')
        print(options.output + '.hg19_converted' + '.gz.tbi (index file)')
        print(options.output + '.hg19_converted' + '.txt (list of transcripts)')
        os.remove(os.path.join(options.output_dir, options.output + '.hg19_converted'))

    # Removing uncompressed output file
    os.remove(os.path.join(options.output_dir, options.output))

    print('')
    print('CAVA refseq_db successfully finished: ', datetime.datetime.now())
    print("---------------------------------------------------------------------------------------\n")
