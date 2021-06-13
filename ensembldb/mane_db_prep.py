import datetime
import gzip
import os
import sys
from operator import itemgetter

import pybedtools
import requests
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
    # print(f'FAILED: {failed_conversions["GENE"]}, failed_conversions['ENST']')
    # raise Exception(f"Messed up: {transcript.GENE}")


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
            if not 'START' in dir(exon) or not 'END' in dir(exon): continue
            ret += exon.END - exon.START
        return ret

    # Get protein length of the transcript
    def getProteinLength(self):
        codingdna = 0
        if self.STRAND == '1':
            for exon in self.EXONS:
                if not 'START' in dir(exon) or not 'END' in dir(exon): continue
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
                if not 'START' in dir(exon) or not 'END' in dir(exon): continue
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
        if not (self.GENETYPE == 'protein_coding' and self.TRANSTYPE == 'protein_coding'):
            return False
        # return (self.CODING_START is not None and self.CODING_END is not None) and self.isComplete
        return (self.CODING_START > -1 and self.CODING_END > -1) and self.isComplete

    # Output transcript
    def output(self, outfile, outfile_list):
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
                if not 'START' in dir(exondata) or not 'END' in dir(exondata): continue
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
                if not 'START' in dir(exondata) or not 'END' in dir(exondata): continue
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
        assert start > 0
        assert end > 0
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
                transcript.output(outfile, outfile_list)


#######################################################################################################################
def write_temp(output_name, options, candidates, genesdata):
    outfile = open('temp.txt', 'w')

    # Initialize output list file if needed
    outfile_list = open(output_name, 'w')

    outfile_list.write(
        '# Created by CAVA' + options.version + ' based on MANE release ' + options.ensembl + '\n')
    outfile_list.write('GENEID\tSYMBOL\tTranscript\tProtein\n')

    # Output transcripts of each gene
    for ensg, gene in genesdata.items():
        try:
            gene.output(outfile, outfile_list, options.select, candidates)
        except:
            print(f'Write_temp failed {gene.SYMBOL}, {gene.TRANSCRIPTS}')

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


def parse_GTF(filename='', genesdata=None):
    first = True
    prevenst = ''
    transcript = None

    if filename.endswith('gz'):
        opener = gzip.open
    else:
        opener = open
    print(f'Parsing {filename}', end="...")
    tx_to_prot_dict = build_tx_to_prot_dict(opener, filename)

    for line in opener(filename, 'rt'):

        line = line.strip().replace('chr', '')
        if line.startswith('#'): continue
        cols = line.split('\t')

        # Only consider transcripts on the following chromosomes
        if cols[0] not in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17',
                           '18', '19', '20', '21', '22', '23', 'MT', 'X', 'Y']: continue

        # Consider only certain types of lines
        if cols[2] not in ['exon', 'transcript', 'start_codon', 'stop_codon']: continue

        # Annotation tags
        tags = cols[8].split(';')
        # Retrieve transcript ID
        enst = getValue(tags, 'transcript_id')

        # Finalize and output transcript object
        if not enst == prevenst:

            # Finalize transcript and add to Gene object if candidate
            if transcript is not None:
                genesdata = finalize_last_tx(transcript, genesdata)

            # Initialize new Transcript object
            transcript = Transcript()
            transcript.ENST = enst
            transcript.GENE = getValue(tags, 'gene_id')
            transcript.ENSG = getValue(tags, 'gene')
            try:
                transcript.PROT = tx_to_prot_dict[enst]
            except KeyError:
                print(f'enst {enst} not in database ')
                transcript.PROT = ''

            transcript.CHROM = cols[0]
            if cols[6] == '+':
                transcript.STRAND = '1'
            else:
                transcript.STRAND = '-1'

            # Retrieve gene biotype and transcript biotype
            transcript.GENETYPE = getValue(tags, 'gene_type')
            if transcript.GENETYPE is None: transcript.GENETYPE = getValue(tags, 'gene_biotype')
            transcript.TRANSTYPE = getValue(tags, 'transcript_type')
            if transcript.TRANSTYPE is None: transcript.TRANSTYPE = getValue(tags, 'transcript_biotype')
            if transcript.TRANSTYPE is None: transcript.TRANSTYPE = cols[1]

        # If line represents an exon
        style = None
        if cols[2] == 'exon':
            idx = 0
            for x in tags:
                x = x.strip()
                if x.startswith('exon_number'):
                    try:
                        # Ensemble Entry
                        idx = int(x.split()[1])
                        style = 'ensembl'
                        transcript.GENE = getValue(tags, 'gene_name')

                    except ValueError:
                        # RefSeq Entry
                        s = x[x.find('\"') + 1:]
                        idx = int(s[:s.find('\"')]) - 1
                        style = 'refseq'
                    break
            start = int(cols[3]) - 1
            end = int(cols[4])
            assert type(start) is int
            assert type(end) is int

            if idx >= len(transcript.EXONS):
                for _ in range(len(transcript.EXONS), idx): transcript.EXONS.append(None)
            if style == 'ensembl':
                if idx >= len(transcript.EXONS):
                    for _ in range(len(transcript.EXONS), idx): transcript.EXONS.append(None)
                transcript.EXONS[idx - 1] = Exon(start, end)
            else:
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
            transcript.CCDS = True

        prevenst = enst
        if first: first = False
    return transcript, prevenst, first, genesdata


def sort_tmpfile(f):
    # Sort temporary output file
    data = dict()
    counter = 0
    for line in open(f, 'r'):
        if not line.startswith('ENST') and not line.startswith('NM_'): continue
        counter += 1
        line.rstrip()
        record = line.split('\t')
        record[6] = int(record[6])
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
              '20', '21', '22', '23', 'MT', 'X', 'Y']
    for i in range(len(chroms)):
        chrom = chroms[i]
        if chrom in list(records.keys()):
            records[chrom] = sorted(records[chrom], key=itemgetter(idx1, idx2))
    for i in range(len(chroms)):
        chrom = chroms[i]
        if chrom in list(records.keys()):
            for record in records[chrom]: ret.append(record)
    return ret


# Write records to file
def writeToFile(sortedRecords, filename):
    print(f"writeToFile to {filename}")
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
    0
    0
    ref_records_hg19 = 0
    enst_records_hg19 = 0
    ######################################################################
    # Download data if necessary
    source_compressed_gtf_ens = 'MANE.GRCh38.v' + options.ensembl + '.select_ensembl_genomic.gtf.gz'
    source_compressed_gtf_ref = 'MANE.GRCh38.v' + options.ensembl + '.select_refseq_genomic.gtf.gz'
    source_compressed_gtf_ens = os.path.join('data', source_compressed_gtf_ens)
    source_compressed_gtf_ref = os.path.join('data', source_compressed_gtf_ref)
    download_gtf(source_compressed_gtf_ens, options.ensembl)
    download_gtf(source_compressed_gtf_ref, options.ensembl)

    # Parse the entries & return the number of transcripts found
    genesdata = dict()
    transIDs = None
    enst_records = parse_gtf_loop(source_compressed_gtf_ens, options, genesdata, transIDs)
    print(f'Completed {enst_records} ENSEMBL records')

    genesdata = dict()
    transIDs = None
    ref_records = parse_gtf_loop(source_compressed_gtf_ref, options, genesdata, transIDs)
    print(f'Completed {ref_records} RefSeq records')
    ################################################################
    # Use crossmap to get hg19 if desired
    #################################################################
    if options.no_hg19 is not False:
        crossmap(source_compressed_gtf_ens)
        crossmap(source_compressed_gtf_ref)

        failed_conversions['GENE'] = set()
        failed_conversions['GENETYPE'] = set()
        failed_conversions['TRANSTYPE'] = set()
        failed_conversions['ENST'] = set()
        # ################################################################
        # Begin converted GTF conversion
        # ################################################################
        sys.stdout.write('Extracting transcript data for hg19 version...\n')
        sys.stdout.flush()
        genesdata = dict()
        transIDs = None

        enst_records_hg19 = parse_gtf_loop(source_compressed_gtf_ens.replace('.gtf.gz', '.hg19_converted.gtf'), options,
                                           genesdata, transIDs)
        print(f'Completed {enst_records_hg19} hg19-converted ENSEMBL records')

        assert enst_records_hg19 > 0, "Uh oh.  I didn't get any ensemble transcripts converted. Check your inputs"
        genesdata = dict()
        transIDs = None

        ref_records_hg19 = parse_gtf_loop(source_compressed_gtf_ref.replace('.gtf.gz', '.hg19_converted.gtf'), options,
                                          genesdata, transIDs)
        print(f'Completed {ref_records_hg19} hg19-converted RefSeq records')
        assert ref_records_hg19 > 0, "Uh oh.  I didn't get any refseq transcripts converted. Check your inputs"

    # ################################################################
    # END converted GTF conversion
    # ################################################################

    # Remove temporary files
    sys.stdout.write('OK\n')
    sys.stdout.write('Removing temporary files... ')
    sys.stdout.flush()
    os.remove('temp.txt')
    # os.remove(source_compressed_gtf)

    print(f"Failed {failed_conversions['GENE'].__len__()} Genes and {failed_conversions['ENST'].__len__()} transcripts")

    # Return sorted records
    return enst_records, ref_records, ref_records_hg19, enst_records_hg19


# Use Tabix to index output file
def indexFile(f):
    sys.stdout.write(f'Compressing output file {f}... ')
    sys.stdout.flush()
    assert os.path.exists(f), f"{f} does not exist"
    pysam.tabix_compress(f, f + '.gz', force=True)
    sys.stdout.write('OK\n')
    if os.path.exists(f):
        os.remove(f)
    sys.stdout.write(f'Indexing output file {f}.gz... ')
    sys.stdout.flush()
    pysam.tabix_index(f + '.gz', seq_col=4, start_col=6, end_col=7, meta_char='#', force=True)
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
    if options.ensembl is None:
        print('\nError: no release specified. Use option -h to get help!\n')
        quit()

    # Genome build
    # genome_build = options.genome
    genome_build = 'GRCh38'

    # Printing out version.py information
    print("\n---------------------------------------------------------------------------------------")
    print('CAVA ' + options.version + ' transcript database preparation tool (ensembl_db) is now running')
    print('Started: ', datetime.datetime.now(), '\n')

    # Print info
    print('Version:  ' + options.ensembl)
    print('Reference genome: ' + genome_build)

    # Creating compressed output file
    enst_parsed, ref_parsed, ens_lifted, ref_lifted = process_data(options)

    report_summary(enst_parsed, options, hg19=False, enst=True)
    report_summary(ref_parsed, options, hg19=False, enst=False)
    report_summary(ens_lifted, options, hg19=True, enst=True)
    report_summary(ref_lifted, options, hg19=True, enst=False)

    print('')
    print('CAVA ensembl_db successfully finished: ', datetime.datetime.now())
    print("---------------------------------------------------------------------------------------\n")


def download_gtf(source_compressed_gtf, version):
    if not os.path.exists(source_compressed_gtf):
        sys.stdout.write(f'Downloading {os.path.basename(source_compressed_gtf)}... ')
        sys.stdout.flush()
        if 'ensembl' in source_compressed_gtf:
            url = 'ftp://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_' + version + '/MANE.GRCh38.v' + version + '.select_ensembl_genomic.gtf.gz'
        else:
            url = 'ftp://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_' + version + '/MANE.GRCh38.v' + version + '.select_refseq_genomic.gtf.gz'
        try:
            wget.download(url)
            os.rename(os.path.basename(source_compressed_gtf), source_compressed_gtf)
        except Exception as e:
            print('\n\nCannot connect to FTP site. No internet connection?\n')
            print(f'{e}\n{url}')
            quit()
    print('')
    sys.stdout.flush()


def crossmap(source_compressed_gtf):
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

    sys.stdout.write('\nMaking a hg19-conveterted GTF file\n')
    mapTree, targetChromSizes, sourceChromSizes = read_chain_file(os.path.join('data', 'hg38ToHg19.over.chain.gz'))

    converted_gtf = source_compressed_gtf.replace('.gtf.gz', '.hg19_converted.gtf')
    crossmap_gff_file(mapTree, source_compressed_gtf, converted_gtf)

    # Note this file is not sorted!
    a = pybedtools.BedTool(converted_gtf)
    a.sort().remove_invalid().saveas('tmp.txt')
    os.rename('tmp.txt', converted_gtf)


def finalize_last_tx(transcript, genesdata):
    # Finalize last transcript and add to Gene object if candidate
    if transcript is not None:
        try:
            transcript.finalize()
        except:
            warn(transcript)
            print(f'WARNING: Trancript failed: {transcript.GENE}, {transcript.ENST}, {transcript.ENSG}')

        if transcript.ENSG not in list(genesdata.keys()): genesdata[transcript.ENSG] = Gene(transcript.GENE,
                                                                                            transcript.ENSG)
        genesdata[transcript.ENSG].TRANSCRIPTS[transcript.ENST] = transcript

    return genesdata


def write_out(source_compressed_gtf, options, transIDs, genesdata):
    # write the text file and the sorted list (in temp.txt)
    if os.path.exists('temp.txt'): os.remove('temp.txt')
    write_temp(source_compressed_gtf.replace('gtf', 'txt').replace('.gz', ''), options, transIDs, genesdata)
    sortedRecords = sort_tmpfile('temp.txt')
    writeToFile(sortedRecords, source_compressed_gtf.replace('gtf.gz', 'db').replace('gtf', 'db'))
    return sortedRecords


def parse_gtf_loop(source_compressed_gtf, options, genesdata, transIDs):
    transcript, prevenst, first, genesdata = parse_GTF(filename=source_compressed_gtf,
                                                       genesdata=genesdata)
    sys.stdout.write('Done\n')
    sys.stdout.flush()

    # Finalize last transcript and add to Gene object
    if transcript is not None:
        genesdata = finalize_last_tx(transcript, genesdata)

    # If no transcript ID from the input file was found in the  release
    if len(genesdata) == 0:
        print('\n\nNo transcripts from found in release.')
        print('\nNo transcript database created.')
        print("-----------------------------------------------------------------\n")
        quit()

    records = write_out(source_compressed_gtf, options, transIDs, genesdata)
    return len(records)


def report_summary(enst_parsed, options, hg19=False, enst=True):
    options.ensembl
    if hg19 is True:
        hg19 = '.hg19_converted'
    else:
        hg19 = ''
    if enst:
        enst = '.select_ensembl_genomic'
    else:
        enst = '.select_refseq_genomic'
    full_name = os.path.join(options.output_dir, 'MANE.GRCh38.v' + options.ensembl + enst + hg19)
    if enst_parsed > 0:
        print('\n######################################################################')
        print('A total of ' + str(enst_parsed) + ' transcripts have been retrieved\n')
        # Indexing output file with Tabix
        outfile = full_name + '.db'
        indexFile(outfile)
        """
        # Printing out summary information
        print('')
        print('---------------------')
        print('Output files created:')
        print('---------------------')
        print(outfile+ '.gz (transcript database)')
        print(outfile + '.gz.tbi (index file)')
        print(outfile.replace('.db','') + '.txt (list of transcripts)')
        """
