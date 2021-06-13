import datetime
import gzip
import os
import pickle
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

    outfile_list.write(
        '# Created by CAVA ensembl_db ' + options.version + ' based on Ensembl release ' + options.ensembl + '\n')
    outfile_list.write('#GENEID\tSYMBOL\tTranscript\tProtein\n')

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
    print(f'Parsing {filename}', end="...")
    tx_to_prot_dict = build_tx_to_prot_dict(opener, filename)

    for line in opener(filename, 'rt'):

        line = line.strip()
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

        # Do not consider transcript if it is not on the custom transcript list
        if options.input is not None and enst not in transIDs: continue

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
            transcript.GENE = getValue(tags, 'gene_name')
            transcript.ENSG = getValue(tags, 'gene_id')
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
            transcript.isComplete = True
            transcript.CCDS = True
        prevenst = enst
        if first: first = False
    return transcript, prevenst, first, genesdata


def sort_tmpfile(f):
    # Sort temporary output file
    data = dict()
    counter = 0
    for line in open(f, 'r'):
        if not line.startswith('ENST'): continue
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
              '20', '21', '22', '23', 'M', 'MT', 'X', 'Y']
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
def process_data(options, genome_build):
    # Dictionary of Gene objects
    genesdata = dict()

    # Load custom transcript IDs
    transIDs = None
    if options.input is not None:
        transIDs = readTranscriptIDs(options.input)
        print('\nOnly ' + str(len(transIDs)) + ' transcripts read from ' + options.input + ' are considered\n')
    else:
        print('\nAll transcripts from the Ensembl release are considered\n')

    # Load candidate and CCDS data for Ensembl <75
    candidates = dict()
    if int(options.ensembl) < 75:
        datadir = os.path.dirname(os.path.realpath(__file__)) + '/data'
        for line in open(datadir + '/info' + options.ensembl + '.txt'):
            line = line.strip()
            if line == '': continue
            cols = line.split('\t')
            if cols[0] not in list(candidates.keys()): candidates[cols[0]] = dict()
            candidates[cols[0]][cols[1]] = int(cols[2])

    ######################################################################

    # Download Ensembl data if necessary
    source_compressed_gtf = 'Homo_sapiens.' + genome_build + '.' + options.ensembl + '.gtf.gz'
    source_compressed_gtf = os.path.join('data', source_compressed_gtf)
    if not os.path.exists(source_compressed_gtf):
        sys.stdout.write('Downloading Ensembl database... ')
        sys.stdout.flush()

        url = 'ftp://ftp.ensembl.org/pub/release-' + options.ensembl + '/gtf/homo_sapiens/Homo_sapiens.' + genome_build + '.' + options.ensembl + '.gtf.gz'
        try:
            wget.download(url)
            os.rename('Homo_sapiens.' + genome_build + '.' + options.ensembl + '.gtf.gz', source_compressed_gtf)
        except Exception as e:
            print('\n\nCannot connect to Ensembl FTP site. No internet connection?\n')
            print(f'{e}\n{url}')
            quit()

    ################################################################
    # Use crossmap to get hg19 if desired
    #################################################################
    if options.no_hg19 is not False:
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
        converted_gtf = os.path.join('data', 'Homo_sapiens.hg19_converted' + options.ensembl + '.gtf')
        crossmap_gff_file(mapTree, source_compressed_gtf, converted_gtf)

        # Note this file is not sorted!
        a = pybedtools.BedTool(converted_gtf)
        a.sort().remove_invalid().saveas('tmp.txt')
        os.rename('tmp.txt', converted_gtf)

    ################################################################
    #
    #################################################################
    # Iterate through the lines in the ensembl data file
    sys.stdout.write('Extracting transcript data from Ensembl...')

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
        print('\n\nNo transcripts from ' + options.input + ' found in Ensembl release.')
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
    if options.no_hg19 is not False:
        sys.stdout.write('Extracting transcript data for hg19 version...')
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

        # If no transcript ID from the input file was found in the Ensembl release
        if len(genesdata) == 0:
            print('\n\nNo transcripts from ' + options.input + ' found in Ensembl release.')
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
    # os.remove(source_compressed_gtf)

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
    if options.ensembl is None:
        print('\nError: no Ensembl release specified. Use option -h to get help!\n')
        quit()
    if not is_number(options.ensembl):
        print('\nError: Ensembl release specified is not an integer. Use option -h to get help!\n')
        quit()
    if options.output is None:
        print('\nError: no output file name specified. Use option -h to get help!\n')
        quit()

    # Must use Ensembl release >= 70
    if not (int(options.ensembl) >= 70 or int(options.ensembl) == 65):
        print('\nError: This version.py works with Ensembl v65 or >= v70.\n')
        quit()

    # Genome build
    # genome_build = options.genome
    genome_build = 'GRCh37' if int(options.ensembl) <= 75 else 'GRCh38'

    # Printing out version.py information
    print("\n---------------------------------------------------------------------------------------")
    print('CAVA ' + options.version + ' transcript database preparation tool (ensembl_db) is now running')
    print('Started: ', datetime.datetime.now(), '\n')

    # Print info
    print('Ensembl version.py:  ' + options.ensembl)
    print('Reference genome: ' + genome_build)

    # Creating compressed output file
    enst_parsed, ens_lifted = process_data(options, genome_build)
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
    print('CAVA ensembl_db successfully finished: ', datetime.datetime.now())
    print("---------------------------------------------------------------------------------------\n")
