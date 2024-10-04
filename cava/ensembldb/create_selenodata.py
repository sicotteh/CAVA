
from urllib.request import urlopen
from time import sleep
import re
import xml.etree.ElementTree as ET
import os
import sys
from optparse import OptionParser


# Read transcript IDs from file
def readIDs(inputfile):
    ret = []
    if inputfile is not None and len(inputfile)>0:
        try:
            for line in open(inputfile):
                if not line.startswith("#"):
                    ret.append(line.strip())
        except FileNotFoundError as e:
            sys.stderr.write(e)
    return ret

# Annotation for SECIS elements only exists at NCBI, so we have to
# query NCBI by refsef id. .. which is why we need the mappings.
def load_ensembl_mappings(ensembl_refseq_file):
    refseq2ensembl = None
    refseqids = None
    try:
        f = open(ensembl_refseq_file)
        refseq2ensembl = dict()
        refseqids = []
        for line in f:
            l = line.strip().split("\t")
            if len(l)>=2:
                if l[0].startswith("ENST"):
                    refseqids.append(l[1])
                    if l[1] in refseq2ensembl:
                        refseq2ensembl[l[1]] = refseq2ensembl[l[1]] +"," + l[0]
                    else:
                        refseq2ensembl[l[1]] = l[0]

    except FileNotFoundError as e:
        sys.stderr.write(e)

    return refseq2ensembl,refseqids


def getRefseqUid(refseq_version):
    SECISelem = []
    urls=["https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term="+refseq_version+"[ACCN]&retmax=1000",
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term="+refseq_version+"[ACCN]&retmax=1000",
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term="+refseq_version+"&retmax=1000"]

    matchids = None
    iurls = 0
    url = urls[iurls]
    while(url is not None and matchids is None):
        sleep(0.151)
        sys.stderr.write("posting to get old refseq: " + url + "\n")
        page = urlopen(url)
        html_bytes = page.read()
        html = html_bytes.decode("utf-8")


        matchids = re.search(r'<IdList>\n(.*)\n</IdList>',html,re.DOTALL)
        if matchids is None:
            if iurls<len(urls)-1:
                iurls+=1
                url = urls[iurls]
    if matchids is None:
        return None
    allidsstr = matchids.group(1)
    allidstags = allidsstr.splitlines()
#   look for uids
    allids = []
    for elem in allidstags:
        matchid = re.match(r'<Id>([0-9]+)</Id>',elem)
        if matchid is not None:
            newid = matchid.group(1)
            if newid not in allids:
                allids.append(newid)
    if len(allids) == 0:
        sys.stderr.write("WARNING: no uids from "+str(len(allidstags))+" tags, for query="+url+"\n")
        return None
    return allids
def getDataFromNCBI(allids):
    SECIS = dict()
    for uid in allids:
        sleep(0.151)
        iurl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id="+uid+"&rettype=gb&retmode=xml"
        sys.stdout.write('posting ' + iurl + '\n')
        ipage = urlopen(iurl)
        ihtml_bytes = ipage.read()
        ihtml = ihtml_bytes.decode("utf-8")
        root = ET.fromstring(ihtml)
        nmseq = ""
        nmlen = ""
        nmelem = root.find('GBSeq/GBSeq_locus')
        if nmelem is not None:
            nmseq = nmelem.text
            if not nmseq.startswith("NM_"):
                nmseq = ""
        if len(nmseq) == 0 and nmelem is not None:
            nmelem = root.find('GBSeq/GBSeq_primary-accession')
            if nmelem is not None:
                nmseq = nmelem.text

        nmlen_elem = root.find('GBSeq/GBSeq_length')
        if nmlen_elem is not None:
            nmlen = nmlen_elem.text

        nmelemv = root.find('GBSeq/GBSeq_accession-version')
        nmseqv = ""
        if nmelemv is not None:
            nmseqv = nmelemv.text
        else:
            nmelemv = root.find('GBSeq/GBSeq_other-seqids')
            if nmelemv is not None:
                for gbseqid in nmelemv:
                    idsv = gbseqid.text.split('|')
                    if len(idsv) ==2 and idsv[0] == 'ref' and idsv[1].startswith('NM_'):
                        if re.match('^NM_[0-9]+\.[0-9]+$',idsv[1]):
                            nmseqv = idsv[1]


        feature_table = root.find('GBSeq/GBSeq_feature-table')
        foundSECIS = False
        gene = ""
        startpos = ""
        endpos = ""
        accnvpos = ""
        secpos = ""
        secpos2 = ""
        cds_start = ""
        cds_end = ""
        for feat in feature_table.findall('GBFeature'):
            fkey = feat.find('GBFeature_key')
            if fkey is not None:
                if fkey.text == 'regulatory':
                    gbint = feat.find('GBFeature_intervals/GBInterval')
                    if gbint is not None:
                        gbquals = feat.find('GBFeature_quals')
                        found_this_SECIS = False
                        for gbqual in gbquals:
                            gbqual_name = gbqual.find('GBQualifier_name').text
                            gbqual_value = gbqual.find('GBQualifier_value').text
                            if gbqual_name == 'gene':
                                gene = gbqual_value
                            elif gbqual_name == 'note' and gbqual_value == 'SECIS_element':
                                found_this_SECIS = True
                        if found_this_SECIS is True:
                            foundSECIS = True
                            newstartpos = gbint.find('GBInterval_from').text
                            newendpos = gbint.find('GBInterval_to').text
                            newaccnvpos = gbint.find('GBInterval_accession').text
                            if len(startpos) > 0:
                                if len(newstartpos) > 0:
                                    if int(newstartpos) > int(startpos):  # further toward the end of the UTR has biggest reach/
                                        startpos = newstartpos
                                        endpos = newendpos
                                        accnvpos = newaccnvpos
                            else:
                                startpos = newstartpos
                                endpos = newendpos
                                accnvpos = newaccnvpos
                elif fkey.text == "CDS":
                    cds_loc_elem = feat.find('GBFeature_location')
                    if cds_loc_elem is not None:
                        cds_loc = cds_loc_elem.text
                        mloc = re.match(r"([0-9]+)\.\.([0-9]+)", cds_loc)
                        if mloc is not None:
                            cds_start = mloc.group(1)
                            cds_end = mloc.group(2)

                    gbquals = feat.find('GBFeature_quals')
                    for gbqual in gbquals:
                        gbqual_name_elem = gbqual.find('GBQualifier_name')
                        if gbqual_name_elem is not None:
                            gbqual_name = gbqual_name_elem.text
                            if gbqual_name == 'transl_except':
                                gbqual_value = gbqual.find('GBQualifier_value').text
                                mtrans = re.match(r"\(pos:([0-9]+)\.\.([0-9]+),aa:Sec",gbqual_value)
                                if mtrans is not None:
                                    newsecpos = mtrans.group(1)
                                    newsecpos2 = mtrans.group(2)
                                    if len(secpos2)>0:
                                        if int(newsecpos2)> int(secpos2): #latest UGA recoded guarantees SECIS is working that far
                                            secpos = newsecpos
                                            secpos2 = newsecpos2
                                    else:
                                        secpos = newsecpos
                                        secpos2 = newsecpos2

            if foundSECIS is True:
                if accnvpos is not None:
                    if accnvpos not in SECIS:
                        SECIS[accnvpos] = [gene,startpos,endpos,accnvpos,secpos,cds_start,cds_end,nmlen]
                elif nmseqv is not None:
                    if nmseqv not in SECIS:
                        SECIS[nmseqv] = [gene,startpos,endpos,nmseqv,secpos,cds_start,cds_end,nmlen]
                elif nmseq is not None:
                    if nmseq not in SECIS:
                        SECIS[nmseq] = [gene,startpos,endpos,nmseq,secpos,cds_start,cds_end,nmlen]

# UGA codons ending 52 bp or greater BEFORE SECIS element can recode as Selenocysterin.
    return SECIS

def getSel(outfile,organism,selenogenes,refseqids,refseq2ensembl):

    # scan for refseq transcripts with SECIS elements
    org = organism.replace('_','%20')
    url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=SECIS%5bAll%20Fields%5d+AND+srcdb_refseq%5bPROP%5d+AND+%22'+org+'%22%5bORGANISM%5d+AND+%22biomol%20mrna%22%5bProperties%5d&retmax=1000'
    page = urlopen(url)
    html_bytes = page.read()
    html = html_bytes.decode("utf-8")


    matchids = re.search(r'<IdList>\n(.*)\n</IdList>',html,re.DOTALL)
    allidsstr = matchids.group(1)
    allidstags = allidsstr.splitlines()
#   look for uids
    allids = []
    for elem in allidstags:
        matchid = re.match(r'<Id>([0-9]+)</Id>',elem)
        if matchid is not None:
            newid = matchid.group(1)
            if newid not in allids:
                allids.append(newid)

    # scan for refseq transcripts in known selenogenes (ideally, this should not return any additional transcripts)
    for gene in  selenogenes:
        sleep(0.201)  # max rate is 10 request/min . so sleep at least 0.1 secs
        gurl = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=SECIS%5bAll%20Fields%5d+AND+"+gene+"%5bGENE%5d+AND+srcdb_refseq%5bPROP%5d+AND+%22'+org+'%22%5bORGANISM%5d+AND+%22biomol%20mrna%22%5bProperties%5d&retmax=1000'
        gpage = urlopen(gurl)
        sys.stdout.write('posting '+gurl+'\n')
        ghtml_bytes = gpage.read()
        ghtml = ghtml_bytes.decode("utf-8")
        matchids = re.search(r'<IdList>\n(.*)\n</IdList>', ghtml, re.DOTALL)
        if matchids  is not None:
            allidsstr = matchids.group(1)
            allidstags = allidsstr.splitlines()
            for elem in allidstags:
                matchid = re.match(r'<Id>([0-9]+)</Id>', elem)
                if matchid is not None:
                    newid = matchid.group(1)
                    if newid not in allids:
                        allids.append(newid)
                else:
                    sys.stdout.write("No ids in : "+elem)
# get all SECIS info for transcripts (NM.most_curent_version) from NCBI
    SECIS = getDataFromNCBI(allids)
    # If user requested specific refseqids (versions)
    #  .. then the current SECIS list identifies which refsew id have SECIS elements .
    # ... so can fetch a different "version" of the transcript if needed for
    #  .. different builds.
    if refseqids is not None and len(refseqids) > 0:
        refseqnv_to_refseqv = dict()
        for ridv in refseqids:
            refseqnv_to_refseqv[ridv.split(".")[0]] = ridv
        accnv_is_SECIS = []
        accn_is_SECIS = []
        new_SECIS = dict()
        for accnv in SECIS:
            if re.match('^NM_[0-9]+\.[0-9]+$',accnv):
                accnv_is_SECIS.append(accnv)
                accn = accnv.split('.')[0]
            else: # no period in accnv
                accn_is_SECIS.append(accnv)
                accn = accnv
            if accn in refseqnv_to_refseqv:
                wanted_refseqv = refseqnv_to_refseqv[accn]
                if accnv == wanted_refseqv:
                    new_SECIS[accnv] = SECIS[accnv]
                else:
                    sleep(0.15)
                    refsequids = getRefseqUid(wanted_refseqv)
                    nrs = 0
                    if refsequids is not None and len(refsequids)>0:
                        SECISelems = getDataFromNCBI(refsequids)
                        if SECISelems is not None:
                            for uid in SECISelems:
                                SECISelem = SECISelems[uid]
                                if SECISelem is not None:
                                    new_SECIS[wanted_refseqv] = SECISelem
                                    nrs+=1
                    if nrs == 0: # This is risky .. assume that old NM transcript is the same..
                        sys.stderr.write("WARNING: No SECIS annotation found for "+wanted_refseqv+
                                        " for SECIS element, using annotation of "+accnv+ " instead\n")
                        new_SECIS[wanted_refseqv] = SECIS[accnv]
        SECIS = new_SECIS

    fout = open(outfile,mode='w')
    fout.write("gene\tSECIS_maxstart\tSECIS_maxend\taccn\tlast_sec_pos\tcds_start\tcds_end\tmRNA_len\n")
    for accnv in SECIS:
        elemv = SECIS[accnv]
        if refseq2ensembl is not None:
            if accnv in refseq2ensembl:
                ensemblids = refseq2ensembl[accnv].split(",")
                for eids in ensemblids:
                    elemv[3] = eids
                    fout.write("\t".join(elemv) + "\n")
        else:
            fout.write("\t".join(elemv)+"\n")
    fout.close()

# search for /regulatory POS..POS
# # /regulatory_class="recoding_stimulatory_region"
# note SECIS element
#/note="SECIS_element"fse



def main():
# Command line argument parsing
    usage_str = "\nExample usage: Run after catalog file were created: Without any option, will build for latest Refseqs for Human.\n python3 create_selenodata.py \n"
    parser = OptionParser()
    parser.usage = usage_str

    parser.add_option('-e', "--ensembl_refseq", dest='ensembl_refseq', default = None, action='store', help="tab-delimited file of ensembl and matching refseqs.versiom. If specified, overrides the --refseqs option")
    parser.add_option('-r', "--refseqs", dest='refseqs', action='store', help="file with refseq.version, the version number can make the refseqs build-=specific.")
    parser.add_option('-D', "--outdir", dest='output_dir', action='store', default='data', help="Output directory")
    parser.add_option('-g', "--genes", dest='genes', action='store', default='', help="file containes list of known selenocysteine genes (organism-specific)")
    parser.add_option('-o', "--organism", dest='organism', action='store', default='homo_sapiens', help="Organism to use for the entrez query (defaults to homo_sapiens). \nNo spaces in name (e.g. mus_musculus,sus_scrofa)")
    parser.add_option('-t', "--tag", dest='build_tag', action='store', help="arbitrary tag to put in output filenames. Defaults  to refseq or ensembl.")# organism should be the same syntax as on the ftp site for ensembl.
# https://ftp.ensembl.org/pub/release-112/regulation/ e.g. homo_sapient mus_musculus, etc...
    (options, args) = parser.parse_args()
    options.select = False

    if not os.path.exists(options.output_dir):
        os.mkdir(options.output_dir)
#    getSel(os.path.join(options.output_dir,'SECIS_in_refseq_pos.txt').)
    id_type = "refseq"
    if options.refseqs and options.ensembl_refseq:
        sys.stderr.write("ERROR: create_selenodata specified both refseqs and ensembl_refseq. Select only 1")
        exit(1)
    elif options.ensembl_refseq:
        id_type = "ensembl"


    if options.genes:
        selenogenes = readIDs(options.genes)
    else:
        selenogenes = [     'DIO1', '1733', 'TXDI1', 'THMA2',
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
                            'SELENOU', 'not human',
                            'SELENOV', '348303', 'SELV',
                            'SELENOW',  '6415', 'SEPW1', 'SELW',
                            'SELENOP1', 'not in human', 'SELENOPZ', 'SEPP1',
                            'SELENOP2', 'not in human', 'SELPB', 'SEPP1L', 'SEPP2',
                            'MSRB1',  '51734', 'SELENOX', 'SELENOR', 'SELR', 'SELX', 'SEPX1', 'SEPR', 'HSOC270',
                            'SEPHS2', '22928',  'SPS2', 'SPS2B',
                            'TXNRD1',  '7296', 'TXNR', 'GRIM-12', 'TRXR1', 'TXNR1', 'TR', 'TR1', 'TRXR1',
                            'TXNRD2', '10587', 'SELZ', 'TR',  'TRXR2',  'TR3', 'TXNR2', 'GCCD5',  'TR',  'TR-BETA',
                            'TXNRD3', '114112',  'TXNRD3NB', 'TXNRD3IT1', 'TR2', 'TRXR3', 'TGR', 'TR2IT1', 'TXNRD3NT1', 'TXNR3', 'TGR']


    if options.ensembl_refseq:
        refseq2ensembl,refseqids = load_ensembl_mappings(options.ensembl_refseq)
    elif options.refseqs:
        refseqids = readIDs(options.refseqs)
        refseq2ensembl = None
    else: #Load most recent refseqs
        refseqids = None
        refseq2ensembl = None
    if options.build_tag is not None and len(options.build_tag)>0:
        outfile = 'SECIS_in_refseq_pos.' + options.organism + '.' + options.build_tag + '.txt'
    else:
        outfile = 'SECIS_in_refseq_pos.' + options.organism + '.' + id_type + '.txt'

    getSel(os.path.join(options.output_dir,outfile),options.organism,
           selenogenes,refseqids,refseq2ensembl)

if __name__ == '__main__':
    main()


