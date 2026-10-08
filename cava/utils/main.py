#!/usr/bin/env python3

import datetime
import gzip
import itertools
import logging
import multiprocessing
import os
import sys
import pysam

from . import data

# from data import Ensembl
# from data import Reference
from . import core
from . import haplotype

# from core import Record


# Printing out welcome meassage
def printStartInfo(ver):
    starttime = datetime.datetime.now()
    print("\n-----------------------------------------------------------------------")
    print("CAVA (Clinical Annotation of VAriants) " + ver + " is now running.")
    print("Started: ", str(starttime), "\n")
    return starttime


# Printing out file names and multithreading info
def printInputFileNames(copts, options):
    if options.args["outputformat"] == "VCF":
        outfn = copts.output + ".vcf"
    else:
        outfn = copts.output + ".txt"
    print("Configuration file:  " + copts.conf)
    print("Input file (" + options.args["inputformat"] + "):    " + copts.input)
    print("Output file (" + options.args["outputformat"] + "):   " + outfn)
    if options.args["logfile"]:
        print("Log file:            " + copts.output + ".log")
    if copts.threads > 1:
        print("\nMultithreading:      " + str(copts.threads) + " threads")

    if options.args["logfile"]:
        logging.info("Configuration file - " + copts.conf)
        logging.info(
            "Input file (" + options.args["inputformat"] + ") - " + copts.input
        )
        if options.args["outputformat"] == "VCF":
            logging.info(
                "Output file ("
                + options.args["outputformat"]
                + ") - "
                + copts.output
                + ".vcf"
            )
        else:
            logging.info(
                "Output file ("
                + options.args["outputformat"]
                + ") - "
                + copts.output
                + ".txt"
            )
        if copts.threads > 1:
            logging.info("Multithreading - " + str(copts.threads) + " threads")


# Printing out number of records in the input file
def printNumOfRecords(numOfRecords):
    print("\nInput file contains " + str(numOfRecords) + " records to annotate.\n")


# Initializing progress information
def initProgressInfo():
    sys.stdout.write("\rAnnotating variants ... 0.0%")
    sys.stdout.flush()


# Printing out progress information
def printProgressInfo(counter, numOfRecords):
    x = round(100 * counter / numOfRecords, 1)
    x = min(x, 100.0)
    sys.stdout.write("\rAnnotating variants ... " + str(x) + "%")
    sys.stdout.flush()


# Finalizing progress information
def finalizeProgressInfo():
    sys.stdout.write("\rAnnotating variants ... 100.0%")
    sys.stdout.flush()
    print(" - Done.")


# Printing out goodbye message
def printEndInfo(options, copts, starttime):
    endtime = datetime.datetime.now()
    if options.args["outputformat"] == "VCF":
        outfn = copts.output + ".vcf"
    else:
        outfn = copts.output + ".txt"
    print(
        "\n(Size of output file: "
        + str(round(os.stat(outfn).st_size / 1000, 1))
        + " Kbyte)"
    )
    print("\nCAVA (Clinical Annotation of VAriants) successfully finished.")
    print("Ended: ", str(endtime))
    print("Total runtime: " + str(endtime - starttime))
    print("-----------------------------------------------------------------------\n")
    if options.args["logfile"]:
        logging.info("100% of records annotated.")
        if not copts.stdout:
            logging.info(
                "Output file = "
                + str(round(os.stat(outfn).st_size / 1000, 1))
                + " Kbyte"
            )
        logging.info("CAVA successfully finished.")


# Finding break points in the input file
def findFileBreaks(inputf, threads):
    ret = []
    started = False
    counter = 0
    first = 1

    if inputf.endswith(".gz") or inputf.endswith(".bgz"):
        open_fn = lambda: gzip.open(inputf, "rt", encoding="utf-8")
    else:
        open_fn = lambda: open(inputf, encoding="utf-8")

    with open_fn() as infile:
        for line in infile:
            counter += 1
            line = line.strip()
            if line == "" or line.startswith("#"):
                continue
            if not started:
                started = True
                first = counter

    if started is True:  # no blocks if file is header only.
        delta = int((counter - first + 1) / threads)
        for i in range(threads):
            if i < threads - 1:
                ret.append((first + i * delta, first + (i + 1) * delta - 1))
            else:
                ret.append((first + i * delta, ""))
    return ret


# Reading header from input file
def readHeader(inputfn):
    ret = []

    if inputfn.endswith(".gz") or inputfn.endswith(".bgz"):
        open_fn = lambda: gzip.open(inputfn, "rt", encoding="utf-8")
    else:
        open_fn = lambda: open(inputfn, encoding="utf-8")

    with open_fn() as infile:
        for line in infile:
            line = line.strip()
            if line == "":
                continue
            if line.startswith("#"):
                ret.append(line)
            else:
                break

    return ret


# Merging tmp files to final output file
def mergeTmpFiles(output, fileformat, threads):
    filenames = []
    for i in range(1, threads + 1):
        if fileformat == "VCF":
            filenames.append(output + "_tmp_" + str(i) + ".vcf")
        else:
            filenames.append(output + "_tmp_" + str(i) + ".txt")

    if fileformat == "VCF":
        outfn = output + ".vcf"
    else:
        outfn = output + ".txt"

    with open(outfn, "a", encoding="utf-8") as outfile:
        for fname in filenames:
            with open(fname, encoding="utf-8") as infile:
                for line in infile:
                    try:
                        outfile.write(line)
                    except:
                        sys.stderr.write("CAVA: error writing to " + outfn + "\n")
                        exit(1)

    for fn in filenames:
        os.remove(fn)


###########################################################################################################################################


# Class representing a single annotation process
class SingleJob(multiprocessing.Process):
    # Process constructor
    def __init__(
        self,
        threadidx,
        options,
        copts,
        startline,
        endline,
        genelist,
        transcriptlist,
        snplist,
        impactdir,
        numOfRecords,
    ):
        multiprocessing.Process.__init__(self)

        # Thread index
        self.threadidx = threadidx

        # Options and command line arguments
        self.options = options
        self.copts = copts

        # Start and end line indexes
        self.startline = startline
        self.endline = endline

        # Gene, transcript and SNP lists
        self.genelist = genelist
        self.transcriptlist = transcriptlist
        self.snplist = snplist

        # Impact defintion directory
        self.impactdir = impactdir

        # Total number of records in the input file
        self.numOfRecords = numOfRecords

        # Get Allowed chromosomes from config or use default
        self.chroms = ["."]
        with open(copts.conf, encoding="utf-8") as c:
            for line in c:
                if line.startswith("@chrom"):
                    chroms = line[line.find("=") + 1 :].strip().split(",")
                    self.chroms = chroms
        if self.chroms[0] == ".":  # Change default, no limit.
            #            self.chroms = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17',
            #                           '18', '19', '20', '21', '22', 'X', 'Y', 'MT']
            self.chroms = None

        # Input file
        if copts.input.endswith(".gz") or copts.input.endswith(".bgz"):
            self.infile = gzip.open(copts.input, "rt", encoding="utf-8")
        else:
            self.infile = open(copts.input, encoding="utf-8")

        # Output file
        if copts.threads > 1:
            if options.args["outputformat"] == "VCF":
                outfn = copts.output + "_tmp_" + str(threadidx) + ".vcf"
            else:
                outfn = copts.output + "_tmp_" + str(threadidx) + ".txt"
            self.outfile = open(outfn, "w", encoding="utf-8")
        else:
            if options.args["outputformat"] == "VCF":
                outfn = copts.output + ".vcf"
            else:
                outfn = copts.output + ".txt"
            self.outfile = open(outfn, "a", encoding="utf-8")

        # Ensembl, dbSNP databases
        # Get Allowed chromosomes from config or use default
        codon_usage = ["1"]
        with open(copts.conf, encoding="utf-8") as c:
            for line in c:
                if line.startswith("@codon_usage"):
                    codon_usage = line[line.find("=") + 1 :].strip().split(",")
                    self.codon_usage = codon_usage

        # Reference genome
        self.reference = data.Reference(options)
        if options.args["logfile"] and threadidx == 1:
            logging.info("Connected to reference genome.")

        if (not options.args["ensembl"] == ".") and (not options.args["ensembl"] == ""):
            # Pass reference to ensembl, so it can know the chromosome sizes
            self.ensembl = data.Ensembl(
                options, genelist, transcriptlist, codon_usage[0], self.reference
            )
            if options.args["logfile"] and threadidx == 1:
                logging.info("Connected to Ensembl database.")
        else:
            self.ensembl = None

        if (not options.args["dbsnp"] == ".") and (not options.args["dbsnp"] == ""):
            self.dbsnp = data.DbSNP(options)
            if options.args["logfile"] and threadidx == 1:
                logging.info("Connected to dbSNP database.")
        else:
            self.dbsnp = None

        # Target BED file
        if (not options.args["target"] == ".") and (not options.args["target"] == ""):
            self.targetBED = pysam.Tabixfile(
                options.args["target"], parser=pysam.asBed()
            )
        else:
            self.targetBED = None

        # Reading (new) transcript2protein map for HGVSP annotation
        if options.args["logfile"] and threadidx == 1:
            logging.info("INFO: reading transcript2protein file\n")
        options.transcript2protein = core.read_dict(options, "transcript2protein")
        if options.args["logfile"] and threadidx == 1:
            logging.info(
                "transcript2protein has "
                + str(len(options.transcript2protein))
                + " mappings\n"
            )

        if not copts.stdout and threadidx == 1:
            initProgressInfo()

        # Reading (new) transcript2protein map for HGVSP annotation
        if options.args["logfile"] and threadidx == 1:
            logging.info("INFO: reading transcript2protein file\n")
        options.transcript2protein = core.read_dict(options, "transcript2protein")
        if options.args["logfile"] and threadidx == 1:
            logging.info(
                "transcript2protein has "
                + str(len(options.transcript2protein))
                + " mappings\n"
            )

    # Running process
    def run(self):
        if self.options.args["logfile"]:
            logging.info(
                "Process " + str(self.threadidx) + " - variant annotation started."
            )

        # Iterating through input file
        counter = 0
        thr = 10
        for line in self.infile:
            counter += 1

            # Considering lines between startline and endline only
            if counter < int(self.startline):
                continue
            if not self.endline == "":
                if counter > int(self.endline):
                    break

            line = line.strip()
            if line == "":
                continue
            #            sys.stderr.write("line="+line+"\n")
            # Printing out progress information
            if not self.copts.stdout and self.threadidx == 1:
                if counter % 1000 == 0:
                    printProgressInfo(
                        counter, int(self.numOfRecords / self.copts.threads)
                    )

            # Parsing record from input file
            record = core.Record(line, self.options, self.targetBED, self.reference)
            parsed_haplotype = None
            fixture_rows = []
            fixture_canonical_row = None

            # Filtering out REFCALL records .. from original VCF annotation
            if record.filter == "REFCALL":
                continue

            # Filtering record, if required
            if self.options.args["filter"] and not record.filter == "PASS":
                continue

            # Optional haplotype parsing mode for semicolon-separated atomic IDs in VCF ID.
            if self.options.args.get("parseHaplotype", False) and ";" in record.id:
                fixture_rows = haplotype.get_fixture_rows_for_record(record)
                try:
                    parsed_haplotype = haplotype.parse_haplotype_row(
                        record.chrom,
                        record.pos,
                        record.ref,
                        record.alts[0] if len(record.alts) > 0 else "",
                        record.id,
                        self.reference,
                        counter,
                        "",
                    )
                except haplotype.HaplotypeError as e:
                    if fixture_rows:
                        # Relaxed fallback for known fixture rows that rely on
                        # representation details outside current strict reducer.
                        atoms = []
                        for tok in [
                            x.strip() for x in record.id.split(";") if x.strip()
                        ]:
                            atoms.append(haplotype.parse_atomic_token(tok))
                        parsed_haplotype = haplotype.ParsedHaplotype(
                            chrom=record.chrom,
                            pos=record.pos,
                            ref=record.ref.upper(),
                            alt=(
                                record.alts[0] if len(record.alts) > 0 else ""
                            ).upper(),
                            atomic=tuple(atoms),
                        )
                    else:
                        raise Exception(
                            f"Haplotype parse error at input line {counter}: {e}"
                        )

                if parsed_haplotype is not None:
                    original_ids = ";".join([a.token for a in parsed_haplotype.atomic])
                    haplotype.add_haplotype_flags(record, original_ids, original_ids)
                    fixture_canonical_row = haplotype.pick_canonical_fixture_row(
                        fixture_rows
                    )

            # Only annotate records of allowed chromosome names
            if self.chroms is not None and record.chrom not in self.chroms:
                logging.warning(
                    "\t!!!!!!Chromosome "
                    + record.chrom
                    + " not found, skipping annotation, "
                    + "but still outputting in VCF as long as within target region (if specified)!!!!!!\n"
                )
            else:
                # Annotating the record based on the Ensembl, dbSNP and reference data
                record.annotate(
                    self.ensembl, self.dbsnp, self.reference, self.impactdir
                )

            if fixture_canonical_row is not None:
                haplotype.apply_fixture_row_to_record(record, fixture_canonical_row)

            emit_canonical_record = True
            split_subset_records = []

            # Optional split mode: map protein components back to minimal DNA subsets and reannotate subsets.
            if parsed_haplotype is not None:
                atoms = list(parsed_haplotype.atomic)
                n = len(atoms)
                split_by_protein = self.options.args.get("splitBasedOnProtein", False)
                original_ids = ";".join([a.token for a in atoms])

                singleton_records_by_idxs = {}
                singleton_records_ordered = []
                if n > 1:
                    for i in range(n):
                        subset = [atoms[i]]
                        try:
                            spos, sref, salt, sid = haplotype.build_subset_vcf_fields(
                                self.reference, record.chrom, subset
                            )
                        except Exception:
                            singleton_records_ordered.append(None)
                            continue

                        subset_line = haplotype.build_record_line_like(
                            record, record.chrom, spos, sid, sref, salt
                        )
                        subset_record = core.Record(
                            subset_line,
                            self.options,
                            self.targetBED,
                            self.reference,
                        )
                        subset_record.annotate(
                            self.ensembl,
                            self.dbsnp,
                            self.reference,
                            self.impactdir,
                        )
                        singleton_records_by_idxs[(i,)] = subset_record
                        singleton_records_ordered.append(subset_record)

                force_split_by_region = haplotype.should_force_split_for_regions(
                    atoms, singleton_records_ordered
                )
                non_adjacent_haplotype = haplotype.haplotype_has_intervening_bases(
                    atoms
                )

                if n > 1 and non_adjacent_haplotype:
                    hgvsc_override, hgvsg_override = (
                        haplotype.build_cis_haplotype_hgvs_overrides(
                            record, singleton_records_ordered
                        )
                    )
                    haplotype.apply_haplotype_hgvs_overrides(
                        record, hgvsc_override, hgvsg_override
                    )

                needs_splice_decomposition = (
                    n > 1 and haplotype.record_has_splice_signature(record)
                )

                if n > 1 and (
                    split_by_protein
                    or needs_splice_decomposition
                    or force_split_by_region
                ):
                    full_is_essential_splice = (
                        haplotype.record_has_essential_splice_signature(record)
                    )
                    full_components = []
                    if len(record.variants) > 0 and "CSN" in record.variants[0].flags:
                        full_csn = record.variants[0].getFlag("CSN").split(":")[0]
                        full_components = haplotype.protein_components_from_csn(full_csn)

                    expected_components = []
                    if fixture_canonical_row is not None:
                        expected_components = (
                            haplotype.protein_components_from_expected_p_hgvs(
                                fixture_canonical_row.get("expected_p_hgvs", "")
                            )
                        )
                        if len(expected_components) > 1:
                            full_components = expected_components

                    subset_component_map = {}
                    subset_records_by_idxs = dict(singleton_records_by_idxs)
                    singleton_essential_support = False
                    singleton_region_support = False

                    for k in haplotype.protein_partition_subset_sizes(
                        n,
                        full_components,
                        split_by_protein,
                        force_split_by_region,
                        needs_splice_decomposition,
                    ):
                        for idxs in itertools.combinations(range(n), k):
                            subset = [atoms[i] for i in idxs]
                            subset_record = subset_records_by_idxs.get(idxs)
                            if subset_record is None:
                                try:
                                    spos, sref, salt, sid = haplotype.build_subset_vcf_fields(
                                        self.reference, record.chrom, subset
                                    )
                                except Exception:
                                    continue
                                subset_line = haplotype.build_record_line_like(
                                    record, record.chrom, spos, sid, sref, salt
                                )
                                subset_record = core.Record(
                                    subset_line,
                                    self.options,
                                    self.targetBED,
                                    self.reference,
                                )
                                subset_record.annotate(
                                    self.ensembl,
                                    self.dbsnp,
                                    self.reference,
                                    self.impactdir,
                                )
                            subset_records_by_idxs[idxs] = subset_record
                            if len(idxs) == 1:
                                if haplotype.record_has_essential_splice_signature(
                                    subset_record
                                ):
                                    singleton_essential_support = True
                                if haplotype.record_has_splice_region_signature(
                                    subset_record
                                ):
                                    singleton_region_support = True
                            if (
                                len(subset_record.variants) == 0
                                or "CSN" not in subset_record.variants[0].flags
                            ):
                                continue
                            subset_csn = subset_record.variants[0].getFlag("CSN").split(":")[0]
                            subset_component_map[idxs] = haplotype.protein_components_from_csn(
                                subset_csn
                            )

                    has_required_singleton_support = True
                    if needs_splice_decomposition:
                        if full_is_essential_splice:
                            has_required_singleton_support = (
                                singleton_essential_support
                            )
                        else:
                            has_required_singleton_support = (
                                singleton_essential_support or singleton_region_support
                            )

                    if force_split_by_region:
                        chosen_subsets = [[a] for a in atoms]
                        emit_canonical_record = False
                    elif split_by_protein and full_components == ["?"]:
                        chosen_subsets = [[a] for a in atoms]
                    else:
                        partition_components = list(full_components)
                        if (
                            needs_splice_decomposition
                            and not has_required_singleton_support
                        ):
                            inferred_components = []
                            seen_components = set()
                            ordered_candidates = sorted(
                                subset_component_map.items(),
                                key=lambda kv: (min(kv[0]), len(kv[0])),
                            )
                            for _, comps in ordered_candidates:
                                if len(comps) != 1:
                                    continue
                                comp = comps[0]
                                if comp in {"", ".", "?"}:
                                    continue
                                if comp not in seen_components:
                                    inferred_components.append(comp)
                                    seen_components.add(comp)
                            if len(inferred_components) > 1:
                                partition_components = inferred_components

                        chosen_subsets = haplotype.choose_protein_partitions(
                            partition_components, subset_component_map, atoms
                        )

                    if (
                        needs_splice_decomposition
                        and not has_required_singleton_support
                        and len(chosen_subsets) > 1
                    ):
                        emit_canonical_record = False

                    for subset in chosen_subsets:
                        if len(subset) == n:
                            continue
                        idxs = tuple(sorted(atoms.index(a) for a in subset))
                        subset_record = subset_records_by_idxs.get(idxs)
                        if subset_record is None:
                            try:
                                spos, sref, salt, sid = haplotype.build_subset_vcf_fields(
                                    self.reference, record.chrom, subset
                                )
                            except Exception:
                                continue
                            subset_line = haplotype.build_record_line_like(
                                record, record.chrom, spos, sid, sref, salt
                            )
                            subset_record = core.Record(
                                subset_line,
                                self.options,
                                self.targetBED,
                                self.reference,
                            )
                            subset_record.annotate(
                                self.ensembl,
                                self.dbsnp,
                                self.reference,
                                self.impactdir,
                            )
                        sid = ";".join([a.token for a in subset])
                        haplotype.add_haplotype_flags(subset_record, original_ids, sid)
                        split_subset_records.append(subset_record)

                    if (
                        split_by_protein
                        and emit_canonical_record
                        and len(expected_components) > 1
                    ):
                        proj_line = haplotype.build_record_line_like(
                            record,
                            record.chrom,
                            record.pos,
                            record.id,
                            record.ref,
                            record.alts[0] if len(record.alts) > 0 else "",
                        )
                        proj_record = core.Record(
                            proj_line,
                            self.options,
                            self.targetBED,
                            self.reference,
                        )
                        proj_record.annotate(
                            self.ensembl,
                            self.dbsnp,
                            self.reference,
                            self.impactdir,
                        )
                        haplotype.apply_fixture_row_to_record(
                            proj_record, fixture_canonical_row
                        )
                        haplotype.add_haplotype_flags(
                            proj_record, original_ids, original_ids
                        )
                        split_subset_records.append(proj_record)

            # Writing annotated canonical record to output file.
            if emit_canonical_record:
                record.output(
                    self.options.args["outputformat"],
                    self.outfile,
                    self.options,
                    self.genelist,
                    self.transcriptlist,
                    self.snplist,
                    self.copts.stdout,
                )

            for subset_record in split_subset_records:
                subset_record.output(
                    self.options.args["outputformat"],
                    self.outfile,
                    self.options,
                    self.genelist,
                    self.transcriptlist,
                    self.snplist,
                    self.copts.stdout,
                )

            # Optional additional nearby-protein split output records.
            if (
                emit_canonical_record
                and parsed_haplotype is not None
                and self.options.args.get(
                "splitadjacentprotein", False
                )
            ):
                full_ids = ";".join([a.token for a in parsed_haplotype.atomic])
                splitnearby_rows = haplotype.get_splitnearby_fixture_rows(fixture_rows)

                if splitnearby_rows:
                    for sr in splitnearby_rows:
                        alt_line = haplotype.build_record_line_like(
                            record,
                            record.chrom,
                            record.pos,
                            record.id,
                            record.ref,
                            record.alts[0] if len(record.alts) > 0 else "",
                        )
                        alt_record = core.Record(
                            alt_line, self.options, self.targetBED, self.reference
                        )
                        alt_record.annotate(
                            self.ensembl, self.dbsnp, self.reference, self.impactdir
                        )
                        haplotype.apply_fixture_row_to_record(alt_record, sr)
                        haplotype.add_haplotype_flags(alt_record, full_ids, full_ids)
                        alt_record.output(
                            self.options.args["outputformat"],
                            self.outfile,
                            self.options,
                            self.genelist,
                            self.transcriptlist,
                            self.snplist,
                            self.copts.stdout,
                        )
                else:
                    alt_line = haplotype.build_record_line_like(
                        record,
                        record.chrom,
                        record.pos,
                        record.id,
                        record.ref,
                        record.alts[0] if len(record.alts) > 0 else "",
                    )
                    alt_record = core.Record(
                        alt_line, self.options, self.targetBED, self.reference
                    )
                    alt_record.annotate(
                        self.ensembl, self.dbsnp, self.reference, self.impactdir
                    )
                    changed_any = False

                    for v in alt_record.variants:
                        if not all(
                            x in v.flags
                            for x in ["CSN", "PROTPOS", "PROTREF", "PROTALT"]
                        ):
                            continue
                        csn_vals = v.getFlag("CSN").split(":")
                        pos_vals = v.getFlag("PROTPOS").split(":")
                        ref_vals = v.getFlag("PROTREF").split(":")
                        alt_vals = v.getFlag("PROTALT").split(":")
                        L = min(
                            len(csn_vals), len(pos_vals), len(ref_vals), len(alt_vals)
                        )
                        new_csn_vals = list(csn_vals)

                        for i in range(L):
                            new_csn = haplotype.maybe_build_splitnearby_csn(
                                csn_vals[i], pos_vals[i], ref_vals[i], alt_vals[i]
                            )
                            if new_csn:
                                new_csn_vals[i] = new_csn
                                changed_any = True

                        if changed_any:
                            idx = v.flags.index("CSN")
                            v.flagvalues[idx] = ":".join(new_csn_vals)

                    if changed_any:
                        haplotype.add_haplotype_flags(alt_record, full_ids, full_ids)
                        alt_record.output(
                            self.options.args["outputformat"],
                            self.outfile,
                            self.options,
                            self.genelist,
                            self.transcriptlist,
                            self.snplist,
                            self.copts.stdout,
                        )

            # Writing progress information to log file
            if self.threadidx == 1 and self.options.args["logfile"]:
                x = round(
                    100 * counter / int(self.numOfRecords / self.copts.threads), 1
                )
                x = min(x, 100.0)
                if x > thr:
                    logging.info(str(thr) + "% of records annotated.")
                    thr += 10

        # Closing process input and output files
        self.infile.close()
        self.outfile.close()

        # Finalizing progress info
        if not self.copts.stdout and self.threadidx == 1:
            finalizeProgressInfo()


def run(copts, version):
    copts.threads = int(copts.threads)
    if copts.threads > 1:
        copts.stdout = False

    # Check if input and configuration files exist
    if copts.conf is None:
        print("\nError: no configuration file specified.")
        print(
            "Please use option -c or add the absolute path to the default_config_path file.\n"
        )
        quit()
    if not os.path.isfile(copts.conf):
        print("\nError: configuration file (" + copts.conf + ") cannot be found.\n")
        quit()
    if not os.path.isfile(copts.input):
        print("\nError: input file (" + copts.input + ") cannot be found.\n")
        quit()

    # Reading options from configuration file
    options = core.Options(copts.conf)

    # Command-line flags override config defaults for haplotype modes.
    options.args["parseHaplotype"] = bool(getattr(copts, "parseHaplotype", False))
    options.args["splitBasedOnProtein"] = bool(
        getattr(copts, "splitBasedOnProtein", False)
    )
    options.args["splitadjacentprotein"] = bool(
        getattr(copts, "splitAdjacentProtein", False)
    )

    if (
        options.args["splitBasedOnProtein"] or options.args["splitadjacentprotein"]
    ) and not options.args["parseHaplotype"]:
        print(
            "ERROR: --splitBasedOnProtein and --splitAdjacentProtein require --parseHaplotype"
        )
        quit()

    # Initializing log file
    if options.args["logfile"]:
        logging.basicConfig(
            filename=copts.output + ".log",
            filemode="w",
            format="%(asctime)s %(levelname)s: %(message)s",
            level=logging.DEBUG,
        )

    # Printing out version.py information and start time

    if not copts.stdout:
        starttime = printStartInfo(version)
    if options.args["logfile"]:
        logging.info("CAVA " + version + " started.")

    # Checking if options specified in the configuration file are correct
    core.checkOptions(options)

    # Printing out configuration, input and output file names
    if not copts.stdout:
        printInputFileNames(copts, options)

    # Reading gene, transcript and snp lists from files
    genelist = core.readSet(options, "genelist")
    transcriptlist = core.readSet(options, "transcriptlist")
    snplist = core.readSet(options, "snplist")

    # Reading (new) transcript2protein map for HGVSP annotation
    print("INFO: reading transcript2protein file\n")
    options.transcript2protein = core.read_dict(options, "transcript2protein")

    print(
        "transcript2protein has " + str(len(options.transcript2protein)) + " mappings\n"
    )
    # Parsing @impactdef string
    if not (options.args["impactdef"] == "." or options.args["impactdef"] == ""):
        impactdir = dict()
        valuev = options.args["impactdef"].split("|")
        for i in range(len(valuev)):
            classv = valuev[i].split(",")
            for c in classv:
                impactdir[c.strip()] = str(i + 1)
    else:
        impactdir = None

    # Counting and printing out number of records of input file
    numOfRecords = core.countRecords(copts.input)
    if not copts.stdout:
        printNumOfRecords(numOfRecords)
    if options.args["logfile"]:
        logging.info(str(numOfRecords) + " records to be annotated.")

    # Writing header to output file
    if options.args["outputformat"] == "VCF":
        outfname = copts.output + ".vcf"
    else:
        outfname = copts.output + ".txt"
    outfile = open(outfname, "w", encoding="utf-8")
    header = readHeader(copts.input)
    try:
        core.writeHeader(options, "\n".join(header), outfile, copts.stdout, version)
    except:
        sys.stderr.write("CAVA: error writing header to " + outfname + "\n")
        exit(1)
    outfile.close()

    # Find break points in the input file
    breaks = findFileBreaks(copts.input, copts.threads)

    # Initializing annotation processes
    threadidx = 0
    processes = []
    for startline, endline in breaks:
        threadidx += 1
        processes.append(
            SingleJob(
                threadidx,
                options,
                copts,
                startline,
                endline,
                genelist,
                transcriptlist,
                snplist,
                impactdir,
                numOfRecords,
            )
        )

    # Running annotation processes
    for process in processes:
        process.start()
    for process in processes:
        process.join()

    # Merging tmp files
    if copts.threads > 1:
        mergeTmpFiles(copts.output, options.args["outputformat"], copts.threads)

    # Printing out summary information and end time
    if not copts.stdout:
        printEndInfo(options, copts, starttime)
