#!/usr/env python3

"""This script takes a CRAM or BAM file and imputes the reference genome and sample type (exome or genome)
by computing the genotypes of common SNVs on the X chromosome.
"""
import argparse
import collections
import hashlib
import os
import pysam
import pandas as pd
import re
import subprocess
import tqdm


NUM_VARIANTS_NEEDED_FOR_IMPUTATION = 100

COMMON_VARIANT_POSITIONS = {
    "hg19_exome":  "hg19_exome.common_variants.tsv",
    "hg19_genome": "hg19_genome.common_variants.tsv",
    "hg38_exome":  "hg38_exome.common_variants.tsv",
    "hg38_genome": "hg38_genome.common_variants.tsv",
}


def parse_args():
    """Define and then parse command-line args"""

    parser = argparse.ArgumentParser()

    g = parser.add_argument_group(title="cram parsing", description="A reference fasta is often necessary for "
        "parsing CRAMs, so in that case we need to first parse the reference version from the CRAM file header.")
    g.add_argument("--hg19-fasta", help="hg19 fasta path")
    g.add_argument("--hg38-fasta", help="hg38 fasta path")

    parser.add_argument("-o", "--output-tsv", help="Optional output tsv file path")
    parser.add_argument("-p", "--show-progress-bar", action="store_true", help="Whether to show a progress bar")
    parser.add_argument("-v", "--verbose", action="store_true", help="Whether to print extra details during the run")
    parser.add_argument("cram_or_bam_path", nargs="+", help="One or more CRAM or BAM file paths")

    args = parser.parse_args()

    # make sure the input files exist
    for cram_or_bam_path in args.cram_or_bam_path:
        if not os.path.isfile(cram_or_bam_path):
            parser.error(f"File not found: {cram_or_bam_path}")

    # define the output_tsv if it wasn't specified
    if not args.output_tsv:
        if len(args.cram_or_bam_path) > 1:
            input_paths_hash = hashlib.sha256("|".join(args.cram_or_bam_path).encode('utf-8')).hexdigest().upper()
            args.output_tsv = f"imputed_sample_metadata.{len(args.cram_or_bam_path)}_samples.{input_paths_hash[:10]}.tsv"
        else:
            filename_prefix, _ = get_filename_prefix_and_file_type(args.cram_or_bam_path[0])
            args.output_tsv = f"{filename_prefix}.imputed_sample_metadata.tsv"

    if args.verbose:
        print("Input args:")
        print(f"    --output-tsv: {os.path.abspath(args.output_tsv)}")
        print(f"    CRAMS or BAMS:", ", ".join(map(os.path.abspath, args.cram_or_bam_path)))

    return args


def get_genome_version_and_sample_type_from_bam_or_cram_header(bam_or_cram_path, verbose=False):
    # get genome version from file header
    header_text = subprocess.check_output(
        f"samtools view -H {bam_or_cram_path}", shell=True, encoding="UTF-8", stderr=subprocess.DEVNULL)
    genome_version = None
    if "AS:GRCh37" in header_text or "Homo_sapiens_assembly19.fasta" in header_text:
        genome_version = 19
    elif "AS:GRCh38" in header_text or "Homo_sapiens_assembly38.fasta" in header_text:
        genome_version = 38
    else:
        if verbose:
            print(f"WARNING: unable to determine genome version from {bam_or_cram_path} file header: {header_text}")

    header_text = header_text.lower()
    if "exome" in header_text:
        sample_type = "exome"
    elif "genome" in header_text:
        sample_type = "genome"
    else:
        sample_type = None
    
    return genome_version, sample_type


def count_nucleotides_at_position(alignment_file, chrom, pos_1based):
    """Count the number of A, C, G, and T's found at the given genomic position within the given read data.

    Args:
        alignment_file (pysam.AlignmentFile): The pysam AlignmentFile object representing the input BAM or CRAM file.
        chrom (str): Chromosome of the genomic position.
        pos_1based (int): 1-based genomic position where to count nucleotides.

    Return:
        dict: The keys are nucleotides "A", "C", "G", "T", and the values are counts representing the number of times
            a read contained that nucleotide at the given position.
    """

    pos_0based = pos_1based - 1
    nucleotide_counts = {"A": 0, "C": 0, "G": 0, "T": 0}

    # handle chr prefix
    chrom = chrom.replace("chr", "")
    if f"chr{chrom}" in alignment_file.references:
        chrom = f"chr{chrom}"
    elif chrom not in alignment_file.references:
        raise ValueError(f"Chromosome not found: {chrom}")

    for pileup_column in alignment_file.pileup(
            region=f"{chrom}:{pos_0based}-{pos_1based}",
            stepper="samtools",
            ignore_overlaps=True,
            ignore_orphans=False,
            min_mapping_quality=0,
            min_base_quality=13,
            truncate=True,
            multiple_iterators=False):

        if pileup_column.pos < pos_0based:
            continue

        if pileup_column.pos != pos_0based:
            raise ValueError(f"Unexpected pileup position: {chrom}:{pileup_column.pos}. "
                             f"Expecting {chrom}:{pos_0based}")

        # iterate over the reads in the pileup
        for base in pileup_column.get_query_sequences():
            if not base:
                continue

            base = base.upper()
            if base in nucleotide_counts:
                nucleotide_counts[base] += 1
            else:
                print(f"WARNING: Unexpected base '{base}' found at {chrom}:{pos_1based}")
        break

    return nucleotide_counts


def get_filename_prefix_and_file_type(cram_or_bam_path):
    """Returns the filename prefix and file type suffix.

    Args:
        cram_or_bam_path (str): Input CRAM or BAM path.
        output_row (dict): dictionary of fields to be written to the .tsv output row for the current sample
    """
    filename_pattern_match = re.match("(.*)[.](cram|bam)$", os.path.basename(cram_or_bam_path))
    if not filename_pattern_match:
        raise ValueError(f"File path doesn't end with have '.cram' or '.bam': {cram_or_bam_path}")

    filename_prefix = filename_pattern_match.group(1)
    file_type = filename_pattern_match.group(2)

    return filename_prefix, file_type


def set_sample_id(alignment_file, output_row):
    """Try reading the sample_id from the alignment file header by looking for a read group (@RG) with a sample (SM)
    field.

    Args:
        alignment_file (pysam.AlignmentFile): The pysam AlignmentFile object representing the input BAM or CRAM file.
        output_row (dict): dictionary of fields to be written to the .tsv output row for the current sample
    """
    output_row["sample_id"] = ""
    for read_group in alignment_file.header.get("RG", []):
        if "SM" in read_group:
            output_row["sample_id"] = read_group["SM"]
            break


def compute_het_hom_or_missing(ref_allele, nucleotide_counts):
    total = sum(nucleotide_counts.values())
    if total < 5:
        return "MISSING"

    alleles_detected = [base for base, count in nucleotide_counts.items() if count >= 2]
    if len(alleles_detected) == 1:
        if alleles_detected[0].upper() == ref_allele.upper():
            return "HOM_REF"
        else:
            return "HOM_ALT"
    elif len(alleles_detected) > 1:
        return "HET"
    else:
        raise ValueError(f"Unexpected number of alleles detected: {alleles_detected}: {nucleotide_counts}")


def impute_metadata(output_row):
    """Compute the reference genome version, whether it's an exome or a genome , and whether it's a male or a female.

    Args:
        output_row (dict): dictionary of fields to be written to the .tsv output row for the current sample
    """
    for column in "imputed_reference_genome", "imputed_sample_type":  # ,  "imputed_sex"
        if not output_row.get(column):
            output_row[column] = "unknown"   # "hg19" or "hg38"

    if output_row[f"hg38_exome:HET_OR_HOM_ALT"] > 0.15 * NUM_VARIANTS_NEEDED_FOR_IMPUTATION:
        genome_version = output_row["imputed_reference_genome"] = "hg38"   # "hg19" or "hg38"
    else:
        genome_version = output_row["imputed_reference_genome"] = "hg19"   # "hg19" or "hg38"

    if output_row["reference_genome_in_header"] != "unknown" and output_row["imputed_reference_genome"] != "unknown":
        if output_row["reference_genome_in_header"] != output_row["imputed_reference_genome"]:
            print(f"WARNING: reference genome in file header (", output_row["reference_genome_in_header"],
                ") doesn't match imputed reference genome based on common SNVs (",
                output_row["reference_genome_in_header"], ").")
        genome_version = output_row["reference_genome_in_header"]
    else:
        print("WARNING: unable to impute reference genome version")
        return

    if output_row[f"{genome_version}_genome:MISSING"] > 0.5 * NUM_VARIANTS_NEEDED_FOR_IMPUTATION:
        sample_type = "exome"
    else:
        sample_type = "genome"

    output_row["imputed_sample_type"] = sample_type

    #if output_row[f"{genome_version}_exome:HET"] < 0.15 * NUM_VARIANTS_NEEDED_FOR_IMPUTATION:
    #    output_row["imputed_sex"] = "male"
    #elif output_row[f"{genome_version}_exome:HET"] > 0.25 * NUM_VARIANTS_NEEDED_FOR_IMPUTATION:
    #    output_row["imputed_sex"] = "female"
    #else:
    #    print("WARNING: unable to impute sample sex")
    #    return


def main():
    args = parse_args()

    for label, filename in COMMON_VARIANT_POSITIONS.items():
        COMMON_VARIANT_POSITIONS[label] = []
        with open(filename, "rt") as f:
            for line in f:
                COMMON_VARIANT_POSITIONS[label].append(line.strip())

    # process the input BAM or CRAM files
    output_rows = []
    for cram_or_bam_path in args.cram_or_bam_path:
        output_row = {}
        output_row["filename_prefix"], output_row["file_type"] = get_filename_prefix_and_file_type(cram_or_bam_path)

        genome_version, sample_type = get_genome_version_and_sample_type_from_bam_or_cram_header(
            cram_or_bam_path, verbose=args.verbose)

        if genome_version is not None:
            print(f"Parsed genome version hg{genome_version} from file header: {cram_or_bam_path}")
            output_row["reference_genome_in_header"] = f"hg{genome_version}"
        else:
            print(f"Couldn't parse genome version from file header: {cram_or_bam_path}")
            output_row["reference_genome_in_header"] = f"unknown"

        reference_fasta = None
        if genome_version is not None and cram_or_bam_path.endswith(".cram") and (args.hg19_fasta or args.hg38_fasta):
            if genome_version == 38:
                reference_fasta = args.hg38_fasta
            elif genome_version == 19:
                reference_fasta = args.hg19_fasta

        if sample_type is not None:
            print(f"Parsed sample_type '{sample_type}' from file header: {cram_or_bam_path}")
            output_row["sample_type_in_header"] = sample_type
        else:
            print(f"Couldn't parse sample_type from file header: {cram_or_bam_path}")
            output_row["sample_type_in_header"] = f"unknown"    
            
        sample_genotype_counters = collections.defaultdict(int)
        print(f"Reading {cram_or_bam_path} (using reference: {reference_fasta})")
        with pysam.AlignmentFile(cram_or_bam_path, reference_filename=reference_fasta) as alignment_file:
            set_sample_id(alignment_file, output_row)

            for label, variants in COMMON_VARIANT_POSITIONS.items():
                variants = variants[:NUM_VARIANTS_NEEDED_FOR_IMPUTATION]
                if args.show_progress_bar:
                    variants = tqdm.tqdm(variants, unit=" common variants")
                for variant in variants:
                    chrom, pos, ref, alt = variant.split("-")
                    nucleotide_counts = count_nucleotides_at_position(alignment_file, chrom, int(pos))
                    genotype = compute_het_hom_or_missing(ref, nucleotide_counts)
                    sample_genotype_counters[f"{label} {genotype}"] += 1

        for label in COMMON_VARIANT_POSITIONS:
            for genotype in "MISSING", "HET", "HOM_REF", "HOM_ALT":
                output_row[f"{label}:{genotype}"] = sample_genotype_counters[f"{label} {genotype}"]
            output_row[f"{label}:HET_OR_HOM_ALT"] = output_row[f"{label}:HET"] + output_row[f"{label}:HOM_ALT"]
            output_row[f"{label}:MISSING_OR_HOM_REF"] = output_row[f"{label}:MISSING"] + output_row[f"{label}:HOM_REF"]

        impute_metadata(output_row)

        output_rows.append(output_row)

    # write results to .tsv
    df = pd.DataFrame(output_rows)
    if args.verbose:
        for i, (_, row) in enumerate(df.iterrows()):
            print("----")
            print(f"Output row #{i+1}:")
            for column in sorted(df.columns, reverse=True):
                print(f"        {column:35s} {row[column]}")

    df[sorted(df.columns, reverse=True)].to_csv(args.output_tsv, sep='\t', header=True, index=False)
    print(f"Wrote {len(output_rows)} rows to {os.path.abspath(args.output_tsv)}")


if __name__ == "__main__":
    main()
