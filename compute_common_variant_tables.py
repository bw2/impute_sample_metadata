"""This hail query script computes lists of common variants for hg38 and hg19 in exome and genome callsets.
These variants are later used to impute the reference version and sample type (exome or genome) of unknown samples.
"""

import os
import hail as hl
hl.init(log="/dev/null")


OUTPUT_DIR = "gs://bw-proj/temp/"

def get_sites(ht, label, previous_sites, max_loci=1000):
    """Takes a hail table and returns a list of sites that are common variants and are not in in the previous_sites
    list.

    Args:
        ht (HailTable): A hail table with 4 fields: locus, alleles, AF, AN
        label (str): A name used for the checkpoint filename
        previous_sites: A collection of sites to exclude from the results list.
        max_loci (int): The maximum number of loci to return.
    Return:
        list: strings representing common variants like ["chrX-12345-A-G", ...]
    """
    #ht = ht.filter(ht.locus.in_x_nonpar(), keep=True)

    ht = ht.checkpoint(f"gs://bw2-delete-after-5-days/{label}.ht", overwrite=False, _read_if_exists=True)

    mean_AN = ht.aggregate(hl.agg.mean(ht.AN))
    print(f"Mean AN: {mean_AN}")
    ht = ht.filter(
        (ht.AF > 0.3) &
        (ht.AF < 0.7) &
        (ht.AN > 0.5 * mean_AN), keep=True)

    print("Found", ht.count(), "common variant sites in", label)

    ht = ht.select(variant=hl.str("-").join([
        ht.locus.contig, hl.str(ht.locus.position), ht.alleles[0], ht.alleles[1],
    ]))
    if len(previous_sites) > 0:
        ht = ht.filter(hl.set(previous_sites).contains(ht.variant), keep=False)
        print(ht.count(), "total sites in", label, f"that are not in the {len(previous_sites)} previous sites.")

    results = ht.variant.take(max_loci)
    print(f"Downloaded", len(results), "sites from", label)

    return results


def main():
    sites = {}
    previous_sites = []

    for label, path in [
        ("hg38_exomes", "gs://seqr-datasets/v02/GRCh38/RDG_WES_Broad_Internal/v16/RDG_WES_Broad_Internal.mt"),
        ("hg38_genomes", "gs://seqr-datasets/v02/GRCh38/RDG_WGS_Broad_Internal/v30/RDG_WGS_Broad_Internal.mt"),
        ("hg19_exomes", "gs://gcp-public-data--gnomad/release/2.1.1/ht/exomes/gnomad.exomes.r2.1.1.sites.ht"),
        ("hg19_genomes", "gs://gcp-public-data--gnomad/release/2.1.1/ht/genomes/gnomad.genomes.r2.1.1.sites.ht"),
    ]:
        print(f"Reading {path}")
        if path.endswith(".ht"):
            ht = hl.read_table(path)
            ht = ht.key_by()
            ht = ht.select(AF=ht.freq.AF[0], AN=ht.freq.AN[0], locus=ht.locus, alleles=ht.alleles)
        elif path.endswith(".mt"):
            ht = hl.read_matrix_table(path).rows()
            ht = ht.key_by()
            ht = ht.select(AF=ht.AF, AN=ht.AN, locus=ht.locus, alleles=ht.alleles)

        sites[label] = get_sites(ht, label, previous_sites)
        previous_sites.extend(sites[label])

        output_path = os.path.join(OUTPUT_DIR, f"{label}.common_variants.tsv")
        with hl.hadoop_open(output_path, "w") as f:
            for variant in sites[label]:
                f.write(f"{variant}\n")
        print(f"Wrote {len(sites[label])} to {output_path}")


if __name__ == "__main__":
    main()