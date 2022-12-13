"""Hail Batch (https://hail.is/docs/batch/service.html) pipeline that runs imputed_sample_metadata.py in parallel on many samples."""

import hashlib
import os
import pandas as pd
import sys

from step_pipeline import pipeline, Backend, Localize, Delocalize, all_outputs_exist

DOCKER_IMAGE = "weisburd/impute_sample_metadata@sha256:9cd73f4253079782d8a7352bed7088c73c3d868194e8a41f93eb0d9261096e22"

OUTPUT_FILENAME_PREFIX = "imputed_sample_metadata"

HG19_REFERENCE_FASTA = "gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta"
HG38_REFERENCE_FASTA = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"

def parse_args(batch_pipeline):
    """Define and parse command-line args"""

    arg_parser = batch_pipeline.get_config_arg_parser()
    group = arg_parser.add_argument_group("imputed_sample_metadata general settings")
    group.add_argument("-o", "--output-dir", default="gs://bw-proj/impute_sample_metadata",
        help="Output directory where to copy the imputed_sample_metadata output .tsv file(s).")
    group.add_argument("-s", "--sample-id", action="append",
        help="If specified, only this sample id will be processed from the input table (useful for testing).")
    group.add_argument("-n", "--num-samples-to-process", type=int,
        help="If specified, only this many samples will be processed from the input table (useful for testing).")

    group.add_argument("sample_table",
        help="Path of tab-delimited table containing sample ids along with their BAM or CRAM file paths "
             "and other metadata. The table should contain at least the following columns: "
             "'sample_id' or 'individual_id', "
             "'cram_path' or 'bam_path' or 'reads', "
             "'bai_path' or 'crai_path' or 'index'. ")

    group.add("--sample-id-column",
              help="Optionally specify the name of input table column that contains the sample id")
    group.add("--cram-or-bam-path-column",
              help="Optionally specify the name of input table column that contains the CRAM or BAM path")
    group.add("--crai-or-bai-path-column",
              help="Optionally specify the name of input table column that contains the CRAI or BAI path")
    args = batch_pipeline.parse_known_args()

    return args


def parse_sample_table(batch_pipeline):
    """Parse and validate the sample input table which contains paths and metadata for samples to process.

    Return:
        2-tuple (pandas.DataFrame, args): The input table and command line args.
    """
    args = parse_args(batch_pipeline)

    df = pd.read_table(args.sample_table, dtype=str)

    # check table columns
    arg_parser = batch_pipeline.get_config_arg_parser()
    if not args.sample_id_column:
        if "individual_id" in df.columns:
            args.sample_id_column = "individual_id"
        elif "sample_id" in df.columns:
            args.sample_id_column = "sample_id"
        else:
            arg_parser.error(f"{args.sample_table} must have one of these columns: 'individual_id', 'sample_id'")
    elif args.sample_id_column not in df.columns:
        arg_parser.error(f"{args.sample_table} doesn't have a '{args.sample_id_column}' column")

    if not args.cram_or_bam_path_column:
        if "cram_path" in df.columns:
            args.cram_or_bam_path_column = "cram_path"
        elif "bam_path" in df.columns:
            args.cram_or_bam_path_column = "bam_path"
        elif "reads" in df.columns:
            args.cram_or_bam_path_column = "reads"
        else:
            arg_parser.error(f"I{args.sample_table} must have one of these columns: 'cram_path', 'bam_path', 'reads'")
    elif args.cram_or_bam_path_column not in df.columns:
        arg_parser.error(f"{args.sample_table} doesn't have a '{args.cram_or_bam_path_column}' column")

    if not args.crai_or_bai_path_column:
        if "crai_path" in df.columns:
            args.crai_or_bai_path_column = "crai_path"
        elif "bai_path" in df.columns:
            args.crai_or_bai_path_column = "bai_path"
        elif "index" in df.columns:
            args.crai_or_bai_path_column = "index"
        else:
            arg_parser.error(f"{args.sample_table} must have one of these columns: 'crai_path', 'bai_path', 'index'")
    elif args.crai_or_bai_path_column not in df.columns:
        arg_parser.error(f"{args.sample_table} doesn't have a '{args.crai_or_bai_path_column}' column")

    # filter table to rows that have CRAM or BAM paths
    df = df[~df[args.cram_or_bam_path_column].isna() & ~df[args.crai_or_bai_path_column].isna()]
    df = df.drop_duplicates(subset=[args.cram_or_bam_path_column, args.crai_or_bai_path_column])
    df = df.sort_values(args.sample_id_column)

    # apply --sample-id and --num-samples-to-process args if they were specified
    if args.sample_id:
        df = df[df[args.sample_id_column].isin(args.sample_id)]
        print(f"Found {len(df)} out of {len(args.sample_id)} requested sample ids")
        if len(df) == 0:
            sys.exit(1)
        elif len(df) < len(args.sample_id):
            print(f"WARNING: Couldn't find sample ids: {set(args.sample_id) - set(df.sample_id)}")

    if args.num_samples_to_process:
        df = df.iloc[:args.num_samples_to_process]

    print(f"Parsed {len(df)} rows from {args.sample_table}")

    return df, args


def main():
    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE)

    df, args = parse_sample_table(bp)

    bp.set_name(f"imputed_sample_metadata: {len(df)} samples")

    # compute a hash of the sample ids being processed
    analysis_id = ", ".join(df[args.sample_id_column])
    analysis_id = hashlib.md5(analysis_id.encode('UTF-8')).hexdigest().upper()
    analysis_id = analysis_id[:10]  # shorten

    if not args.force:
        bp.precache_file_paths(os.path.join(args.output_dir, "**", f"{OUTPUT_FILENAME_PREFIX}*"))

    steps = []
    print(f"Processing {len(df)} samples")
    for _, row in df.iterrows():
        row_sample_id = row[args.sample_id_column]
        row_cram_or_bam_path = row[args.cram_or_bam_path_column]
        row_crai_or_bai_path = row[args.crai_or_bai_path_column]

        # step1: run imputed_sample_metadata.py
        s1 = bp.new_step(
            f"Impute Sample Metadata pipeline: {row_sample_id}",
            arg_suffix="step1",
            image=DOCKER_IMAGE,
            cpu=0.25,
            memory="standard",
            delocalize_by=Delocalize.COPY,
        )
        s1.switch_gcloud_auth_to_user_account()
        s1.regions("us-central1")
        s1.command("set -euxo pipefail")

        # process input files
        hg19_fasta = s1.input(HG19_REFERENCE_FASTA, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)
        hg38_fasta = s1.input(HG38_REFERENCE_FASTA, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)
    
        #cram_or_bam_input = s1.input(row_cram_or_bam_path, localize_by=Localize.HAIL_BATCH_CLOUDFUSE_VIA_TEMP_BUCKET)  # HAIL_BATCH_CLOUDFUSE_VIA_TEMP_BUCKET
        #crai_or_bai_input = s1.input(row_crai_or_bai_path, localize_by=Localize.HAIL_BATCH_CLOUDFUSE_VIA_TEMP_BUCKET)
        cram_or_bam_input = s1.input(row_cram_or_bam_path, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)  # HAIL_BATCH_CLOUDFUSE_VIA_TEMP_BUCKET
        crai_or_bai_input = s1.input(row_crai_or_bai_path, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)

        s1.command(f"ls -lh {cram_or_bam_input}")

        # create symlinks in the same directory to work around cases when they are in different directories in the cloud
        s1.command(f"ln -s {cram_or_bam_input} /{cram_or_bam_input.filename}")
        s1.command(f"ln -s {crai_or_bai_input} /{crai_or_bai_input.filename}")

        output_tsv_path = f"/io/{OUTPUT_FILENAME_PREFIX}.{row_sample_id}.tsv"
        s1.command("cd /")
        s1.command(f"python3 /impute_sample_metadata.py /{cram_or_bam_input.filename} "
                   f"--hg19-fasta {hg19_fasta} "
                   f"--hg38-fasta {hg38_fasta} "
                   f"--verbose "
                   f"-o {output_tsv_path}")

        # add cram path column
        s1.command(f"""python3 << CODE
import pandas as pd
df = pd.read_table('{output_tsv_path}')
df.loc[:, 'cram_path'] = '{row_cram_or_bam_path}'
df.to_csv('{output_tsv_path}', sep='\t', index=False, header=True)
CODE
""")
        s1.command(f"cat {output_tsv_path}")
        s1.command("ls")

        # delocalize the output tsv
        destination_dir = os.path.join(args.output_dir, os.path.dirname(row_cram_or_bam_path).replace("gs://", ""), os.path.basename(output_tsv_path))
        s1.output(output_tsv_path, destination_dir)
        steps.append(s1)

    # step2: combine tables from step1 into a single table
    s2 = bp.new_step(
        f"Combine {len(steps)} tables",
        image=DOCKER_IMAGE,
        cpu=1,
        memory="standard",
        output_dir=args.output_dir,
        delocalize_by=Delocalize.COPY,
        arg_suffix="step2",
    )
    s2.regions("us-central1")
    s2.command("set -euxo pipefail")

    combined_output_tsv_filename = f"combined_results.{len(df)}_samples.{analysis_id}.tsv"
    for i, step in enumerate(steps):
        #if args.skip_step1 and not all_outputs_exist(step):
        #    print(f"WARNING: skipping {step}")
        #    continue
        s2.depends_on(step)
        tsv_input = s2.use_previous_step_outputs_as_inputs(step, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)
        if i == 0:
            s2.command(f"head -n 1 {tsv_input} > {combined_output_tsv_filename}")
        s2.command(f"tail -n +2 {tsv_input} >> {combined_output_tsv_filename}")

    s2.command(f"gzip {combined_output_tsv_filename}")
    combined_output_tsv_filename = f"{combined_output_tsv_filename}.gz"

    s2.output(combined_output_tsv_filename, delocalize_by=Delocalize.COPY)

    bp.run()

    # download the output table from step2 and merge it with the input table given to this pipeline.
    os.system(f"gsutil -m cp {os.path.join(args.output_dir, combined_output_tsv_filename)} .")
    result_df = pd.read_table(combined_output_tsv_filename)
    result_df.loc[:, "sample_id_or_filename_prefix"] = result_df.sample_id.where(
        result_df.sample_id.isin(set(df[args.sample_id_column])), result_df.filename_prefix)

    df = df.drop_duplicates(subset=[args.sample_id_column], keep="first")
    df_with_metadata = pd.merge(result_df, df, how="left", left_on="cram_path", right_on=args.cram_or_bam_path_column)
    df_with_metadata.to_csv(combined_output_tsv_filename, sep="\t", header=True, index=False)
    print(f"Wrote {len(df_with_metadata)} rows to {combined_output_tsv_filename}")


if __name__ == "__main__":
    main()
