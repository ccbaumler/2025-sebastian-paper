#! /usr/bin/env python

import polars as pl
import argparse
import os


def check_file(file):
    if not os.path.isfile(file):
        raise FileNotFoundError(f"Mapping stats file not found: {file}")
    return file

def fix_numeric_col_type(df):
    # Attempt to cast all columns that look numeric to Float64
    for col in df.columns:
        try:
            df = df.with_column(pl.col(col).cast(pl.Float64, strict=False))
        except:
            pass  # If it can't cast, leave as is
    return df

def main():
    p = argparse.ArgumentParser(description="Merge mapping and coverage stats by contig.")
    p.add_argument(
        "-m", "--metadata-file",
        help="The metadata file containing the diagnosis and read_count metrics"
    )
    p.add_argument(
        "-i", "--input-directory", default=None,
        help="The directory that contains the csv files of the merged_stats.py"
    )
    p.add_argument(
        "-o", "--output", default=None,
        help="Output file path. If not provided, prints to stdout."
    )

    args = p.parse_args()
    indir = args.input_directory

    mf = pl.read_csv(args.metadata_file, separator='\t')

    csv_files = [os.path.join(indir, f) 
                 for f in os.listdir(indir) 
                 if f.endswith(".csv")]

    if not csv_files:
        print("No CSV files found in the directory.")
        return

    dfs = []
    for fp in csv_files:
        df = pl.read_csv(check_file(fp))
        if df.is_empty() or df.height == 0:
            print(f"Skipping empty file: {fp}")
            continue
        df = fix_numeric_col_type(df)
        dfs.append(df)

    if not dfs:
        print("No non-empty CSV files found.")
        return

    contigs_df = pl.concat(dfs)
    print(contigs_df)

    print(mf)
    mf_less = mf.select(["run", "diagnosis", "read_count"])
    joined_df = contigs_df.join(mf_less, on="run", how="left")
    df = joined_df.select(["contig", "run", "diagnosis", "read_count"] + [col for col in joined_df.columns if col not in {"contig", "run", "diagnosis", "read_count"}])

    if args.output:
        df.write_csv(args.output, separator=",")
        print(f"Output written to {args.output}")
    else:
        print(df)

if __name__ == "__main__":
    main()

