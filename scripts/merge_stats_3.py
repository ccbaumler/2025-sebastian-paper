#! /usr/bin/env python

import polars as pl
import argparse
import os

def main():
    p = argparse.ArgumentParser(description="Merge mapping and coverage stats by contig.")
    p.add_argument(
        "run_accession",
        help="Run accession used in filename (e.g. SRR1211428) without suffix. Will look for *_mapping_stats.tsv and *_average_coverage_stats.tsv"
    )
    p.add_argument(
        "-d", "--directory", default=None,
        help="The directory to look for *_mapping_stats.tsv and *_average_coverage_stats.tsv"
    )
    p.add_argument(
        "-o", "--output", default=None,
        help="Output file path. If not provided, prints to stdout."
    )
    p.add_argument(
        "-f", "--filter", action='store_true',
        help="Filter the columns containing no information out of Output file."
    )

    args = p.parse_args()
    run_a = args.run_accession

    if args.directory:
        try:
            base = os.path.join(args.directory, run_a)
        except ValueError:
            raise ValueError("Basename must start with accession E.g. SRRxxxxxxx")
    else:
        print("No directory argument used. Search current directory...")
        try:
            base = run_a
        except ValueError:
            raise ValueError("Basename must be of the form SRRxxxxxxx")
 
    map_file = f"{base}_mapping_stats.tsv"
    cov_file = f"{base}_average_coverage_stats.tsv"

    if not os.path.isfile(map_file):
        raise FileNotFoundError(f"Mapping stats file not found: {map_file}")
    if not os.path.isfile(cov_file):
        raise FileNotFoundError(f"Coverage stats file not found: {cov_file}")

    mapping_stats = pl.read_csv(map_file, separator="\t")
    coverage_stats = pl.read_csv(cov_file, separator="\t")

    df = mapping_stats.join(coverage_stats, on="contig", how="inner")

    df = df.with_columns([
        pl.lit(run_a).alias("run"),
    ])

    df = df.select(["contig", "run"] + [col for col in df.columns if col not in {"contig", "run"}])

    if df.is_empty():
        # Create a new DataFrame with one row
        new_row = pl.DataFrame({
            "contig": ["None"],
            "run": [run_a],
            **{col: [0] for col in df.columns if col not in {"contig", "run"}}
        })

        # Assign the new DataFrame
        df = new_row
    df = df.drop('length_right')
    print(df)
    cols_to_check = [col for col in df.columns if col not in {"contig", "run", "length"}]
    condition = ~pl.all_horizontal([pl.col(c) == 0 for c in cols_to_check])
    df_filter = df.filter(condition)
    print(df_filter)

    if args.output:
        if args.filter:
            df_filter.write_csv(args.output, separator=",")
        else:
            df.write_csv(args.output, separator=",")
        print(f"Output written to {args.output}")
    else:
        print(df)


if __name__ == "__main__":
    main()

