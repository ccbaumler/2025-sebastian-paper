#! /usr/bin/env python

import csv
import argparse

def load_blast_ids(csv_file):
    idents = set()
    with open(csv_file, newline='') as fp:
        reader = csv.reader(fp)
        for row in reader:
            if len(row) >= 2:
                idents.add(row[1].strip())
                print(row)
                print(row[1])
    return idents

def load_cds_ids(header_line):
    for h in header_line.split():
        return h.strip('>')


def extract_sequences(fasta_file, idents, output):
    with open(fasta_file) as fasta_in, open(output, 'w') as fasta_out:
        write = False
        header = ""
        sequence = []

        for line in fasta_in:
            if line.startswith(">"):
                if write and header and sequence:
                    fasta_out.write(header)
                    fasta_out.writelines(sequence)

                header = line
                cds_id = load_cds_ids(header)

                sequence = []
                write = cds_id in idents
            else:
                if write:
                    sequence.append(line)

        if write and header and sequence:
            fasta_out.write(header)
            fasta_out.writelines(sequence)

if __name__ == "__main__":
    p = argparse.ArgumentParser()

    p.add_argument('blast_csv', help='The blast output in csv format (outfmt 10)')
    p.add_argument('cds_fasta', help='The concatinated Coding Sequence Fasta file from NCBI')
    p.add_argument('output', help='The output')

    args = p.parse_args()

    idents = load_blast_ids(args.blast_csv)
    extract_sequences(args.cds_fasta, idents, args.output)

