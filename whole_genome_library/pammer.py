from Bio import SeqIO
import pandas as pd
import argparse
import os

def extract_pam(row, genome_record):
    start, end = map(int, row['coord'].split('-'))
    if row['q_dir'] == 'R':
        pam = genome_record.seq[start-4:start-1]
    else:
        pam = genome_record.seq[end+1:end+4].reverse_complement()
    return str(pam)

def main(genbank_file, tsv_file):
    genome_record = SeqIO.read(genbank_file, "genbank")
    data = pd.read_csv(tsv_file, sep='\t')
    data['pam'] = data.apply(extract_pam, axis=1, genome_record=genome_record)
    output_file = os.path.splitext(tsv_file)[0] + "_with_pam.tsv"
    data.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract PAM sequences from genomic coordinates")
    parser.add_argument("genbank_file", help="Path to the GenBank file")
    parser.add_argument("tsv_file", help="Path to the TSV file containing genomic coordinates")
    args = parser.parse_args()
    main(args.genbank_file, args.tsv_file)
