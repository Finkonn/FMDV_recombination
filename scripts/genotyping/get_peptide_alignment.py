import pandas as pd
import argparse
from Bio import SeqIO

def extract_sequences(input_fasta, input_table, product):
    gapped_coords = pd.read_csv(input_table)

    row = gapped_coords[gapped_coords['Product'] == product].iloc[0]
    start_coordinate = int(row['Start'])
    end_coordinate = int(row['End'])

    output_fasta = f"{product}.fasta"

    with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            extracted_seq = record.seq[start_coordinate - 1:end_coordinate]
            outfile.write(f">{record.id}\n")
            outfile.write(f"{extracted_seq}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract sequences from a FASTA file based on coordinates from a CSV file for a specified product.')
    parser.add_argument('-f', '--input_fasta', required=True, help='Path to the input FASTA file')
    parser.add_argument('-t', '--input_table', required=True, help='Path to the CSV file containing the coordinates')
    parser.add_argument('-p', '--product', default="VP1", help='Product name to extract sequences for')

    args = parser.parse_args()

    extract_sequences(args.input_fasta, args.input_table, args.product)
