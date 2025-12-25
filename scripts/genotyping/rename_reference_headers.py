import pandas as pd
from Bio import SeqIO
import argparse

def load_annotations_from_excel(excel_file):
    df = pd.read_excel(excel_file)

    df = df.dropna(subset=['Serotype', 'Topotype', 'Isolate name', 'Accession no.', 'Lineage'])

    annotations = {}
    for _, row in df.iterrows():
        accession = row['Accession no.'].strip()
        serotype = row['Serotype'].strip()
        topotype = row['Topotype'].strip() if pd.notnull(row['Topotype']) else '-'
        lineage = row['Lineage'].strip() if pd.notnull(row['Lineage']) else '-'

        new_header = f"{accession}/{serotype}/{topotype}/{lineage}".replace(" ", "")
        annotations[accession] = new_header

    return annotations


def update_fasta_headers(fasta_file, annotations, output_file):
    with open(output_file, 'w') as fout:
        for record in SeqIO.parse(fasta_file, "fasta"):
            accession = record.id.split('.')[0]  # Extract accession number (before the dot)
            if accession in annotations:
                new_header = annotations[accession]
                record.id = new_header
                record.description = ""
            SeqIO.write(record, fout, "fasta")

def main():
    parser = argparse.ArgumentParser(description="Update FASTA headers based on Excel annotations.")
    parser.add_argument('-f', "--fasta", required=True, help="Path to the input FASTA file.")
    parser.add_argument('-t', "--excel", required=True, help="Path to the Excel file containing the annotations.")
    parser.add_argument('-o', "--output", required=True, help="Path to the output FASTA file.")
    
    args = parser.parse_args()
    
    # Load annotations from the Excel file
    annotations = load_annotations_from_excel(args.excel)
    
    # Update FASTA headers
    update_fasta_headers(args.fasta, annotations, args.output)
    
    print(f"Headers successfully updated in {args.output}")

if __name__ == "__main__":
    main()
