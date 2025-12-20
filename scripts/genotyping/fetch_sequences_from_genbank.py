import argparse
from Bio import Entrez
import pandas as pd

def fetch_sequences(file_path, output_fasta, email):
    df = pd.read_excel(file_path)
    accession_list = df['Accession no.'].dropna().tolist()

    Entrez.email = email

    with open(output_fasta, "w") as fasta_file:
        for accession in accession_list:
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
            record = handle.read()
            handle.close()

            fasta_file.write(record)
            print(f"Fetched and saved sequence for accession: {accession}")

    print(f"All sequences have been fetched and saved to: {output_fasta}")

def main():
    parser = argparse.ArgumentParser(
        description='Fetch sequences from NCBI using accession numbers from an Excel file and save them to a FASTA file.'
    )
    parser.add_argument('-f', '--file', required=True, help='Path to the Excel file containing accession numbers.')
    parser.add_argument('-o', '--output', required=True, help='Path to the output FASTA file where sequences will be saved.')
    parser.add_argument('-e', '--email', required=True, help='Your email address (required by NCBI Entrez).')

    args = parser.parse_args()

    fetch_sequences(args.file, args.output, args.email)

if __name__ == '__main__':
    main()
