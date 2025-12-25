import argparse
import pandas as pd
from Bio import SeqIO

def extract_and_process_peptides(genbank_file, fasta_file, peptide_output_file, gapped_output_file):
    # Extract mature peptide coordinates from the GenBank file
    records = SeqIO.parse(genbank_file, "genbank")
    mat_peptides = []

    for record in records:
        for feature in record.features:
            if feature.type == "mat_peptide":
                location = feature.location
                product = feature.qualifiers.get('product', [''])[0]
                mat_peptides.append({
                    'Start': location.start + 1,
                    'End': location.end,
                    'Product': product
                })

    peptide_coordinates = pd.DataFrame(mat_peptides)
    peptide_coordinates.to_csv(peptide_output_file, index=False)
    print(f"peptide coordinates saved to: {peptide_output_file}")

    # Get the first sequence from the FASTA file
    def get_first_sequence(fasta_file):
        with open(fasta_file, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                return str(record.seq)

    reference_sequence = get_first_sequence(fasta_file)

    def find_gapped_coordinates(sequence, ungapped_coordinates):
        gapped_positions = []
        ungapped_index = 0

        for i, letter in enumerate(sequence):
            if letter != "-":
                ungapped_index += 1
                if ungapped_index in ungapped_coordinates:
                    gapped_positions.append(i + 1)

        if not gapped_positions:
            return None, None

        return gapped_positions[0], gapped_positions[-1]

    def process_table(table, sequence):
        results = []

        for _, row in table.iterrows():
            start = row['Start']
            end = row['End']
            product = row['Product']

            ungapped_coordinates = range(start, end + 1)
            gapped_start, gapped_end = find_gapped_coordinates(sequence, ungapped_coordinates)

            results.append({
                'Start': gapped_start,
                'End': gapped_end,
                'Product': product
            })

        return pd.DataFrame(results)

    results = process_table(peptide_coordinates, reference_sequence)
    results.to_csv(gapped_output_file, index=False)
    print(f"Gapped peptide coordinates saved to: {gapped_output_file}")

def main():
    parser = argparse.ArgumentParser(description='Extract mat_peptide information from a GenBank file and map gapped coordinates to a reference sequence.')
    parser.add_argument('-g', '--genbank', required=True, help='Input GenBank file containing mat_peptide features.')
    parser.add_argument('-f', '--fasta', required=True, help='Input FASTA file with the reference sequence.')
    parser.add_argument('-p', '--peptide_output', required=True, help='Output CSV file to save peptide coordinates.')
    parser.add_argument('-o', '--gapped_output', required=True, help='Output CSV file to save gapped peptide coordinates.')

    args = parser.parse_args()

    extract_and_process_peptides(args.genbank, args.fasta, args.peptide_output, args.gapped_output)

if __name__ == '__main__':
    main()
