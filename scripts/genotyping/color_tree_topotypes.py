import argparse
import pandas as pd
import re
import matplotlib.colors as mcolors

def generate_colors(n):
    return [
        mcolors.to_hex(c)
        for c in mcolors.hsv_to_rgb([(i / n, 1.0, 1.0) for i in range(n)])
    ]

def color_nexus(classification_xlsx, input_nexus, output_nexus, max_colors=58):
    classification_df = pd.read_excel(classification_xlsx)
    classification_df = classification_df[classification_df['Accession no.'].notna()]

    unique_combinations = (
        classification_df[['Serotype', 'Topotype']]
        .drop_duplicates()
        .reset_index(drop=True)
    )

    colors = generate_colors(max_colors)

    combination_to_color = {
        (row['Serotype'], row['Topotype']): colors[i % len(colors)]
        for i, row in unique_combinations.iterrows()
    }

    accession_to_combination = (
        classification_df
        .set_index('Accession no.')[['Serotype', 'Topotype']]
        .to_dict(orient='index')
    )

    def color_sequence(seq_name):
        accession = seq_name.split('_')[0]
        if accession in accession_to_combination:
            serotype = accession_to_combination[accession]['Serotype']
            topotype = accession_to_combination[accession]['Topotype']
            color = combination_to_color.get((serotype, topotype), 'black')
            return f"{seq_name}[&!color={color}]"
        return seq_name

    with open(input_nexus, 'r') as file:
        nexus_data = file.read()

    taxa_block_match = re.search(
        r'begin taxa;.*?taxlabels(.*?);',
        nexus_data,
        re.DOTALL | re.IGNORECASE
    )

    if not taxa_block_match:
        raise RuntimeError("TAXA / TAXLABELS block not found in NEXUS file")

    taxa_block = taxa_block_match.group(1)
    taxa_lines = taxa_block.strip().splitlines()

    modified_lines = []
    for i, line in enumerate(taxa_lines):
        line = line.strip()
        if line:
            colored_sequence = color_sequence(line)
            if i == 0:
                modified_lines.append(colored_sequence)
            else:
                modified_lines.append(f"\t{colored_sequence}")
        else:
            modified_lines.append(line)

    modified_taxa_block = "\n" + "\n".join(modified_lines)

    nexus_data = re.sub(
        r'(begin taxa;.*?taxlabels)(.*?)(;)',
        r'\1' + modified_taxa_block + r'\3',
        nexus_data,
        flags=re.DOTALL | re.IGNORECASE
    )

    with open(output_nexus, 'w') as file:
        file.write(nexus_data)

    print("Topotypes colored successfully!")
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
    description='Color taxa in a Nexus tree file according to topotype information provided in a classification table.')
    parser.add_argument('-t', '--table', required=True, help='Input Excel file containing accession-to-topotype mapping.')
    parser.add_argument('-i', '--input_tree',required=True,help='Input Nexus tree file.')
    parser.add_argument('-o', '--output_tree',required=True,help='Output Nexus tree file with colored taxa.')
    parser.add_argument("--max_colors",type=int,default=58,help="Maximum number of distinct colors (default: 58)")

    args = parser.parse_args()

    color_nexus(
        args.classification, 
        args.input_nexus, 
        args.output_nexus, 
        args.max_colors
    )
