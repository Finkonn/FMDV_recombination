import pandas as pd
import re
import argparse
import matplotlib.colors as mcolors

def generate_colors(n):
    return [
        mcolors.to_hex(c)
        for c in mcolors.hsv_to_rgb([(i / n, 1.0, 1.0) for i in range(n)])
    ]

def color_tree_by_lineage(classification_table, input_tree, output_tree):
    classification_df = pd.read_excel(classification_table)
    classification_df = classification_df[classification_df['Accession no.'].notna()]

    unique_lineages = classification_df['Lineage'].dropna().unique()
    colors = generate_colors(len(unique_lineages))

    lineage_to_color = {
        lineage: colors[i % len(colors)]
        for i, lineage in enumerate(unique_lineages)
    }

    accession_to_lineage = (
        classification_df
        .set_index('Accession no.')['Lineage']
        .to_dict()
    )

    def color_sequence(seq_name):
        accession = seq_name.split('_')[0]
        if accession in accession_to_lineage:
            lineage = accession_to_lineage[accession]
            color = lineage_to_color.get(lineage, 'black')
            return f"{seq_name}[&!color={color}]"
        return seq_name

    with open(input_tree, 'r') as file:
        nexus_data = file.read()

    taxa_block_match = re.search(
        r'begin taxa;.*?taxlabels(.*?);',
        nexus_data,
        re.DOTALL
    )

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
        rf'\1{modified_taxa_block}\3',
        nexus_data,
        flags=re.DOTALL
    )

    with open(output_tree, 'w') as file:
        file.write(nexus_data)

    print("Lineages colored successfully!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
    description='Color taxa in a Nexus tree file according to lineage information provided in a classification table.')
    parser.add_argument('-t', '--table', required=True, help='Input Excel file containing accession-to-lineage mapping.')
    parser.add_argument('-i', '--input_tree',required=True,help='Input Nexus tree file.')
    parser.add_argument('-o', '--output_tree',required=True,help='Output Nexus tree file with colored taxa.')

    args = parser.parse_args()

    color_tree_by_lineage(
        args.table,
        args.input_tree,
        args.output_tree
    )

