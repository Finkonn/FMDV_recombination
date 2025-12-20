# Recombination dynamics of Foot-and-Mouth Disease virus
### Mate Malichava, Alexander Lukashev and Yulia Aleshina

The purpose of this work was to analyze the dynamics of natural recombination in Foot-and-Mouth Disease virus (FMDV, genus *Aphthovirus*). All data and scripts used in the study *"Recombination dynamics of Foot-and-Mouth Disease virus"* are freely available in this repository.

## Required python packages
- biopython>=1.78
- pandas==2.3.3
- matplotlib==3.10.8
- beautifulsoup4==4.12.2

## Required R packages (for visualisation)
- randomcoloR
- colorspace
- RColorBrewer
- ggtree
- dplyr
- ggnewscale
- phytools
- ape
- ggplot2

## Parsing a GenBank file and generating a FASTA file with complete genome sequences

A jupyter notebook (*parse_genbank.ipynb*) is provided that filters and annotates GenBank records from the .gb file based on selected maps (country and host annotation; isolation sources and references filtration), and extracts coding sequences into a FASTA file with annotated headers (**>Accession/country/host/year_of_isolation/serotype**).

## Genotyping serotypes, topotypes and lineages

First, a classification table is required for genotyping topotypes and lineages. The script *parse_wrlfmd.py* parses the WRLFMD prototype strains page (https://www.wrlfmd.org/fmdv-genome/fmd-prototype-strains) and creates a classification table with reference sequences (serotypes, topotypes and some lineages annotated). Sequences in this table can be downloaded using *fetch_sequences_from_genbank.py*. These sequences can then be added to a VP1 alignment to construct a NEXUS tree.

Serotype-colored NEXUS tree can be manually colored based on clades correponding to each serotype using FigTree v1.4.4. A script for obtaining VP1 coordinates in the alignment based on a reference sequence is available at *get_VP1_coordinates_in_alignment.py*. A script for making a VP1-only alignment is available at *get_VP1_alignment.py* (uses the coordinates table obtained from the previous script). 

Topotype-colored NEXUS tree can be manually colored based on reference topotypes using FigTree v1.4.4. Reference topotype sequences can be automatically colored using a script *color_tree_topotypes.py* (classification table is required).

Lineage-colored NEXUS tree can be manually colored based on reference lineages using FigTree v1.4.4. Reference lineage sequences can be automatically colored using a script *color_tree_lineages.py* (classification table is required).

Now, script *genotyping.py* can be used to annotate sequences in a FASTA file based on a table with colors and corresponding serotypes/topotypes or lineages. For the script to work, a colored tree in NEXUS format is required. 
This is a modified version, original available at https://github.com/v-julia/GenAlignment.

### parse_wrlfmd.py

Parse prototype strain information from the WRLFMD website and save it as a classification table in .xlsx format.

#### Usage
```
-u, --url: URL to be parsed. Default = "https://www.wrlfmd.org/fmdv-genome/fmd-prototype-strains"
-o, --output: output file in .xlsx format. (Required)
```
```bash
python parse_wrlfmd.py -o <output.xlsx>
```

### fetch_sequences_from_genbank.py

Download sequences in the classification table from GenBank using Entrez.

#### Usage
```
-f, --file: path to the classification table. (Required)
-e, --email: your email address. (Required)
-o, --output: output FASTA. (Required)
```
```bash
python fetch_sequences_from_genbank.py -f <input_table> -o <output_fasta>
```

### get_VP1_coordinates_in_alignment.py

This script extracts mature peptide (mat_peptide) coordinates from a GenBank file, maps these coordinates to a reference sequence from a FASTA file, and accounts for gaps in the reference sequence. It outputs two CSV files: one with the original protein coordinates and another with the gapped coordinates.

#### Usage  
```
-g, --genbank: path to the input GenBank file containing mat_peptide features. (Required)
-f, --fasta: path to the input FASTA file containing the reference sequence. (Required)
-p, --protein_output: path to save the output CSV file with protein coordinates. (Required)
-o, --gapped_output: path to save the output CSV file with gapped protein coordinates. (Required)
```
```bash
python get_VP1_coordinates_in_alignment.py -g <genbank_file> -f <fasta_file> -p <protein_output_file> -o <gapped_output_file>
```

### get_VP1_alignment.py

This script creates a VP1-only alignment based on gapped coordinates of the complete genome alignment.

#### Usage
```
-i, --input_fasta: path to the input FASTA alignment. (Required)
-t, --input_table: path to the CSV file with gapped coordinates. (Required)
-p, --product: protein to make alignment of. Default = VP1
```
```bash
python get_VP1_alignment.py -i <input_fasta> -t <input_table>
```

### color_tree_topotypes.py

Color topotype reference sequences in a NEXUS tree file according to classification table.

#### Usage
```
-t, --table: classification table. (Required)
-i, --input_tree: input NEXUS tree file. (Required)
-o, --output_tree: colored output NEXUS tree file. (Required)
```
```bash
python color_tree_topotypes.py -t <input_table> -i <input_tree> -o <output_tree>
```

### color_tree_lineages.py

Color lineage reference sequences in a NEXUS tree file according to classification table.

#### Usage
```
-t, --table: classification table. (Required)
-i, --input_tree: input NEXUS tree file. (Required)
-o, --output_tree: colored output NEXUS tree file. (Required)
```
```bash
python color_tree_lineages.py -t <input_table> -i <input_tree> -o <output_tree>
```

### genotyping.py

Map colors from NEXUS tree to serotypes, topotypes or lineages. Change 'Unknown' serotypes to new ones and check for conflicts. 

#### Usage 
```
-in_rep, --input_rep_fasta: input directory with files in FASTA format (Required)
-in_tree, --input_file_tree: input colored NEXUS tree. (Required)
-in_csv, --input_file_csv: input table in csv format with colors (in HEX) and genotypes. (Required)
```
```bash
python genotyping.py -in_rep <input_directory> -in_tree <input_tree> -in_csv <input_table>
```