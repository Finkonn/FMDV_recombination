# Recombination dynamics of Foot-and-Mouth Disease virus
### Mate Malichava, Alexander Lukashev and Yulia Aleshina

The purpose of this work was to analyze the dynamics of natural recombination in Foot-and-Mouth Disease virus (FMDV, genus *Aphthovirus*). All data and scripts used in the study *"Recombination dynamics of Foot-and-Mouth Disease virus"* are freely available in this repository.

Bayesian and maximum likelihood phylogenetic trees, genome alignments, tables and mappings are located in `data/` directory.

- `data/all_records.gb` - all publicly available complete genome sequences of FMDV (N = 1870) with exclusion of synthetic constructs downloaded from the GenBank database (accessed 05.03.2024)

- Alignments of complete genome sequences are available in `data/alignments/`
    - `All_serotypes_complete_genomes.fasta` - alignment of all complete genome sequences
    - `Serotype_*_complete_genomes.fasta` - alignment of complete genome sequences specific to each serotype, where `*` is A, Asia1, C, O, SAT1, SAT2 or SAT3

- Alignments of recombination-free regions are available in `data/alignments/recombinant_free_regions/`
    - `Lpro.fasta`, `P1(VP2-VP3-VP1).fasta`, `P2.fasta`, `P3.fasta` - alignments of recombination-free regions, as determined by recombination analysis of all complete genome sequences using RDP4
    - `A/`, `Asia1/`, `C/`, `O/`, `SAT1/`, `SAT2/`, `SAT3/` - folders containing alignments of recombination-free regions of each serotype, as determined by RDP4 analysis of complete genome sequences of each serotype (coordinates in file names, i.e. `A_1-511.fasta` represent ungapped coordinates in the complete genome alignment)

- Maximum likelihood phylogenetic trees of recombination-free regions are available in `data/ml_trees/recombinant_free_regions/`
    - `Lpro.treefile`, `P1(VP2-VP3-VP1).treefile`, `P2.treefile`, `P3.treefile` - phylogenetic trees of recombination-free regions, as determined by recombination analysis of all complete genome sequences using RDP4
    - `A/`, `Asia1/`, `C/`, `O/`, `SAT1/`, `SAT2/`, `SAT3/` - folders containing phylogenetic trees of recombination-free regions of each serotype, as determined by RDP4 analysis of complete genome sequences of each serotype (coordinates in file names, i.e. `A_1-511.treefile` represent ungapped coordinates in the complete genome alignment)

- Examples of maximum likelihood VP1 NEXUS trees manually colored by serotypes and semi-automatically colored by topotypes and lineages are available in `data/ml_trees/`

- Bayesian phylogenetic trees of recombination-free areas in P1 region of serotypes A, Asia1, SAT1 and 4 samples of serotype O are available at `data/bayes_trees/` in `.tree` format

- Tables and mappings required for scripts and data preparation are available in `data/tables_and_maps/`
    - `country_map.csv` - comma-separated mapping table with 2 columns: first column with names of countries to be replaced and second column with ISO-3166 codes of these countries to be replaced with. Used in the script `scripts/parse_genbank.ipynb` (not required) to replace country names with ISO-3166 codes for clean FASTA headers
    - `host_map.csv` - mapping table with regular expressions of different host species names and simplified, standardized host names for replacement. Used in the script `scripts/parse_genbank.ipynb` (not required) to replace host names with simplified and standardized names for clean FASTA headers and overall interpretability. If editing the file, make sure to keep the same format, i.e. `\b[Ss]mall\s+ruminant(?:s)?\b,sheep` - make sure the replacement name is after a comma and without spaces
    - `isolation_sources_to_remove.txt` - a text file with line-separated isolation sources. Used in the script `scripts/parse_genbank.ipynb` (not required) to filter out GenBank records with experimental sequences and vaccines
    - `references_to_remove.txt` - a text file with line-separated references. Used in the script `scripts/parse_genbank.ipynb` (not required) to filter out GenBank records from specific experiments
    - `metadata_with_gradient_tree_colors.csv` - 

Python and R scripts used in the analysis are located in `scripts/` directory.

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

A jupyter notebook `parse_genbank.ipynb` is provided that filters and annotates GenBank records from the .gb file based on selected maps (country and host annotation; isolation sources and references filtration), and extracts coding sequences into a FASTA file with annotated headers (**>Accession/country/host/year_of_isolation/serotype**). It is not required, but recommended to use `country_map.csv` and `host_map.csv` (may be edited as needed) for clear, standardized annotation of sequences, and `references_to_remove.txt` and `isolation_sources_to_remove.txt` for removing unneeded records.
With default parameters, the output of the notebook is as follows:
- `CDS.fasta` and `metadata.csv` - raw FASTA file and corresponding metadata table without any mappings or filtration
- `CDS_mapped.fasta` and `metadata_mapped.csv` - FASTA file and corresponding metadata table only with host and country mappings
- `CDS_mapped_and_filtered` and `metadata_mapped_and_filtered.csv` - FASTA file and corresponding metadata table with host and country mappings filtered by `references_to_remove.txt` and `isolation_sources_to_remove.txt`

## Genotyping serotypes, topotypes and lineages - pipeline

First, a classification table is required for genotyping topotypes and lineages. The script `parse_wrlfmd.py` parses the WRLFMD prototype strains page (https://www.wrlfmd.org/fmdv-genome/fmd-prototype-strains) and creates a classification table with reference sequences (serotypes, topotypes and some lineages annotated). Sequences in this table can be downloaded using `fetch_sequences_from_genbank.py`. These sequences can then be added to a VP1 alignment to construct a NEXUS tree.

Serotype-colored NEXUS tree can be manually colored based on clades correponding to each serotype using FigTree v1.4.4. A script for obtaining VP1 coordinates in the alignment based on a reference sequence is available at `get_VP1_coordinates_in_alignment.py`. A script for making a VP1-only alignment is available at `get_VP1_alignment.py` (uses the coordinates table obtained from the previous script). 

Topotype-colored NEXUS tree can be manually colored based on reference topotypes using FigTree v1.4.4. Reference topotype sequences can be automatically colored using a script `color_tree_topotypes.py` (classification table is required).

Lineage-colored NEXUS tree can be manually colored based on reference lineages using FigTree v1.4.4. Reference lineage sequences can be automatically colored using a script `color_tree_lineages.py` (classification table is required).

Now, script `genotyping.py` can be used to annotate sequences in a FASTA file based on a table with colors and corresponding serotypes/topotypes or lineages. For the script to work, a colored tree in NEXUS format is required. 
This is a modified version, original available at https://github.com/v-julia/GenAlignment.

### General usage steps

1. Launch the script `parse_wrlfmd.py` to create an updated, most recent classification table with representative strains for each FMDV topotype as available at https://www.wrlfmd.org/fmdv-genome/fmd-prototype-strains. Table contains serotypes, topotypes, lineages, sub-lineage, isolate names, accessions of genome sequences and corresponding references.
```bash
python parse_wrlfmd.py -o classification_table.xlsx
```
2. Launch the script `fetch_sequences_from_genbank.py` to download available genome sequences of representative FMDV strains from GenBank using the classification table. This might take a few minutes depending on your download speed.
```bash
python fetch_sequences_from_genbank.py -f classification_table.xlsx -o references.fasta -e your@email.com
```
3. If you don't have a sequence alignment, use the script `trans_alignment.py` available at https://github.com/v-julia/alignment_of_orfs/blob/master/trans_alignment.py to make it. Following script will only work if your alignment has a reference sequence with "mat_peptide" feature (i.e. X00429). Make sure the reference sequence is the first sequence in the alignment and download a GenBank record file with that sequence. Now you can launch the script `get_VP1_coordinates_in_alignment.py` to get coordinates of mature VP1 protein in your alignment.
```bash
python 
```

4. 
## Description and usage of scripts

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