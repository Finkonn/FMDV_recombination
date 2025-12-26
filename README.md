# Recombination dynamics of Foot-and-Mouth Disease virus
### Mate Malichava, Alexander Lukashev and Yulia Aleshina

The purpose of this work was to analyze the dynamics of natural recombination in Foot-and-Mouth Disease virus (FMDV, genus *Aphthovirus*). All data and scripts used in the study *"Recombination dynamics of Foot-and-Mouth Disease virus"* are freely available in this repository.

Bayesian and maximum likelihood phylogenetic trees, genome alignments, tables and mappings are located in `data/` directory.

`data/all_records.gb` - all publicly available complete genome sequences of FMDV (N = 1870) with exclusion of synthetic constructs downloaded from the GenBank database (accessed 05.03.2024)

Alignments of complete genome sequences are available in `data/alignments/`
- `All_serotypes_complete_genomes.fasta` - alignment of all complete genome sequences
- `Serotype_*_complete_genomes.fasta` - alignment of complete genome sequences specific to each serotype, where `*` is A, Asia1, C, O, SAT1, SAT2 or SAT3

Alignments of recombination-free regions are available in `data/alignments/recombinant_free_regions/`
- `Lpro.fasta`, `P1(VP2-VP3-VP1).fasta`, `P2.fasta`, `P3.fasta` - alignments of recombination-free regions, as determined by recombination analysis of all complete genome sequences using RDP4
- `A/`, `Asia1/`, `C/`, `O/`, `SAT1/`, `SAT2/`, `SAT3/` - folders containing alignments of recombination-free regions of each serotype, as determined by RDP4 analysis of complete genome sequences of each serotype (coordinates in file names, i.e. `A_1-511.fasta` represent ungapped coordinates in the complete genome alignment)

Maximum likelihood phylogenetic trees of recombination-free regions are available in `data/ml_trees/`
- `Lpro.treefile`, `P1(VP2-VP3-VP1).treefile`, `P2.treefile`, `P3.treefile` - phylogenetic trees of recombination-free regions, as determined by recombination analysis of all complete genome sequences using RDP4
- `A/`, `Asia1/`, `C/`, `O/`, `SAT1/`, `SAT2/`, `SAT3/` - folders containing phylogenetic trees of recombination-free regions of each serotype, as determined by RDP4 analysis of complete genome sequences of each serotype (coordinates in file names, i.e. `A_1-511.treefile` represent ungapped coordinates in the complete genome alignment)

Bayesian phylogenetic trees of recombination-free areas in P1 region of serotypes A, Asia1, SAT1 and 4 samples of serotype O are available at `data/bayes_trees/` in `.tree` format

Tables and mappings required for scripts and data preparation are available in `data/tables_and_maps/`
- `country_map.csv` - comma-separated mapping table with 2 columns: first column with names of countries to be replaced and second column with ISO-3166 codes of these countries to be replaced with. Used in the script `scripts/parse_genbank.ipynb` (not required) to replace country names with ISO-3166 codes for clean FASTA headers
- `host_map.csv` - mapping table with regular expressions of different host species names and simplified, standardized host names for replacement. Used in the script `scripts/parse_genbank.ipynb` (not required) to replace host names with simplified and standardized names for clean FASTA headers and overall interpretability. If editing the file, make sure to keep the same format, i.e. `\b[Ss]mall\s+ruminant(?:s)?\b,sheep` - make sure the replacement name is after a comma and without spaces
- `isolation_sources_to_remove.txt` - a text file with line-separated isolation sources. Used in the script `scripts/parse_genbank.ipynb` (not required) to filter out GenBank records with experimental sequences and vaccines
- `references_to_remove.txt` - a text file with line-separated references. Used in the script `scripts/parse_genbank.ipynb` (not required) to filter out GenBank records from specific experiments
- `metadata_with_gradient_tree_colors.csv` - table with gradient colors of each tip for phylogenetic tree visualisation using `scripts/plot_gradient_trees_with_bars.R` script

Python and R scripts used in the analysis are located in `scripts/` directory.

## Required python packages
- biopython>=1.78
- pandas==2.3.3
- matplotlib==3.10.8
- beautifulsoup4==4.12.2
- openpyxl==3.0.10

## Required R packages (for visualisation)
- randomcoloR
- colorspace
- RColorBrewer
- ggtree==3.16.0
- dplyr
- ggnewscale==0.5.1
- phytools==2.4-4
- ape==5.8-1
- ggplot2==3.5.2

## Parsing a GenBank file and generating a FASTA file with complete genome sequences

A jupyter notebook `parse_genbank.ipynb` is provided that filters and annotates GenBank records from the .gb file based on selected maps (country and host annotation; isolation sources and references filtration), and extracts coding sequences into a FASTA file with annotated headers (**>Accession/country/host/year_of_isolation/serotype**). It is not required, but recommended to use `country_map.csv` and `host_map.csv` (may be edited as needed) for clear, standardized annotation of sequences, and `references_to_remove.txt` and `isolation_sources_to_remove.txt` for removing unneeded records.
With default parameters, the output of the notebook is as follows:
- `CDS.fasta` and `metadata.csv` - raw FASTA file and corresponding metadata table without any mappings or filtration
- `CDS_mapped.fasta` and `metadata_mapped.csv` - FASTA file and corresponding metadata table only with host and country mappings
- `CDS_mapped_and_filtered` and `metadata_mapped_and_filtered.csv` - FASTA file and corresponding metadata table with host and country mappings filtered by `references_to_remove.txt` and `isolation_sources_to_remove.txt`

## Genotyping serotypes, topotypes and lineages - pipeline

### General usage steps (assumes root directory)

1. Launch the script `parse_wrlfmd.py` to create an updated, most recent classification table with representative strains for each FMDV topotype as available at https://www.wrlfmd.org/fmdv-genome/fmd-prototype-strains. Table contains serotypes, topotypes, lineages, sub-lineage, isolate names, accessions of genome sequences and corresponding references.
```bash
python scripts\genotyping\parse_wrlfmd.py -o classification_table.xlsx
```
2. Launch the script `fetch_sequences_from_genbank.py` to download available genome sequences of representative FMDV strains from GenBank using the classification table. This might take a few minutes depending on your download speed.
```bash
python scripts\genotyping\fetch_sequences_from_genbank.py -f classification_table.xlsx -o references.fasta -e your@email.com
```
3. Launch the script `rename_reference_headers.py` to rename the headers of reference sequences to a more readable format >Accession/Serotype/Topotype/Lineage
```bash
python scripts\genotyping\rename_reference_headers.py -f references.fasta -t classification_table.xlsx -o references_renamed.fasta
```
4. If you don't have a sequence alignment, use the script `trans_alignment.py` available at https://github.com/v-julia/alignment_of_orfs/blob/master/trans_alignment.py to make it. Following script will only work if your alignment has a reference sequence with "mat_peptide" feature (i.e. X00429). Make sure the reference sequence is the first sequence in the alignment and download a GenBank record file of that sequence. Now you can launch the script `get_peptide_coordinates_in_alignment.py` to get coordinates of mature peptides in your alignment.
```bash
python scripts\genotyping\get_peptide_coordinates_in_alignment.py -g your_reference_sequence.gb -f your_alignment.fasta -p peptide_coords.csv -o gapped_peptide_coords.csv
```
5. Launch the script `get_peptide_alignment.py` to make a VP1-only (by default) FASTA file of your alignment. You can also choose any other mature peptide using "-p" argument.
```bash
python scripts\genotyping\get_peptide_alignment.py -f your_alignment.fasta -t gapped_peptide_coords.csv
```
6. Add `references_renamed.fasta` to your VP1-only alignment. Align them using any alignment tool (MAFFT v7.520 recommended). Manually excise only VP1 region. For this, it is recommended to use AliView v1.3 (select the region and use the option "Copy (as FASTA)").
7. Construct a phylogenetic tree of that VP1-only alignment of your sequences and reference sequences, and save it in NEXUS format.
8. Launch the script `color_reference_branches.py` to highlight the branches of the reference sequences on your tree. This is not required, but makes next steps significantly easier.
```bash
python scripts\genotyping\color_reference_branches.py -i VP1_with_references_aligned.nexus -o VP1_with_reference_branches_colored.nexus
```
9. Manually color tips of serotypes, topotypes and/or lineages on the VP1 NEXUS tree:
    - For typing serotypes, manually color the clades corresponding to each serotype on the VP1 tree with different colors using, for example, FigTree v1.4.4. Save the colored tree in .nexus format and name it `serotype.nexus`.
    - For typing topotypes, manually color the clades corresponding to each topotype on the VP1 tree with different colors using, for example, FigTree v1.4.4 (this is easier with highlighted branches of reference sequences). Save the colored tree in .nexus format and name it `topotype.nexus`.
    - For typing lineages, manually color the clades corresponding to each lineage on the VP1 tree with different colors using, for example, FigTree v1.4.4 (this is easier with highlighted branches of reference sequences). Save the colored tree in .nexus format and name it `lineage.nexus`.
10. Launch the script `genotyping.py` to annotate your alignment FASTA file and create an annotated metadata table. The script has options to only annotate serotype, topotype or lineage - for this, use argument "-n" with a single colored .nexus file and choose a corresponding "--mode" (serotype, topotype or lineage). Or use argument "-d" with a directory of .nexus files named `serotype.nexus`, `topotype.nexus` and `lineage.nexus` to annotate all three. Examples:
```bash
python scripts\genotyping\genotyping.py -f your_alignment.fasta -o your_alignment_annotated.fasta -m metadata_annotated.xlsx -d directory_with_nexus_files
```
```bash
python scripts\genotyping\genotyping.py -f your_alignment.fasta -o your_alignment_annotated.fasta -m metadata_annotated.xlsx -n your_serotype_colored_tree.nexus --mode serotype
```

## Tree visualisation

Script `plot_gradient_trees_with_bars.R` can be used to visualise maximum likelihood trees. It requires a metadata table with columns:
- "GBAC" - FASTA headers - should be the same as tree tips
- "serotype" - serotype of the sequence 
- "Topotype" - topotype of the sequence 
- "lineage" - lineage of the sequence 
- "color" - gradient colors of tree tips
- "pool" - geographic pool of the sequence

Also requires a script `modified_gradients.R` to be present in the same directory and a directory with phylogenetic trees in .treefile format. Colors of topotypes, lineages and pools can be edited as needed. The script also transforms long topotype names to shorter ones in the legend, i.e. "MIDDLE-EAST-SOUTH-ASIA-ME-SA" ~ "ME-SA". The output is trees in .svg format with gradient colored tip labels (according to serotypes) and color bars of serotypes (S), topotypes (T), lineages (L) and pools (P).

## Description and usage of scripts

### parse_wrlfmd.py

Parse prototype strain information from the WRLFMD website and save it as a classification table in .xlsx format

#### Usage
```
-u, --url: URL to be parsed. Default = "https://www.wrlfmd.org/fmdv-genome/fmd-prototype-strains"
-o, --output: output file in .xlsx format. (Required)
```
```bash
python parse_wrlfmd.py -o <output.xlsx>
```

### fetch_sequences_from_genbank.py

Download sequences in the classification table from GenBank using Entrez

#### Usage
```
-f, --file: path to the classification table. (Required)
-e, --email: your email address. (Required)
-o, --output: output FASTA. (Required)
```
```bash
python fetch_sequences_from_genbank.py -f <input_table> -o <output_fasta> -e <email>
```

### rename_reference_headers.py 

Rename the headers of reference sequences to a more readable format >Accession/Serotype/Topotype/Lineage

#### Usage
```
-f, --fasta: path to the input FASTA file. (Required)
-t, --excel: path to the classification table. (Required)
-o, --output: path to the output FASTA file. (Required)
```
```bash
python rename_reference_headers.py -f <input_fasta> -t <input_table> -o <output_fasta>
```

### get_peptide_coordinates_in_alignment.py

This script extracts mature peptide (mat_peptide) coordinates from a GenBank file, maps these coordinates to the first sequence in a FASTA file and accounts for gaps. It outputs two CSV files: one with the original peptide coordinates and another with the gapped coordinates.

#### Usage  
```
-g, --genbank: path to the input GenBank file containing mat_peptide features. (Required)
-f, --fasta: path to the input FASTA file containing the reference sequence. (Required)
-p, --peptide_output: path to save the output CSV file with peptide coordinates. (Required)
-o, --gapped_output: path to save the output CSV file with gapped peptide coordinates. (Required)
```
```bash
python get_peptide_coordinates_in_alignment.py -g <genbank_file> -f <fasta_file> -p <peptide_output_file> -o <gapped_output_file>
```

### get_peptide_alignment.py

This script creates a VP1-only alignment based on gapped coordinates of the complete genome alignment.

#### Usage
```
-f, --input_fasta: path to the input FASTA alignment. (Required)
-t, --input_table: path to the CSV file with gapped coordinates. (Required)
-p, --product: peptide to make alignment of. Default = VP1
```
```bash
python get_peptide_alignment.py -f <input_fasta> -t <input_table>
```

### color_reference_branches.py

Color the branches of the reference sequences in the NEXUS tree red.

#### Usage
```
-i, --input: input NEXUS tree. (Required)
-o, --output: output NEXUS tree. (Required)
```
```bash
python color_reference_branches.py -i <input_nexus_tree> -o <output_nexus_tree>
```

### genotyping.py

Annotate FASTA headers (>Accession/.../Serotype/Topotype/Lineage) and create an annotated metadata table

#### Usage 
```
-f, --fasta: input FASTA file. (Required)
-o, --output: output FASTA file. (Required)
-m, --metadata: output metadata file.
-n, --nexus: single colored NEXUS tree.
-d, --nexus_dir: directory with serotype.nexus, topotype.nexus and lineage.nexus
--mode: typing mode - serotype, topotype or lineage. Required for single colored NEXUS tree
```
```bash
python genotyping.py -f <input_fasta> -o <output_fasta> -m <output_table> -d <directory_with_trees>
```
```bash
python genotyping.py -f <input_fasta> -o <output_fasta> -m <output_table> -n <nexus_tree> --mode <serotype/topotype/lineage>
```