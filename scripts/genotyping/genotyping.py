import argparse
import os
import re
import csv
from Bio import SeqIO
from openpyxl import Workbook
from openpyxl.styles import PatternFill

def parse_colored_nexus(nexus_file):
    """Parse a colored Nexus tree and return {color_hex: [sequence_ids]}"""
    color_map = {}
    in_taxlabels = False

    with open(nexus_file) as f:
        for line in f:
            if line.strip().startswith("taxlabels"):
                in_taxlabels = True
                continue
            if in_taxlabels:
                if line.strip() == ";":
                    break
                color_match = re.search(r"(#[0-9a-fA-F]{6})", line)
                name_match = re.search(r"([A-Za-z0-9_\-./]+)", line)
                if name_match:
                    name = name_match.group(1)
                    color = color_match.group(1) if color_match else None
                    if color:
                        color_map.setdefault(color, []).append(name)
    return color_map

def extract_reference_annotations(color_map):
    """Identify reference sequences and extract annotations: returns {color: (serotype, topotype, lineage)}"""
    refs = {}
    for color, names in color_map.items():
        for name in names:
            parts = name.split("/")
            if len(parts) == 4:
                _, serotype, topotype, lineage = parts
                refs[color] = (serotype, topotype, lineage)
                break
    return refs

def apply_typing(records, color_map, ref_ann, mode):
    """
    Apply serotype / topotype / lineage typing to FASTA records.
    Missing or uncolored tips are assigned 'Unknown'.
    """
    name_to_color = {
        name.split("/")[0]: color
        for color, names in color_map.items()
        for name in names
    }

    for record in records:
        acc = record.id.split("/")[0]
        parts = record.id.split("/")

        # skip reference sequences
        if len(parts) == 4:
            continue

        # default unknowns
        serotype, topotype, lineage = "Unknown", "Unknown", "Unknown"

        # fill from colored reference if available
        color = name_to_color.get(acc)
        if color and color in ref_ann:
            s, t, l = ref_ann[color]
            serotype = s if s and s != "-" else "Unknown"
            topotype = t if t and t != "-" else "Unknown"
            lineage = l if l and l != "-" else "Unknown"

        # update record id depending on mode
        if mode == "serotype":
            parts[-1] = serotype
        elif mode == "topotype":
            parts.append(topotype)
        elif mode == "lineage":
            parts.append(lineage)

        record.id = "/".join(parts)
        record.description = ""

    return records

def save_metadata(records, output_xlsx):
    wb = Workbook()
    ws = wb.active
    ws.title = "Metadata"

    headers = ["Accession", "Country", "Host", "Isolation_year", "Serotype", "Topotype", "Lineage"]
    ws.append(headers)

    for record in records:
        parts = record.id.split("/")
        accession = parts[0] if len(parts) > 0 and parts[0] != "Unknown" else ""
        country = parts[1] if len(parts) > 1 and parts[1] != "Unknown" else ""
        host = parts[2] if len(parts) > 2 and parts[2] != "Unknown" else ""
        year = parts[3] if len(parts) > 3 and parts[3] != "Unknown" else ""
        serotype = parts[4] if len(parts) > 4 and parts[4] != "Unknown" else ""
        topotype = parts[5] if len(parts) > 5 and parts[5] != "Unknown" else ""
        lineage = parts[6] if len(parts) > 6 and parts[6] != "Unknown" else ""

        ws.append([accession, country, host, year, serotype, topotype, lineage])

    wb.save(output_xlsx)

def run_single_mode(fasta_in, fasta_out, nexus_file, mode, metadata_out=None):
    records = list(SeqIO.parse(fasta_in, "fasta"))
    color_map = parse_colored_nexus(nexus_file)
    ref_ann = extract_reference_annotations(color_map)
    records = apply_typing(records, color_map, ref_ann, mode)
    SeqIO.write(records, fasta_out, "fasta")
    if metadata_out:
        save_metadata(records, metadata_out)

def run_batch_mode(fasta_in, fasta_out, nexus_dir, metadata_out=None):
    temp_fasta = fasta_in
    temp_files = []

    for mode in ["serotype", "topotype", "lineage"]:
        nexus_file = os.path.join(nexus_dir, f"{mode}.nexus")
        temp_out = fasta_out + f".{mode}.tmp"
        temp_files.append(temp_out)
        run_single_mode(temp_fasta, temp_out, nexus_file, mode)
        temp_fasta = temp_out

    os.rename(temp_fasta, fasta_out)

    # remove other temp files
    for f in temp_files[:-1]:
        if os.path.exists(f):
            os.remove(f)

    # save metadata
    if metadata_out:
        records = list(SeqIO.parse(fasta_out, "fasta"))
        save_metadata(records, metadata_out)

def main():
    parser = argparse.ArgumentParser(description="Reference-based genotyping from colored Nexus trees")
    parser.add_argument("-f", "--fasta", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file")
    parser.add_argument("-m", "--metadata", help="Output metadata file")
    parser.add_argument("-n", "--nexus", help="Single colored Nexus tree")
    parser.add_argument("-d", "--nexus_dir", help="Directory with serotype.nexus, topotype.nexus, lineage.nexus")
    parser.add_argument("--mode", choices=["serotype", "topotype", "lineage"], help="Typing mode (required for single nexus)")

    args = parser.parse_args()

    if args.nexus_dir:
        run_batch_mode(args.fasta, args.output, args.nexus_dir, metadata_out=args.metadata)
    else:
        run_single_mode(args.fasta, args.output, args.nexus, args.mode, metadata_out=args.metadata)

if __name__ == "__main__":
    main()
