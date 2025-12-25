import argparse

REFERENCE_COLOR = "[&!color=#ff0000]"

def is_reference(label: str) -> bool:
    return label is not None and len(label.split("/")) == 4

def color_references_in_newick(newick_str: str) -> str:
    import re

    # This regex matches all tip labels in single quotes (e.g., 'AB079061/JPN/Unknown/2000/O')
    pattern = r"'([^']+)'"
    
    def replacer(match):
        label = match.group(1)
        if is_reference(label):
            return f"'{label}'{REFERENCE_COLOR}"
        else:
            return f"'{label}'"

    return re.sub(pattern, replacer, newick_str)

def main():
    parser = argparse.ArgumentParser(description="Color reference sequence branches in a NEXUS tree (text-based)")
    parser.add_argument("-i", "--input", required=True, help="Input NEXUS tree")
    parser.add_argument("-o", "--output", required=True, help="Output colored NEXUS tree")
    args = parser.parse_args()

    # Read the NEXUS tree as text
    with open(args.input, "r") as f:
        content = f.read()

    # Find the tree block (assumes single tree)
    import re
    match = re.search(r"(tree\s+\w+\s*=\s*)([^;]+);", content, flags=re.IGNORECASE)
    if not match:
        raise ValueError("No tree found in NEXUS file")
    
    tree_prefix = match.group(1)
    tree_newick = match.group(2)

    # Color the references
    colored_newick = color_references_in_newick(tree_newick)

    # Replace the old tree with colored tree
    content_colored = content[:match.start(2)] + colored_newick + content[match.end(2):]

    # Write back to output
    with open(args.output, "w") as f:
        f.write(content_colored)

if __name__ == "__main__":
    main()
