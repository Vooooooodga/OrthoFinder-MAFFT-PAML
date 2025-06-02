#!/usr/bin/env python3

import argparse
import os
import pandas as pd
from ete3 import Tree

def _mark_tree_and_save(tree_path, sociality_map, target_sociality, paml_marker, output_folder):
    """
    Loads a tree, marks nodes based on descendant sociality, and saves the modified tree.
    """
    try:
        # format=0: ETE attempts to parse topology, branch lengths, support values, and names. Robust for reading.
        tree = Tree(tree_path, format=0)
    except Exception as e:
        print(f"Warning: Could not parse tree {os.path.basename(tree_path)}: {e}")
        return

    # Step 1: Annotate leaves with their sociality status based on the sociality_map
    for leaf in tree.get_leaves():
        original_name = leaf.name  # Assume leaf names in the input tree are original
        social_level = sociality_map.get(original_name)
        
        if social_level is None:
            print(f"Warning: Species '{original_name}' in tree '{os.path.basename(tree_path)}' not found in sociality file. Treating as non-target.")
            leaf.add_features(is_target_sociality=False, sociality_known=False)
        elif social_level == target_sociality:
            leaf.add_features(is_target_sociality=True, sociality_known=True)
        else:
            leaf.add_features(is_target_sociality=False, sociality_known=True)

    # Step 2: Traverse and identify nodes whose ALL leaf descendants are of target_sociality
    nodes_to_be_marked = []
    for node in tree.traverse("postorder"):  # Postorder is fine, any order works for this logic
        all_desc_are_target = True
        
        current_node_leaf_descendants = node.get_leaves()

        if not current_node_leaf_descendants: # Node is a leaf itself, or an internal node with no leaf descendants
            if node.is_leaf(): # If it's a leaf, check its own features
                if not (hasattr(node, 'is_target_sociality') and node.is_target_sociality):
                    all_desc_are_target = False
            else: # Internal node with no leaves (e.g., from pruning or unusual tree structure)
                all_desc_are_target = False 
        else: # Node is an internal node with leaf descendants
            for leaf_descendant in current_node_leaf_descendants:
                if not (hasattr(leaf_descendant, 'is_target_sociality') and leaf_descendant.is_target_sociality):
                    all_desc_are_target = False
                    break
        
        # Only mark if all descendants are target AND the node actually has descendants (or is a target leaf)
        if all_desc_are_target and (node.is_leaf() or current_node_leaf_descendants):
            nodes_to_be_marked.append(node)

    # Step 3: Mark the identified nodes by appending the PAML marker to their names
    for node_to_mark in nodes_to_be_marked:
        current_name = node_to_mark.name if node_to_mark.name else ""
        # Avoid double-marking if by some chance the script runs over already marked names,
        # or if a node name naturally ends with the marker string (unlikely for typical markers).
        if not current_name.endswith(paml_marker):
            node_to_mark.name = current_name + paml_marker
        # If node_to_mark.name was None or empty (common for internal nodes), it becomes just the paml_marker.
        # ETE handles this correctly when writing the Newick string, e.g., (...)#1

    # Step 4: Save the modified tree
    base_name = os.path.basename(tree_path)
    name_part, ext_part = os.path.splitext(base_name)
    
    # Ensure the output extension is consistent or preserved
    if not ext_part and base_name.endswith(".treefile"): # e.g. original "name.treefile"
         ext_part = ".treefile"
    elif not ext_part: # e.g. original "name"
         ext_part = ".treefile" # Default to .treefile if no extension

    output_filename = f"{name_part}_marked{ext_part}"
    output_path = os.path.join(output_folder, output_filename)

    try:
        # format=1: Writes Newick with internal node names (if they exist) and branch lengths.
        # This format is generally compatible with PAML.
        tree.write(outfile=output_path, format=1)
        # print(f"Saved marked tree to {output_path}")
    except Exception as e:
        print(f"Error writing marked tree {output_path}: {e}")
        
    # Clean up temporary features from nodes (optional, good practice)
    for node in tree.traverse():
        if hasattr(node, 'is_target_sociality'):
            del node.is_target_sociality
        if hasattr(node, 'sociality_known'):
            del node.sociality_known

def main():
    parser = argparse.ArgumentParser(
        description="Batch add PAML foreground markers to gene tree files based on species sociality.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--tree_folder", required=True, help="Folder containing gene tree files (e.g., with .treefile extension).")
    parser.add_argument("--sociality_file", required=True, help="Path to a CSV or TSV file containing species sociality data.\nMust have 'Organism Name' and 'Social level' columns.")
    parser.add_argument("--output_folder", required=True, help="Folder where marked tree files will be saved.")
    parser.add_argument("--target_sociality", required=True, help="The sociality level to be marked as foreground (e.g., 'Advanced Eusocial').")
    parser.add_argument("--paml_marker", required=True, help="The PAML marker to append to foreground branches (e.g., '#1').")

    args = parser.parse_args()

    # Validate input paths
    if not os.path.isdir(args.tree_folder):
        print(f"Error: Tree folder '{args.tree_folder}' not found.")
        return
    if not os.path.isfile(args.sociality_file):
        print(f"Error: Sociality file '{args.sociality_file}' not found.")
        return

    # Create output folder if it doesn't exist
    if not os.path.isdir(args.output_folder):
        try:
            os.makedirs(args.output_folder, exist_ok=True)
            print(f"Created output folder: {args.output_folder}")
        except OSError as e:
            print(f"Error: Could not create output folder '{args.output_folder}': {e}")
            return

    # Load sociality data
    sociality_map = {}
    try:
        file_ext = os.path.splitext(args.sociality_file)[1].lower()
        if file_ext == ".csv":
            social_df = pd.read_csv(args.sociality_file)
        elif file_ext == ".tsv" or file_ext == ".txt": # Allow .txt for TSV too
            social_df = pd.read_csv(args.sociality_file, sep='\t')
        else:
            print(f"Warning: Unknown sociality file extension '{file_ext}'. Attempting to parse as CSV, then TSV.")
            try:
                social_df = pd.read_csv(args.sociality_file)
            except pd.errors.ParserError:
                try:
                    social_df = pd.read_csv(args.sociality_file, sep='\t')
                except Exception as e_tsv:
                    print(f"Error: Could not parse sociality file '{args.sociality_file}' as CSV or TSV. Error: {e_tsv}")
                    return
            except Exception as e_csv:
                 print(f"Error: Could not parse sociality file '{args.sociality_file}' as CSV. Error: {e_csv}")
                 return


        if "Organism Name" not in social_df.columns or "Social level" not in social_df.columns:
            print(f"Error: Sociality file '{args.sociality_file}' must contain 'Organism Name' and 'Social level' columns.")
            return
        
        for _, row in social_df.iterrows():
            sociality_map[str(row["Organism Name"])] = str(row["Social level"])
        
        if not sociality_map:
            print(f"Warning: No sociality data loaded from '{args.sociality_file}'. This may result in no branches being marked.")
        else:
            print(f"Loaded {len(sociality_map)} species sociality entries from '{args.sociality_file}'.")

    except Exception as e:
        print(f"Error loading or parsing sociality file '{args.sociality_file}': {e}")
        return

    # Process each tree file
    processed_count = 0
    tree_files_found = 0
    for filename in os.listdir(args.tree_folder):
        if filename.lower().endswith(".treefile"): # Case-insensitive check for .treefile
            tree_files_found +=1
            tree_file_path = os.path.join(args.tree_folder, filename)
            print(f"Processing {filename}...")
            _mark_tree_and_save(tree_file_path, sociality_map, args.target_sociality, args.paml_marker, args.output_folder)
            processed_count += 1
    
    if tree_files_found == 0:
        print(f"No '.treefile' files found in '{args.tree_folder}'.")
    else:
        print(f"\nSuccessfully processed {processed_count} out of {tree_files_found} '.treefile' files found.")
        print(f"Marked trees have been saved to: {args.output_folder}")

if __name__ == "__main__":
    main()