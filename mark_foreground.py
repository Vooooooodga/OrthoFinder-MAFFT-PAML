#!/usr/bin/env python3

import argparse
import os
import pandas as pd
from Bio import Phylo # 使用 BioPython

def _mark_tree_and_save_biopython(tree_path, sociality_map, target_sociality, paml_marker, output_folder):
    """
    Loads a tree using BioPython, marks nodes based on descendant sociality, 
    and saves the modified tree.
    """
    try:
        tree = Phylo.read(tree_path, "newick")
    except Exception as e:
        print(f"Warning: Could not parse tree {os.path.basename(tree_path)} with BioPython: {e}")
        return

    # Step 1: Create a mapping from leaf names to their 'is_target_sociality' status
    # This map will store whether a leaf (by its name) belongs to the target sociality group.
    leaf_is_target_status = {}
    for leaf_node in tree.get_terminals(): # get_terminals() gets all leaf nodes
        original_name = leaf_node.name
        if original_name is None: # Should not happen for leaves in valid Newick, but good to check
            print(f"Warning: Leaf node without a name found in tree '{os.path.basename(tree_path)}'. Skipping.")
            continue

        social_level = sociality_map.get(original_name)
        
        if social_level is None:
            print(f"Warning: Species '{original_name}' in tree '{os.path.basename(tree_path)}' not found in sociality file. Treating as non-target.")
            leaf_is_target_status[original_name] = False
        elif social_level == target_sociality:
            leaf_is_target_status[original_name] = True
        else:
            leaf_is_target_status[original_name] = False

    # Step 2: Traverse all clades (nodes) and identify those whose ALL leaf descendants are of target_sociality.
    # A clade (node) is marked if all its descendant leaves are of the target sociality.
    nodes_to_be_marked = []
    # We can iterate through all clades. Order doesn't strictly matter for this check,
    # but postorder (leaves first) is a common way to think about tree properties.
    for clade in tree.find_clades(order="postorder"): 
        all_desc_leaf_are_target = True
        
        # Get all terminal descendants (leaves) of the current clade
        descendant_leaves = list(clade.get_terminals())

        if not descendant_leaves:
            # This case means the clade itself is a leaf (and thus descendant_leaves should contain itself),
            # or it's an internal node with no descendant leaves (e.g. a pruned branch, or a malformed tree part).
            if clade.is_terminal(): # If it's a leaf node
                # Check its own status from the map. Default to False if name somehow missing from map (should not happen).
                if not leaf_is_target_status.get(clade.name, False):
                    all_desc_leaf_are_target = False
            else: # Internal node with no descendant leaves
                all_desc_leaf_are_target = False 
        else: # Clade is an internal node with descendant leaves, or it is a leaf itself.
            for leaf_descendant in descendant_leaves:
                # Check status of each descendant leaf. Default to False if name not in map.
                if not leaf_is_target_status.get(leaf_descendant.name, False):
                    all_desc_leaf_are_target = False
                    break
        
        # Only consider marking if all its descendant leaves are target.
        # Also, the node must either be a leaf itself, or if internal, must have had some leaves to check.
        # An empty internal node (no leaves below it) should not be marked.
        if all_desc_leaf_are_target and (clade.is_terminal() or descendant_leaves):
             nodes_to_be_marked.append(clade)

    # Step 3: Mark the identified clades by appending the PAML marker to their names
    for clade_to_mark in nodes_to_be_marked:
        current_name = clade_to_mark.name if clade_to_mark.name else ""
        
        # Avoid double-marking if the marker is already at the end of the name.
        if not current_name.endswith(paml_marker):
            clade_to_mark.name = current_name + paml_marker
        # If current_name was empty (e.g. unnamed internal node), it becomes just the paml_marker.
        # Bio.Phylo.write handles this correctly for Newick, e.g., (...)#1

    # Step 4: Save the modified tree
    base_name = os.path.basename(tree_path)
    name_part, ext_part = os.path.splitext(base_name)
    
    if not ext_part and base_name.lower().endswith(".treefile"):
         ext_part = ".treefile"
    elif not ext_part:
         ext_part = ".treefile" 

    output_filename = f"{name_part}_marked{ext_part}"
    output_path = os.path.join(output_folder, output_filename)

    try:
        Phylo.write(tree, output_path, "newick")
        # print(f"Saved marked tree to {output_path}")
    except Exception as e:
        print(f"Error writing marked tree {output_path} with BioPython: {e}")

def main():
    parser = argparse.ArgumentParser(
        description="Batch add PAML foreground markers to gene tree files (using BioPython) based on species sociality.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--tree_folder", required=True, help="Folder containing gene tree files (e.g., with .treefile extension).")
    parser.add_argument("--sociality_file", required=True, help="Path to a CSV or TSV file containing species sociality data.\nMust have 'Organism Name' and 'Social level' columns.")
    parser.add_argument("--output_folder", required=True, help="Folder where marked tree files will be saved.")
    parser.add_argument("--target_sociality", required=True, help="The sociality level to be marked as foreground (e.g., 'Advanced Eusocial').")
    parser.add_argument("--paml_marker", required=True, help="The PAML marker to append to foreground branches (e.g., '#1').")

    args = parser.parse_args()

    if not os.path.isdir(args.tree_folder):
        print(f"Error: Tree folder '{args.tree_folder}' not found.")
        return
    if not os.path.isfile(args.sociality_file):
        print(f"Error: Sociality file '{args.sociality_file}' not found.")
        return

    if not os.path.isdir(args.output_folder):
        try:
            os.makedirs(args.output_folder, exist_ok=True)
            print(f"Created output folder: {args.output_folder}")
        except OSError as e:
            print(f"Error: Could not create output folder '{args.output_folder}': {e}")
            return

    sociality_map = {}
    try:
        file_ext = os.path.splitext(args.sociality_file)[1].lower()
        if file_ext == ".csv":
            social_df = pd.read_csv(args.sociality_file)
        elif file_ext in [".tsv", ".txt"]:
            social_df = pd.read_csv(args.sociality_file, sep='\t')
        else:
            print(f"Warning: Unknown sociality file extension '{file_ext}'. Attempting to parse as CSV, then TSV.")
            try:
                social_df = pd.read_csv(args.sociality_file)
            except pd.errors.ParserError:
                social_df = pd.read_csv(args.sociality_file, sep='\t')
            except Exception as e_parse:
                 print(f"Error: Could not parse sociality file '{args.sociality_file}'. Error: {e_parse}")
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

    processed_count = 0
    tree_files_found = 0
    for filename in os.listdir(args.tree_folder):
        if filename.lower().endswith(".treefile"):
            tree_files_found +=1
            tree_file_path = os.path.join(args.tree_folder, filename)
            print(f"Processing {filename}...")
            # Call the BioPython version of the function
            _mark_tree_and_save_biopython(tree_file_path, sociality_map, args.target_sociality, args.paml_marker, args.output_folder)
            processed_count += 1
    
    if tree_files_found == 0:
        print(f"No '.treefile' files found in '{args.tree_folder}'.")
    else:
        print(f"\nSuccessfully processed {processed_count} out of {tree_files_found} '.treefile' files found.")
        print(f"Marked trees have been saved to: {args.output_folder}")

if __name__ == "__main__":
    main()