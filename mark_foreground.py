#!/usr/bin/env python3

import argparse
import os
import pandas as pd
from Bio import Phylo # 使用 BioPython

# Hardcoded list of species names, with underscores, sorted by length descending for longest match first
HARDCODED_SPECIES_LIST_WITH_UNDERSCORES = [
    "Acromyrmex_echinatior",
    "Apis_cerana",
    "Apis_dorsata",
    "Apis_florea",
    "Apis_laboriosa",
    "Apis_mellifera",
    "Atta_cephalotes",
    "Atta_colombica",
    "Bombus_affinis",
    "Bombus_bifarius",
    "Bombus_fervidus",
    "Bombus_flavifrons",
    "Bombus_huntii",
    "Bombus_impatiens",
    "Bombus_pascuorum",
    "Bombus_pyrosoma",
    "Bombus_terrestris",
    "Bombus_vancouverensis_nearcticus",
    "Bombus_vosnesenskii",
    "Camponotus_floridanus",
    "Cardiocondyla_obscurior",
    "Cataglyphis_hispanica",
    "Ceratina_calcarata",
    "Colletes_gigas",
    "Cyphomyrmex_costatus",
    "Dinoponera_quadriceps",
    "Drosophila_melanogaster",
    "Dufourea_novaeangliae",
    "Eufriesea_mexicana",
    "Formica_exsecta",
    "Frieseomelitta_varia",
    "Habropoda_laboriosa",
    "Harpegnathos_saltator",
    "Hylaeus_anthracinus",
    "Hylaeus_volcanicus",
    "Linepithema_humile",
    "Megachile_rotundata",
    "Megalopta_genalis",
    "Monomorium_pharaonis",
    "Nomia_melanderi",
    "Nylanderia_fulva",
    "Odontomachus_brunneus",
    "Ooceraea_biroi",
    "Osmia_bicornis_bicornis",
    "Osmia_lignaria",
    "Pogonomyrmex_barbatus",
    "Polistes_canadensis",
    "Polistes_dominula",
    "Polistes_fuscatus",
    "Polyergus_mexicanus",
    "Prorops_nasuta",
    "Pseudomyrmex_gracilis",
    "Solenopsis_invicta",
    "Temnothorax_curvispinosus",
    "Temnothorax_longispinosus",
    "Temnothorax_nylanderi",
    "Trachymyrmex_cornetzi",
    "Trachymyrmex_septentrionalis",
    "Trachymyrmex_zeteki",
    "Vespa_crabro",
    "Vespa_mandarinia",
    "Vespa_velutina",
    "Vespula_pensylvanica",
    "Vespula_vulgaris",
    "Vollenhovia_emeryi",
    "Wasmannia_auropunctata",
]
HARDCODED_SPECIES_LIST_WITH_UNDERSCORES.sort(key=len, reverse=True)

def _normalize_name_for_matching(name):
    """Converts name to a standard format for matching (lowercase, underscore to space)."""
    if name is None:
        return ""
    return name.replace('_', ' ').lower()

def _extract_species_key_from_gene_id(gene_id):
    """
    Attempts to extract a species key (e.g., 'genus_species') from a gene ID.
    Assumes gene ID might be like 'Genus_species_identifier' or 'Genus_species'.
    Now prioritizes matching against HARDCODED_SPECIES_LIST_WITH_UNDERSCORES.
    """
    if not gene_id:
        return None

    # Prioritize matching against the hardcoded list
    for known_species_key in HARDCODED_SPECIES_LIST_WITH_UNDERSCORES:
        # Check if gene_id starts with the known_species_key followed by an underscore,
        # or if gene_id is exactly the known_species_key (for cases where gene_id is just the species name)
        if gene_id.startswith(known_species_key + "_") or gene_id == known_species_key:
            return known_species_key # Return the version with underscores from the hardcoded list

    # Fallback: if no match from hardcoded list, try the old splitting logic (optional, or remove)
    # For now, let's stick to only hardcoded list matching as requested for more explicit control.
    # If a gene ID doesn't match, it will be treated as unknown later.
    # parts = gene_id.split('_')
    # if len(parts) >= 2:
    #     return f"{parts[0]}_{parts[1]}"
    # elif len(parts) == 1:
    #     return parts[0]
    
    return None # Return None if no match found in the hardcoded list


def _mark_tree_and_save_biopython(tree_path, normalized_sociality_map, target_sociality_normalized, paml_marker, output_folder):
    """
    Loads a tree using BioPython, maps gene IDs to species, marks nodes based on descendant sociality,
    renames leaves to species names, and saves the modified tree.
    """
    try:
        tree = Phylo.read(tree_path, "newick")
    except Exception as e:
        print(f"Warning: Could not parse tree {os.path.basename(tree_path)} with BioPython: {e}")
        return

    # Remove bootstrap values from all clades before any other processing
    for clade in tree.find_clades():
        clade.confidence = None

    # Step 1: Map leaf gene IDs to species names and their sociality status.
    # leaf_social_info stores: {original_gene_id: {'species_name_for_output': str, 'is_target': bool}}
    leaf_social_info = {}
    
    for leaf_node in tree.get_terminals():
        original_gene_id = leaf_node.name
        if not original_gene_id:
            print(f"Warning: Leaf node without a name found in tree '{os.path.basename(tree_path)}'. Skipping.")
            continue

        putative_species_key_from_gene = _extract_species_key_from_gene_id(original_gene_id)
        
        is_target_leaf = False
        # Default to original gene ID if no mapping, or if key extraction fails
        species_name_for_output = original_gene_id 

        if putative_species_key_from_gene:
            normalized_putative_species_key = _normalize_name_for_matching(putative_species_key_from_gene)
            
            social_data = normalized_sociality_map.get(normalized_putative_species_key)
            
            if social_data:
                species_name_for_output = social_data['original_name'].replace(' ', '_') # Use the name from CSV for output, replacing spaces with underscores
                if social_data['social_level_normalized'] == target_sociality_normalized:
                    is_target_leaf = True
                # print(f"Info: Gene '{original_gene_id}' (key: '{normalized_putative_species_key}') mapped to species '{species_name_for_output}', sociality: {social_data['social_level_normalized']}, target: {target_sociality_normalized}, is_target_leaf: {is_target_leaf}")
            else:
                print(f"Warning: Derived species key '{putative_species_key_from_gene}' (normalized: '{normalized_putative_species_key}') from gene '{original_gene_id}' in tree '{os.path.basename(tree_path)}' not found in sociality file. Treating as non-target.")
        else:
            print(f"Warning: Could not derive species key from gene ID '{original_gene_id}' in tree '{os.path.basename(tree_path)}'. Treating as non-target.")
            
        leaf_social_info[original_gene_id] = {
            'species_name_for_output': species_name_for_output,
            'is_target': is_target_leaf
        }

    # Step 2: Traverse all clades (nodes) and identify those whose ALL leaf descendants are of target_sociality.
    nodes_to_be_marked = []
    for clade in tree.find_clades(order="postorder"): 
        all_desc_leaf_are_target = True
        
        current_clade_descendant_leaves = list(clade.get_terminals())

        if not current_clade_descendant_leaves: # Node is an internal node with no leaf descendants (e.g. after pruning) or it's a mis-parsed leaf.
            if clade.is_terminal(): # Actually a leaf, but get_terminals() returned empty. This implies an issue or it's the only node.
                gene_id = clade.name 
                if not leaf_social_info.get(gene_id, {}).get('is_target', False):
                    all_desc_leaf_are_target = False
            else: # Genuinely an internal node with no leaves below it.
                all_desc_leaf_are_target = False
        else: # Node is a leaf itself or an internal node with leaf descendants.
            for leaf_descendant_node in current_clade_descendant_leaves:
                gene_id = leaf_descendant_node.name # At this point, leaf names are still original gene IDs
                if not leaf_social_info.get(gene_id, {}).get('is_target', False):
                    all_desc_leaf_are_target = False
                    break
        
        if all_desc_leaf_are_target and (clade.is_terminal() or current_clade_descendant_leaves):
             nodes_to_be_marked.append(clade)

    # Step 3: Rename ALL leaf nodes in the tree object to their species names.
    # This is done for all leaves, regardless of whether their branch is marked.
    for leaf_node in tree.get_terminals():
        original_gene_id = leaf_node.name 
        if original_gene_id in leaf_social_info: # Check if we have mapping info for this gene ID
            leaf_node.name = leaf_social_info[original_gene_id]['species_name_for_output']
        # If original_gene_id was not in leaf_social_info (e.g. unnamed leaf, or error in processing it),
        # its name remains as is (which might be the original gene ID or None).

    # Step 4: Mark the identified clades by appending the PAML marker.
    # If a marked clade is a leaf, its name has NOW been updated to species_name_for_output.
    for clade_to_mark in nodes_to_be_marked:
        current_name = clade_to_mark.name if clade_to_mark.name else ""
        
        if not current_name.endswith(paml_marker): # Avoid double-marking
            clade_to_mark.name = current_name + paml_marker

    # Step 5: Save the modified tree
    base_name = os.path.basename(tree_path)
    name_part, ext_part = os.path.splitext(base_name)
    
    if not ext_part and base_name.lower().endswith(".treefile"):
         ext_part = ".treefile"
    elif not ext_part: # Default to .treefile if no extension or different one
         ext_part = ".treefile" 

    output_filename = f"{name_part}_marked{ext_part}"
    output_path = os.path.join(output_folder, output_filename)

    try:
        Phylo.write(tree, output_path, "newick", format_branch_length=str)
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

    # Load sociality data
    # normalized_sociality_map stores: {normalized_species_key: {'original_name': str, 'social_level_normalized': str}}
    normalized_sociality_map = {}
    # Normalize target_sociality once for consistent comparison
    target_sociality_normalized = _normalize_name_for_matching(args.target_sociality)

    try:
        file_ext = os.path.splitext(args.sociality_file)[1].lower()
        # Add dtype=str to suppress DtypeWarning and ensure all columns are read as strings initially
        if file_ext == ".csv":
            social_df = pd.read_csv(args.sociality_file, dtype=str)
        elif file_ext in [".tsv", ".txt"]: # Allow .txt for TSV too
            social_df = pd.read_csv(args.sociality_file, sep='\t', dtype=str)
        else:
            print(f"Warning: Unknown sociality file extension '{file_ext}'. Attempting to parse as CSV (treating all columns as strings).")
            try:
                social_df = pd.read_csv(args.sociality_file, dtype=str)
            except pd.errors.ParserError:
                print(f"Warning: Failed to parse '{args.sociality_file}' as CSV. Attempting as TSV (treating all columns as strings).")
                social_df = pd.read_csv(args.sociality_file, sep='\t', dtype=str) # Ensure sep is correctly escaped for TSV
            except Exception as e_parse:
                 print(f"Error: Could not parse sociality file '{args.sociality_file}'. Error: {e_parse}")
                 return

        if "Organism Name" not in social_df.columns or "Social level" not in social_df.columns:
            print(f"Error: Sociality file '{args.sociality_file}' must contain 'Organism Name' and 'Social level' columns.")
            return
        
        for index, row in social_df.iterrows():
            # Explicitly convert to string, handle potential NaN from initial read if dtype=str wasn't perfect
            organism_name_original = str(row["Organism Name"]) if pd.notna(row["Organism Name"]) else ""
            social_level_original = str(row["Social level"]) if pd.notna(row["Social level"]) else ""

            # Skip rows where essential data is missing or effectively empty after conversion
            if not organism_name_original.strip() or not social_level_original.strip():
                # print(f"Warning: Skipping row {index+2} in '{args.sociality_file}' due to missing/empty 'Organism Name' or 'Social level'. Row data: {row.to_dict()}")
                continue
            
            normalized_key = _normalize_name_for_matching(organism_name_original)
            # Also normalize the social level from the file for comparison
            social_level_normalized = _normalize_name_for_matching(social_level_original)

            if normalized_key not in normalized_sociality_map:
                 normalized_sociality_map[normalized_key] = {
                     'original_name': organism_name_original.strip(), # Store the original, stripped name
                     'social_level_normalized': social_level_normalized
                 }
            # else:
                # print(f"Info: Duplicate normalized organism name key '{normalized_key}' (from '{organism_name_original}') found. Using first encountered: '{normalized_sociality_map[normalized_key]['original_name']}'.")

        if not normalized_sociality_map:
            print(f"Warning: No valid sociality data loaded from '{args.sociality_file}'. This may result in no branches being marked if species keys cannot be matched.")
        else:
            print(f"Loaded {len(normalized_sociality_map)} unique-normalized species sociality entries from '{args.sociality_file}'.")
            # print(f"Target sociality for matching (normalized): '{target_sociality_normalized}'")
            # if normalized_sociality_map:
            #     example_key = list(normalized_sociality_map.keys())[0]
            #     print(f"Example entry from sociality map - Key: '{example_key}', Value: {normalized_sociality_map[example_key]}")


    except Exception as e:
        print(f"Error loading or parsing sociality file '{args.sociality_file}': {e}")
        return

    processed_count = 0
    tree_files_found = 0
    for filename in os.listdir(args.tree_folder):
        if filename.lower().endswith(".treefile"): # Case-insensitive check
            tree_files_found +=1
            tree_file_path = os.path.join(args.tree_folder, filename)
            print(f"Processing {filename}...")
            _mark_tree_and_save_biopython(tree_file_path, normalized_sociality_map, target_sociality_normalized, args.paml_marker, args.output_folder)
            processed_count += 1
    
    if tree_files_found == 0:
        print(f"No '.treefile' files found in '{args.tree_folder}'.")
    else:
        print(f"\nSuccessfully processed {processed_count} out of {tree_files_found} '.treefile' files found.")
        print(f"Marked trees have been saved to: {args.output_folder}")

if __name__ == "__main__":
    main()