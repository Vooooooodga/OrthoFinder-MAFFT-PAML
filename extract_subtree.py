#!/usr/bin/env python3

import argparse
import os
import copy
from Bio import Phylo
from Bio import SeqIO

# This hardcoded list is crucial for correctly parsing species names, 
# especially those with subspecies like 'Bombus_vancouverensis_nearcticus'.
# It is sorted by length in descending order to ensure the longest possible match is found first.
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
   # "Drosophila_melanogaster",
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

def _extract_species_key_from_gene_id(gene_id):
    """
    Extracts a species key (e.g., 'genus_species') from a gene ID
    by matching against the hardcoded list of species names.
    """
    if not gene_id:
        return None

    for known_species_key in HARDCODED_SPECIES_LIST_WITH_UNDERSCORES:
        # Check if the gene_id is an exact match or starts with the species name followed by an underscore.
        if gene_id.startswith(known_species_key + "_") or gene_id == known_species_key:
            return known_species_key
    
    return None # Return None if no match is found

def get_species_to_gene_map_from_fasta(fasta_path):
    """
    Scans a single FASTA file and returns a map of {species_key: full_gene_id}.
    If a species has multiple genes, it uses the first one found and prints a warning.
    """
    species_map = {}
    try:
        with open(fasta_path, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                species_key = _extract_species_key_from_gene_id(record.id)
                if species_key:
                    if species_key not in species_map:
                        species_map[species_key] = record.id
                    else:
                        print(f"警告: 物种 '{species_key}' 在文件 '{os.path.basename(fasta_path)}' 中存在多个基因序列。"
                              f"将使用第一个找到的序列 '{species_map[species_key]}'。")
                else:
                    print(f"警告: 无法为 '{os.path.basename(fasta_path)}' 中的序列 '{record.id}' 确定物种。")
    except Exception as e:
        print(f"读取FASTA文件 {fasta_path} 时出错: {e}")
    return species_map

def prune_and_save_subtree(original_tree, species_to_gene_map, output_path):
    """
    Prunes a given tree object, renames leaves to full gene IDs, and saves it.
    This version uses a more robust method of pruning by removing unwanted leaves one by one.
    """
    # Work on a copy to not modify the original tree in the loop
    tree = copy.deepcopy(original_tree)
    
    all_tree_leaves = {leaf.name for leaf in tree.get_terminals()}
    
    species_to_keep = set(species_to_gene_map.keys())

    # Identify which species from the FASTA file are actually present in the tree
    species_in_tree_to_keep = species_to_keep.intersection(all_tree_leaves)
    
    # Report species that were in FASTA files but not in the tree
    species_not_in_tree = species_to_keep.difference(all_tree_leaves)
    if species_not_in_tree:
        print(f"信息: 在FASTA文件中找到但在主树中缺失以下 {len(species_not_in_tree)} 个物种，它们将被忽略:")
        for sp in sorted(list(species_not_in_tree)):
             print(f"- {sp}")

    if not species_in_tree_to_keep:
        print("错误: 此FASTA文件中的所有物种均未在提供的树文件中找到。无法创建子树。")
        return

    print(f"在主树中找到 {len(species_in_tree_to_keep)} 个匹配的物种用于创建子树。")

    # New pruning strategy: iterate through all leaves and remove the ones not in our target list.
    leaves_to_remove = all_tree_leaves - species_in_tree_to_keep
    
    try:
        for leaf_name in leaves_to_remove:
            tree.prune(leaf_name)
    except Exception as e:
        print(f"修剪树时发生错误: {e}")
        return

    # --- New Step: Rename leaves from species name to full gene ID ---
    renamed_count = 0
    for leaf in tree.get_terminals():
        if leaf.name in species_to_gene_map:
            leaf.name = species_to_gene_map[leaf.name]
            renamed_count += 1
        else:
            # This should not happen if pruning was successful, but as a safeguard:
            print(f"警告: 在重命名叶节点时, 树中的叶 '{leaf.name}' 未在物种-基因映射中找到。")
    print(f"已成功将 {renamed_count} 个叶节点重命名为完整的基因ID。")

    # Save the newly pruned and renamed tree to the output file
    try:
        Phylo.write(tree, output_path, "newick")
        print(f"成功创建并保存修剪后的子树到: {os.path.basename(output_path)}")

        # Per user request, forcefully remove wrapping single quotes that Bio.Phylo adds to leaf names.
        # This can potentially create non-standard Newick files if names have special characters.
        with open(output_path, 'r') as file:
            content = file.read()
        
        modified_content = content
        made_changes = False

        # Iterate through the actual leaf names in our final tree
        for leaf in tree.get_terminals():
            leaf_name = leaf.name
            # Bio.Phylo quotes names with special characters by wrapping them in single quotes
            # and escaping internal single quotes by doubling them (e.g., O'Malley -> 'O''Malley').
            quoted_name = "'" + leaf_name.replace("'", "''") + "'"
            
            # We perform a targeted replacement of the quoted version with the original.
            if quoted_name in modified_content:
                modified_content = modified_content.replace(quoted_name, leaf_name)
                made_changes = True
        
        if made_changes:
            with open(output_path, 'w') as file:
                file.write(modified_content)
            print(f"信息: 已根据要求从 {os.path.basename(output_path)} 中移除了叶节点名称的外部单引号。")

    except Exception as e:
        print(f"将修剪后的树写入 {output_path} 时出错: {e}")


def main():
    parser = argparse.ArgumentParser(
        description="为目录中的每个FASTA文件，从一个大型物种树中提取一个对应的子树。",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--tree", required=True, help="输入的物种树文件路径 (Newick格式)。")
    parser.add_argument("--fasta_dir", required=True, help="包含FASTA文件的目录路径。")
    parser.add_argument("--output_dir", required=True, help="输出目录，用于存放所有生成的子树文件。")

    args = parser.parse_args()

    # Create output directory
    if not os.path.isdir(args.output_dir):
        try:
            os.makedirs(args.output_dir)
            print(f"已创建输出目录: {args.output_dir}")
        except OSError as e:
            print(f"错误: 无法创建输出目录 '{args.output_dir}': {e}")
            return

    # Load master tree once
    try:
        master_tree = Phylo.read(args.tree, "newick")
        print(f"主树 '{os.path.basename(args.tree)}' 加载成功。")
    except Exception as e:
        print(f"错误: 无法解析主树文件 {args.tree}: {e}")
        return

    # Loop through fasta files
    processed_files = 0
    total_files = 0
    for filename in os.listdir(args.fasta_dir):
        if filename.lower().endswith(('.fasta', '.fas', '.fa', '.faa', '.fna')):
            total_files += 1
            fasta_path = os.path.join(args.fasta_dir, filename)
            print(f"\n--- 正在处理: {filename} ---")
            
            species_map = get_species_to_gene_map_from_fasta(fasta_path)
            
            if not species_map:
                print(f"未能在 '{filename}' 中找到任何可识别的物种。跳过。")
                continue
            
            # Define output path
            base_name = os.path.splitext(filename)[0]
            output_filename = f"{base_name}_subtree.treefile"
            output_path = os.path.join(args.output_dir, output_filename)

            # Prune, rename, and save
            prune_and_save_subtree(master_tree, species_map, output_path)
            processed_files += 1

    print(f"\n处理完成。在 '{args.fasta_dir}' 中找到 {total_files} 个FASTA文件，并成功为其中 {processed_files} 个生成了子树。")
    print(f"所有子树文件已保存至: {args.output_dir}")

if __name__ == "__main__":
    main() 