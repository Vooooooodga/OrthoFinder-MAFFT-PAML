#!/usr/bin/env python3

import argparse
import os
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

def get_species_from_single_fasta(fasta_path):
    """
    Scans a single FASTA file and returns a unique set of species names.
    """
    species_found = set()
    try:
        with open(fasta_path, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                species_key = _extract_species_key_from_gene_id(record.id)
                if species_key:
                    species_found.add(species_key)
                else:
                    print(f"警告: 无法为 '{os.path.basename(fasta_path)}' 中的序列 '{record.id}' 确定物种。")
    except Exception as e:
        print(f"读取FASTA文件 {fasta_path} 时出错: {e}")
    return species_found

def prune_and_save_subtree(original_tree, species_to_keep, output_path):
    """
    Prunes a given tree object to retain only the specified list of species and saves it.
    """
    # Work on a copy to not modify the original tree in the loop
    tree = original_tree.copy()
    
    all_tree_leaves = {leaf.name for leaf in tree.get_terminals()}
    
    # Identify which species from the FASTA file are actually present in the tree
    species_in_tree_to_keep = list(species_to_keep.intersection(all_tree_leaves))
    
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

    try:
        # Prune the tree to keep only the desired species
        tree.prune(species_in_tree_to_keep)
    except Exception as e:
        # Bio.Phylo.prune can raise an error if the resulting tree is trivial (e.g., only one leaf).
        # We handle this by creating a new tree with just that single node.
        if len(species_in_tree_to_keep) == 1:
             from Bio.Phylo.BaseTree import Tree as BioTree, Clade
             tree = BioTree(root=Clade(name=species_in_tree_to_keep[0]))
        else:
            print(f"修剪树时发生错误: {e}")
            return

    # Save the newly pruned tree to the output file
    try:
        Phylo.write(tree, output_path, "newick")
        print(f"成功创建并保存修剪后的子树到: {os.path.basename(output_path)}")
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
            
            species_from_fasta = get_species_from_single_fasta(fasta_path)
            
            if not species_from_fasta:
                print(f"未能在 '{filename}' 中找到任何可识别的物种。跳过。")
                continue
            
            # Define output path
            base_name = os.path.splitext(filename)[0]
            output_filename = f"{base_name}_subtree.treefile"
            output_path = os.path.join(args.output_dir, output_filename)

            # Prune and save
            prune_and_save_subtree(master_tree, species_from_fasta, output_path)
            processed_files += 1

    print(f"\n处理完成。在 '{args.fasta_dir}' 中找到 {total_files} 个FASTA文件，并成功为其中 {processed_files} 个生成了子树。")
    print(f"所有子树文件已保存至: {args.output_dir}")

if __name__ == "__main__":
    main() 