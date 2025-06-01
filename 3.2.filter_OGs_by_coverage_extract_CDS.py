import csv
import os
import argparse
from Bio import SeqIO

def build_index(nucleotide_db_dir):
    """
    为核苷酸数据库中的CDS序列构建索引。
    索引结构: {species_id: {gene_id: SeqRecord_longest_transcript}}
    """
    index = {} 
    print(f"Starting to build nucleotide index from: {nucleotide_db_dir}")
    if not os.path.isdir(nucleotide_db_dir):
        print(f"Error: Nucleotide DB directory does not exist: {nucleotide_db_dir}")
        return index

    for item_name in os.listdir(nucleotide_db_dir):
        item_path = os.path.join(nucleotide_db_dir, item_name)
        if os.path.isfile(item_path) and (item_name.endswith((".fa", ".fasta", ".fna"))):
            base_name = os.path.splitext(item_name)[0]
            # 从文件名提取物种ID, e.g., "Apis_mellifera" from "Apis_mellifera.cds.fa"
            species_id_from_filename = base_name.split('.')[0] 

            record_dict = {}
            try:
                with open(item_path, "r") as handle:
                    for record in SeqIO.parse(handle, "fasta"):
                        key = record.id # 基因ID，假设FASTA头为 >gene_id
                        if key in record_dict:
                            if len(record.seq) > len(record_dict[key].seq):
                                record_dict[key] = record
                        else:
                            record_dict[key] = record
                if record_dict:
                    index[species_id_from_filename] = record_dict
                # else:
                #     print(f"No records found or indexed for species '{species_id_from_filename}' from file '{item_name}'.")
            except Exception as e:
                print(f"Error parsing or indexing file {item_path}: {e}")
        # else:
            # print(f"Skipping non-FASTA file or directory in nucleotide_db_dir: {item_name}")
    
    if not index:
        print(f"Warning: Nucleotide index is empty. Check content of {nucleotide_db_dir}")
    return index

def get_species_gene_map_from_orthogroups(orthogroups_file, orthofinder_header_species_names_to_exclude):
    """
    解析 Orthogroups.tsv 文件。
    返回:
    1. orthogroup_to_species_genes: dict[str, dict[str, list[str]]]
       OG_ID -> {species_name (from header) -> [gene_id1, gene_id2]}
    2. relevant_species_names_in_header: list[str]
       在 OrthoFinder TSV 文件头中找到的、未被排除的物种名称列表。
    """
    orthogroup_to_species_genes = {}
    relevant_species_column_map = {} # species_name_from_header -> column_index

    if not os.path.isfile(orthogroups_file):
        print(f"Error: Orthogroups file not found: {orthogroups_file}")
        return {}, []

    with open(orthogroups_file, newline='') as csvfile:
        tsv_reader = csv.reader(csvfile, delimiter='\t')
        try:
            header = next(tsv_reader)
        except StopIteration:
            print(f"Error: Orthogroups file is empty or not a valid TSV: {orthogroups_file}")
            return {}, []
        
        # OrthoFinder 输出中非物种的标准列名
        non_species_columns = {"HOGs", "Orthogroup", "Gene Tree Parent Clade"}
        
        for i, col_name in enumerate(header):
            if i == 0 and col_name == "Orthogroup": # 第一列通常是Orthogroup ID
                continue
            elif i == 0 and col_name != "Orthogroup":
                 print(f"Warning: First column of Orthogroups.tsv is '{col_name}', expected 'Orthogroup'. Assuming it's the Orthogroup ID column.")
            
            # 如果列名不是标准非物种列，并且不在排除列表中，则视为相关物种列
            if col_name not in non_species_columns and \
               col_name not in orthofinder_header_species_names_to_exclude:
                relevant_species_column_map[col_name] = i
        
        relevant_species_names_in_header = list(relevant_species_column_map.keys())
        if not relevant_species_names_in_header:
            print("Warning: No relevant species columns identified in Orthogroups.tsv header after exclusions. Check file format and exclusion list.")

        for row_num, row in enumerate(tsv_reader):
            if not row or not row[0]: # 跳过空行或没有Orthogroup ID的行
                continue
            
            orthogroup_id = row[0]
            species_genes_for_current_og = {}
            for species_name, col_idx in relevant_species_column_map.items():
                if col_idx < len(row) and row[col_idx]: # 确保列索引有效且单元格非空
                    genes_str = row[col_idx]
                    # 基因以 ", " 分隔
                    genes = [g.strip() for g in genes_str.split(', ') if g.strip()]
                    species_genes_for_current_og[species_name] = genes
                else:
                    species_genes_for_current_og[species_name] = [] # 物种在该OG中没有基因或列数据缺失
            
            orthogroup_to_species_genes[orthogroup_id] = species_genes_for_current_og
            
    return orthogroup_to_species_genes, relevant_species_names_in_header

def main():
    parser = argparse.ArgumentParser(
        description="Filters Orthogroups based on species single-copy coverage and extracts longest CDS for qualifying genes. "
                    "Assumes gene IDs in Orthogroups.tsv, MSA headers, and CDS FASTA files are consistent or handleable via common OrthoFinder modifications (e.g., '()' to '_').",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--orthogroups_file", required=True, help="Path to the Orthogroups.tsv file from OrthoFinder.")
    parser.add_argument("--msa_dir_from_script1", required=True, help="Directory containing MSA FASTA files (output of script 1), where filenames are Orthogroup IDs (e.g., OG0000001.fa).")
    parser.add_argument("--nucleotide_db_dir", required=True, help="Directory containing nucleotide (CDS) FASTA files. \nEach file should be for one species (e.g., Apis_mellifera.cds.fa), \nand FASTA headers within should be '>gene_id'.")
    parser.add_argument("--output_dir", default='cds_filtered_by_coverage', help="Directory to save extracted CDS sequences for qualified Orthogroups.")
    parser.add_argument("--excluded_species", nargs='*', default=['Drosophila_melanogaster'], help="Species names (matching Orthogroups.tsv header) to exclude from coverage calculation and CDS extraction.")
    parser.add_argument("--single_copy_threshold_percentage", type=float, default=70.0, help="Minimum percentage of (non-excluded) species that must have a single-copy gene in an Orthogroup for it to be processed.")
    parser.add_argument("--subspecies_included", nargs='*', default=['Bombus_vancouverensis_nearcticus', 'Osmia_bicornis_bicornis'], 
                        help="Full names of species (including subspecies, e.g., Genus_species_subspecies) if their MSA headers use these full names. \nThese names must also match the species names in Orthogroups.tsv and the filenames in nucleotide_db_dir if they are to be processed correctly.")
    
    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    print("Step 1: Building nucleotide index from CDS files...")
    nucleotide_index = build_index(args.nucleotide_db_dir)
    if not nucleotide_index:
        print("Failed to build nucleotide index or index is empty. Cannot proceed. Please check --nucleotide_db_dir.")
        return
    print(f"Nucleotide index built for {len(nucleotide_index)} species.")

    print("\nStep 2: Parsing Orthogroups.tsv to map genes to species within orthogroups...")
    orthogroup_to_species_genes_map, tsv_header_species_for_calc = get_species_gene_map_from_orthogroups(args.orthogroups_file, args.excluded_species)
    if not orthogroup_to_species_genes_map:
        print("Failed to parse Orthogroups.tsv or no orthogroups found. Cannot proceed.")
        return
    print(f"Parsed {len(orthogroup_to_species_genes_map)} orthogroups. Using {len(tsv_header_species_for_calc)} species for coverage calculation: {tsv_header_species_for_calc if len(tsv_header_species_for_calc) < 10 else str(tsv_header_species_for_calc[:10]) + '...'}")

    relevant_species_count_for_threshold = len(tsv_header_species_for_calc)
    if relevant_species_count_for_threshold == 0:
        print("Error: No relevant species found for coverage calculation after exclusion. Check --orthogroups_file header and --excluded_species.")
        return
    
    print(f"\nStep 3: Filtering orthogroups based on single-copy species coverage (Threshold: >{args.single_copy_threshold_percentage}% of {relevant_species_count_for_threshold} species)...")
    qualified_orthogroups = {} # OG_ID -> {species_name_from_tsv_header: single_gene_id}

    for og_id, species_to_genes_dict in orthogroup_to_species_genes_map.items():
        single_copy_species_this_og_count = 0
        for species_name, gene_list in species_to_genes_dict.items(): # species_to_genes_dict 只包含未排除的物种
            if len(gene_list) == 1:
                single_copy_species_this_og_count += 1
        
        percentage_single_copy = (single_copy_species_this_og_count / relevant_species_count_for_threshold) * 100 if relevant_species_count_for_threshold > 0 else 0

        if percentage_single_copy > args.single_copy_threshold_percentage:
            single_copy_genes_for_this_qualified_og = {}
            for species_name, gene_list in species_to_genes_dict.items():
                if len(gene_list) == 1:
                    single_copy_genes_for_this_qualified_og[species_name] = gene_list[0]
            qualified_orthogroups[og_id] = single_copy_genes_for_this_qualified_og
    
    print(f"Found {len(qualified_orthogroups)} orthogroups meeting the single-copy coverage criteria.")
    if not qualified_orthogroups:
        print("No orthogroups met the filtering criteria. Exiting.")
        return

    print("\nStep 4: Extracting CDS sequences for qualified orthogroups and their single-copy genes...")
    final_cds_files_created = 0
    msa_files_processed = 0
    for msa_filename in os.listdir(args.msa_dir_from_script1):
        if not (msa_filename.endswith(('.fa', '.fasta'))):
            continue
        msa_files_processed +=1
        current_og_id = os.path.splitext(msa_filename)[0]

        if current_og_id in qualified_orthogroups:
            single_copy_genes_map_for_current_og = qualified_orthogroups[current_og_id]
            
            msa_file_path = os.path.join(args.msa_dir_from_script1, msa_filename)
            output_cds_file_path = os.path.join(args.output_dir, f"{current_og_id}_cds.fa")
            
            sequences_written_to_this_cds_file = 0
            with open(msa_file_path, "r") as msa_handle, open(output_cds_file_path, "w") as cds_output_handle:
                for msa_record in SeqIO.parse(msa_handle, "fasta"):
                    msa_header_id = msa_record.id 
                    
                    parsed_species_name_from_msa = ""
                    is_subspecies_match = False
                    for subspecies_full_name in args.subspecies_included:
                        if msa_header_id.startswith(subspecies_full_name + '_'): 
                            parsed_species_name_from_msa = subspecies_full_name
                            is_subspecies_match = True
                            break
                    
                    if not is_subspecies_match:
                        parts = msa_header_id.split('_')
                        if len(parts) >= 2: 
                            potential_species_name = '_'.join(parts[:2]) # 例如 Apis_mellifera
                            # 检查此解析出的物种名是否存在于我们的目标物种列表或CDS索引中
                            if potential_species_name in single_copy_genes_map_for_current_og or \
                               potential_species_name in nucleotide_index:
                                parsed_species_name_from_msa = potential_species_name
                        
                    if not parsed_species_name_from_msa:
                        # print(f"Warning: Could not determine species name from MSA header '{msa_header_id}' in OG '{current_og_id}'. Skipping this sequence.")
                        continue

                    if parsed_species_name_from_msa in single_copy_genes_map_for_current_og:
                        authoritative_gene_id = single_copy_genes_map_for_current_og[parsed_species_name_from_msa]
                        
                        species_cds_records = nucleotide_index.get(parsed_species_name_from_msa)
                        if species_cds_records:
                            target_cds_seq_record = species_cds_records.get(authoritative_gene_id)
                            
                            if not target_cds_seq_record: # 尝试处理OrthoFinder对基因名的修改 (e.g. "()" -> "_")
                                for original_cds_db_key, record_in_db in species_cds_records.items():
                                    # 模拟OrthoFinder对原始CDS库基因名的修改方式
                                    mangled_cds_db_key = original_cds_db_key.replace('(', '_').replace(')', '_')
                                    if mangled_cds_db_key == authoritative_gene_id or \
                                       mangled_cds_db_key.rstrip('_') == authoritative_gene_id: # 处理末尾可能的额外"_"
                                        target_cds_seq_record = record_in_db
                                        break
                            
                            if target_cds_seq_record:
                                output_fasta_header = msa_record.description if msa_record.description else msa_record.id
                                cds_output_handle.write(f">{output_fasta_header}\n{str(target_cds_seq_record.seq)}\n")
                                sequences_written_to_this_cds_file += 1
                            # else:
                                # print(f"Warning: CDS sequence not found for gene '{authoritative_gene_id}' of species '{parsed_species_name_from_msa}' in OG '{current_og_id}'. MSA header was '{msa_header_id}'.")
                        # else:
                            # print(f"Warning: No CDS records found in nucleotide_index for species '{parsed_species_name_from_msa}' (from MSA header '{msa_header_id}') for OG '{current_og_id}'.")

            if sequences_written_to_this_cds_file > 0:
                final_cds_files_created += 1
            else: # 如果没有序列写入，删除空的输出文件
                if os.path.exists(output_cds_file_path):
                    try:
                        if os.path.getsize(output_cds_file_path) == 0:
                            os.remove(output_cds_file_path)
                    except OSError as e:
                        print(f"Error removing empty file {output_cds_file_path}: {e}")
    
    print(f"\n--- Summary ---")
    print(f"Processed {msa_files_processed} MSA files from '{args.msa_dir_from_script1}'.")
    print(f"Created {final_cds_files_created} CDS FASTA files in '{args.output_dir}'.")
    print("Script execution finished.")

if __name__ == '__main__':
    main() 