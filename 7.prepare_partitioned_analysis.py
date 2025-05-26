import os
import re
import sys
from collections import defaultdict

def get_best_model_from_iqtree_log(iqtree_file_path):
    """
    Parses an IQ-TREE log file (.iqtree) to find the best model.
    Looks for "Best-fit model according to BIC:"
    """
    try:
        with open(iqtree_file_path, 'r') as f:
            for line in f:
                if "Best-fit model according to BIC:" in line:
                    model = line.split("Best-fit model according to BIC:")[1].strip()
                    return model
    except FileNotFoundError:
        print(f"警告: IQ-TREE 日志文件未找到: {iqtree_file_path}", file=sys.stderr)
    except Exception as e:
        print(f"警告: 解析 IQ-TREE 日志文件 '{iqtree_file_path}' 时出错: {e}", file=sys.stderr)
    return None

def read_fasta(fasta_file_path):
    """
    Reads a FASTA file and returns a dictionary of {header: sequence}
    and the length of the first sequence (assuming an alignment).
    Also returns a list of headers in order of appearance.
    """
    sequences = {}
    ordered_headers = []
    seq_len = 0
    try:
        with open(fasta_file_path, 'r') as f:
            header = None
            current_seq = []
            first_seq_done = False
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if header: # Process previous sequence
                        sequences[header] = "".join(current_seq)
                        if not first_seq_done:
                            seq_len = len(sequences[header])
                            first_seq_done = True
                    header = line[1:] # Remove ">"
                    ordered_headers.append(header)
                    current_seq = []
                else:
                    current_seq.append(line)
            if header: # Process the last sequence
                sequences[header] = "".join(current_seq)
                if not first_seq_done:
                    seq_len = len(sequences[header])
    except FileNotFoundError:
        print(f"警告: FASTA 文件未找到: {fasta_file_path}", file=sys.stderr)
        return None, 0, None
    except Exception as e:
        print(f"警告: 读取 FASTA 文件 '{fasta_file_path}' 时出错: {e}", file=sys.stderr)
        return None, 0, None
    return sequences, seq_len, ordered_headers

def main(msa_dir, iqtree_results_dir, concatenated_msa_out, partition_file_out):
    """
    Main function to generate concatenated MSA and partition file.
    """
    print(f"扫描 MSA 目录: {msa_dir}")
    print(f"扫描 IQ-TREE 结果目录: {iqtree_results_dir}")

    og_files = {} # Store paths to msa and iqtree files for each OG

    # Find all potential OGs from MSA directory
    for f_name in os.listdir(msa_dir):
        if f_name.startswith("OG") and f_name.endswith(".clipkit.fasta"):
            og_id = f_name.split(".clipkit.fasta")[0]
            if og_id not in og_files:
                og_files[og_id] = {}
            og_files[og_id]['msa'] = os.path.join(msa_dir, f_name)

    # Match with iqtree results
    for f_name in os.listdir(iqtree_results_dir):
        if f_name.startswith("OG") and f_name.endswith(".clipkit.iqtree"):
            og_id = f_name.split(".clipkit.iqtree")[0]
            if og_id in og_files: # Only consider if MSA exists
                 og_files[og_id]['iqtree_log'] = os.path.join(iqtree_results_dir, f_name)


    concatenated_sequences = defaultdict(str)
    partitions = []
    current_msa_start_pos = 1
    all_species_ordered = None # To maintain species order from the first valid OG
    
    valid_ogs_count = 0

    # Sort OGs for consistent processing order (optional, but good practice)
    sorted_og_ids = sorted(og_files.keys())

    for og_id in sorted_og_ids:
        og_data = og_files[og_id]
        if 'msa' not in og_data or 'iqtree_log' not in og_data:
            print(f"跳过 {og_id}: 缺少 MSA 或 IQ-TREE 日志文件。", file=sys.stderr)
            continue
        
        print(f"处理 {og_id}...")

        # 1. Extract model
        model = get_best_model_from_iqtree_log(og_data['iqtree_log'])
        if not model:
            print(f"跳过 {og_id}: 无法从 {og_data['iqtree_log']} 提取模型。", file=sys.stderr)
            continue

        # 2. Read MSA, get length and sequences
        sequences, msa_len, current_og_species_ordered = read_fasta(og_data['msa'])
        if not sequences or msa_len == 0:
            print(f"跳过 {og_id}: 无法读取 MSA 文件或 MSA 为空: {og_data['msa']}", file=sys.stderr)
            continue
        
        # Initialize species order and concatenated_sequences dictionary with all species from the first valid OG
        if all_species_ordered is None:
            all_species_ordered = current_og_species_ordered
            for sp_header in all_species_ordered:
                concatenated_sequences[sp_header] = "" # Initialize, will be filled with gaps if needed

        # 3. Concatenate sequences
        for sp_header in all_species_ordered:
            if sp_header in sequences:
                concatenated_sequences[sp_header] += sequences[sp_header]
            else:
                # Add gaps if species is missing in this specific MSA
                concatenated_sequences[sp_header] += "-" * msa_len
        
        # 4. Store partition info
        partition_end_pos = current_msa_start_pos + msa_len - 1
        partitions.append(f"{model}, {og_id} = {current_msa_start_pos}-{partition_end_pos}")
        current_msa_start_pos = partition_end_pos + 1
        valid_ogs_count += 1

    if valid_ogs_count == 0:
        print("错误: 未找到任何有效的 OG 进行处理。请检查输入目录和文件命名。", file=sys.stderr)
        return

    # 5. Write concatenated MSA
    print(f"\n正在写入合并后的 MSA 到: {concatenated_msa_out} ...")
    with open(concatenated_msa_out, 'w') as f_out:
        for sp_header in all_species_ordered:
            f_out.write(f">{sp_header}\n")
            # Write sequence in chunks (e.g., 60 chars per line)
            seq = concatenated_sequences[sp_header]
            for i in range(0, len(seq), 60):
                f_out.write(seq[i:i+60] + "\n")
    print("合并后的 MSA 文件写入完成。")

    # 6. Write partition file
    print(f"\n正在写入分区文件到: {partition_file_out} ...")
    with open(partition_file_out, 'w') as f_out:
        for part_info in partitions:
            f_out.write(part_info + "\n")
    print("分区文件写入完成。")
    print(f"\n成功处理了 {valid_ogs_count} 个 Orthogroups。")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("用法: python prepare_partitioned_analysis.py <AA_MSA目录> <IQ-TREE结果目录> <输出合并MSA文件名> <输出分区文件名>")
        print("示例: python prepare_partitioned_analysis.py SCOGs_msa_aa_clipkit gene_trees_with_aa_msa_MFP concatenated_aa.fasta partitions_aa.txt")
        sys.exit(1)
    
    msa_directory = sys.argv[1]
    iqtree_results_directory = sys.argv[2]
    concatenated_msa_output_file = sys.argv[3]
    partition_output_file = sys.argv[4]
    
    main(msa_directory, iqtree_results_directory, concatenated_msa_output_file, partition_output_file)