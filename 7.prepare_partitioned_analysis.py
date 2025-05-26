import os
import re
import sys
from collections import defaultdict

# --- 筛选阈值 ---
MIN_MSA_LEN = 100        # 最小比对长度 (AA)
MAX_BIAS_RATIO = 0.20    # 成分检验失败的最大比例 (例如 0.20 = 20%)
LOW_RATE_PERCENTILE = 5  # 要移除的最慢进化速率的百分比
HIGH_RATE_PERCENTILE = 95 # 要移除的最快进化速率的百分比
# --------------------

def get_best_model_from_iqtree_log(iqtree_file_path):
    """
    解析 IQ-TREE 日志文件 (.iqtree) 以找到最佳模型。
    查找 "Best-fit model according to BIC:"
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
        print(f"警告: 解析 IQ-TREE 日志文件 '{iqtree_file_path}' 时出错 (Model): {e}", file=sys.stderr)
    return None

def get_tree_length_from_iqtree_log(iqtree_file_path):
    """
    解析 IQ-TREE 日志文件 (.iqtree) 以找到总树长。
    """
    try:
        with open(iqtree_file_path, 'r') as f:
            for line in f:
                # 修改: 使用正则表达式以适应更灵活的格式
                match = re.search(r"Total tree length.*:\\s*(\\d+\\.?\\d*)", line)
                if match:
                    try:
                        length = float(match.group(1)) # 从匹配组中提取长度
                        return length
                    except ValueError:
                         print(f"警告: 无法从 '{line.strip()}' 解析树长于 {iqtree_file_path}", file=sys.stderr)
                         return None
    except FileNotFoundError: # 确保文件未找到时有明确提示
        print(f"警告: IQ-TREE 日志文件未找到: {iqtree_file_path}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"警告: 解析 IQ-TREE 日志文件 '{iqtree_file_path}' 时出错 (TreeLen): {e}", file=sys.stderr)
    return None

def get_compositional_bias_stats(iqtree_file_path):
    """
    解析 IQ-TREE 日志文件以找到成分检验失败的序列数。
    查找 "NUMBER OF SEQUENCES FAILED:" 行。
    返回一个元组 (failed_count, total_count) 或 (None, None)。
    """
    try:
        with open(iqtree_file_path, 'r') as f:
            for line in f:
                if "NUMBER OF SEQUENCES FAILED:" in line:
                    match = re.search(r'(\d+)\s*/\s*(\d+)', line)
                    if match:
                        failed = int(match.group(1))
                        total = int(match.group(2))
                        return failed, total
                    else:
                        print(f"警告: 找到了失败行但无法解析: '{line.strip()}' 于 {iqtree_file_path}", file=sys.stderr)
                        return None, None
        # 如果未找到该行
        print(f"警告: 未在 {iqtree_file_path} 中找到 'NUMBER OF SEQUENCES FAILED' 行。", file=sys.stderr)
        return None, None
    except Exception as e:
        print(f"警告: 解析 IQ-TREE 日志文件 '{iqtree_file_path}' 时出错 (Bias): {e}", file=sys.stderr)
    return None, None

def read_fasta(fasta_file_path):
    """
    读取 FASTA 文件并返回 {header: sequence} 字典、
    第一个序列的长度以及按出现顺序排列的头列表。
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
                    if header:
                        sequences[header] = "".join(current_seq)
                        if not first_seq_done:
                            seq_len = len(sequences[header])
                            first_seq_done = True
                    header = line[1:]
                    ordered_headers.append(header)
                    current_seq = []
                else:
                    current_seq.append(line)
            if header:
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

def get_percentile(data, percentile_val):
    """计算给定数据的百分位数（线性插值法）。"""
    if not data:
        return None
    size = len(data)
    sorted_data = sorted(data)
    # k 是 1-based 索引，转换为 0-based
    k = (size - 1) * percentile_val / 100.0
    
    if k <= 0: return sorted_data[0]
    if k >= size - 1: return sorted_data[-1]

    f = int(k) # floor
    c = f + 1 if f < size - 1 else f # ceiling
    
    if f == c: # 整数索引
        return sorted_data[f]
    else: # 插值
        return sorted_data[f] + (sorted_data[c] - sorted_data[f]) * (k - f)


def main(msa_dir, iqtree_results_dir, concatenated_msa_out, partition_file_out):
    """
    主函数，用于收集数据、筛选基因并生成串联 MSA 和分区文件。
    """
    print(f"扫描 MSA 目录: {msa_dir}")
    print(f"扫描 IQ-TREE 结果目录: {iqtree_results_dir}")

    og_files = {}
    for f_name in os.listdir(msa_dir):
        if f_name.startswith("OG") and f_name.endswith(".clipkit.fasta"):
            og_id = f_name.split(".clipkit.fasta")[0]
            if og_id not in og_files: og_files[og_id] = {}
            og_files[og_id]['msa'] = os.path.join(msa_dir, f_name)

    for f_name in os.listdir(iqtree_results_dir):
        if f_name.startswith("OG") and f_name.endswith(".clipkit.iqtree"):
            og_id = f_name.split(".clipkit.iqtree")[0]
            if og_id in og_files:
                og_files[og_id]['iqtree_log'] = os.path.join(iqtree_results_dir, f_name)

    # --- 1. 数据收集 ---
    print("\n--- 1. 正在收集基因数据 ---")
    og_stats = []
    all_species_set = set()
    
    for og_id in sorted(og_files.keys()):
        og_data = og_files[og_id]
        if 'msa' not in og_data or 'iqtree_log' not in og_data:
            print(f"跳过 {og_id}: 缺少 MSA 或 IQ-TREE 日志。", file=sys.stderr)
            continue

        print(f"  正在处理 {og_id}...")
        sequences, msa_len, current_og_species = read_fasta(og_data['msa'])
        if not sequences or msa_len == 0:
            print(f"  -> 跳过 {og_id}: MSA 为空或读取失败。", file=sys.stderr)
            continue

        model = get_best_model_from_iqtree_log(og_data['iqtree_log'])
        if not model:
            print(f"  -> 跳过 {og_id}: 无法提取模型。", file=sys.stderr)
            continue

        tree_len = get_tree_length_from_iqtree_log(og_data['iqtree_log'])
        if tree_len is None:
            print(f"  -> 跳过 {og_id}: 无法提取树长。", file=sys.stderr)
            continue

        failed, total = get_compositional_bias_stats(og_data['iqtree_log'])
        if failed is None:
            print(f"  -> 跳过 {og_id}: 无法提取成分偏好统计。", file=sys.stderr)
            continue

        bias_ratio = failed / total if total > 0 else 1.0

        og_stats.append({
            'id': og_id, 'msa': og_data['msa'], 'iqtree': og_data['iqtree_log'],
            'len': msa_len, 'tree_len': tree_len, 'bias': bias_ratio,
            'model': model
        })
        all_species_set.update(sequences.keys()) # 使用实际读到的物种

    print(f"成功收集了 {len(og_stats)} 个基因的数据。")
    all_species_ordered = sorted(list(all_species_set)) # 确定所有物种的最终顺序
    print(f"共发现 {len(all_species_ordered)} 个物种。")

    # --- 2. 基因筛选 ---
    print("\n--- 2. 正在筛选基因 ---")
    initial_count = len(og_stats)
    
    # 筛选 1: 比对长度
    og_stats_filtered = [og for og in og_stats if og['len'] >= MIN_MSA_LEN]
    print(f"比对长度筛选 ({MIN_MSA_LEN} AA): {initial_count} -> {len(og_stats_filtered)}")

    # 筛选 2: 成分偏好
    og_stats_filtered = [og for og in og_stats_filtered if og['bias'] <= MAX_BIAS_RATIO]
    print(f"成分偏好筛选 (<={MAX_BIAS_RATIO*100}%): {len(og_stats_filtered)} -> {len(og_stats_filtered)}")

    # 筛选 3: 进化速率 (树长)
    if len(og_stats_filtered) > 10: # 只有足够多基因时才进行百分位筛选
        tree_lengths = [og['tree_len'] for og in og_stats_filtered]
        low_threshold = get_percentile(tree_lengths, LOW_RATE_PERCENTILE)
        high_threshold = get_percentile(tree_lengths, HIGH_RATE_PERCENTILE)
        print(f"  计算树长阈值: {low_threshold:.4f} ({LOW_RATE_PERCENTILE}%) - {high_threshold:.4f} ({HIGH_RATE_PERCENTILE}%)")
        
        og_stats_filtered = [og for og in og_stats_filtered if low_threshold <= og['tree_len'] <= high_threshold]
        print(f"进化速率筛选: {len(tree_lengths)} -> {len(og_stats_filtered)}")
    else:
        print("基因数量不足，跳过进化速率筛选。")

    final_ogs = og_stats_filtered
    print(f"最终保留 {len(final_ogs)} 个基因进行串联分析。")

    if not final_ogs:
        print("错误: 筛选后没有剩余基因。请检查阈值或输入文件。", file=sys.stderr)
        return

    # --- 3. 构建串联 MSA 和分区文件 ---
    print("\n--- 3. 正在构建串联 MSA 和分区文件 ---")
    concatenated_sequences = defaultdict(str)
    partitions = []
    current_msa_start_pos = 1

    for og_data in final_ogs:
        og_id = og_data['id']
        msa_path = og_data['msa']
        model = og_data['model']

        print(f"  串联 {og_id}...")
        sequences, msa_len, _ = read_fasta(msa_path) # 重新读取以获取序列
        if not sequences or msa_len == 0:
            print(f"  -> 严重警告: 无法重读 {og_id} 的 MSA，跳过！", file=sys.stderr)
            continue

        for sp_header in all_species_ordered:
            if sp_header in sequences:
                concatenated_sequences[sp_header] += sequences[sp_header]
            else:
                concatenated_sequences[sp_header] += "-" * msa_len

        partition_end_pos = current_msa_start_pos + msa_len - 1
        partitions.append(f"{model}, {og_id} = {current_msa_start_pos}-{partition_end_pos}")
        current_msa_start_pos = partition_end_pos + 1

    # --- 4. 写入输出文件 ---
    print(f"\n正在写入合并后的 MSA 到: {concatenated_msa_out} ...")
    with open(concatenated_msa_out, 'w') as f_out:
        for sp_header in all_species_ordered:
            f_out.write(f">{sp_header}\n")
            seq = concatenated_sequences[sp_header]
            for i in range(0, len(seq), 60):
                f_out.write(seq[i:i+60] + "\n")
    print("合并后的 MSA 文件写入完成。")

    print(f"\n正在写入分区文件到: {partition_file_out} ...")
    with open(partition_file_out, 'w') as f_out:
        for part_info in partitions:
            f_out.write(part_info + "\n")
    print("分区文件写入完成。")
    print(f"\n成功处理并串联了 {len(final_ogs)} 个 Orthogroups。")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("用法: python prepare_partitioned_analysis.py <AA_MSA目录> <IQ-TREE结果目录> <输出合并MSA文件名> <输出分区文件名>")
        print("示例: python prepare_partitioned_analysis.py SCOGs_msa_aa_clipkit gene_trees_with_aa_msa_MFP concatenated_filtered_aa.fasta partitions_filtered_aa.txt")
        sys.exit(1)

    msa_directory = sys.argv[1]
    iqtree_results_directory = sys.argv[2]
    concatenated_msa_output_file = sys.argv[3]
    partition_output_file = sys.argv[4]

    main(msa_directory, iqtree_results_directory, concatenated_msa_output_file, partition_output_file)