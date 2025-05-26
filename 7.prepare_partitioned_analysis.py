import os
import re
import sys
from collections import defaultdict

# --- 用户提供的基本物种名列表 ---
BASE_SPECIES_NAMES = sorted([
    "Acromyrmex_echinatior", "Apis_cerana", "Apis_dorsata", "Apis_florea",
    "Apis_laboriosa", "Apis_mellifera", "Atta_cephalotes", "Atta_colombica",
    "Bombus_affinis", "Bombus_bifarius", "Bombus_fervidus", "Bombus_flavifrons",
    "Bombus_huntii", "Bombus_impatiens", "Bombus_pascuorum", "Bombus_pyrosoma",
    "Bombus_terrestris", "Bombus_vancouverensis_nearcticus", "Bombus_vosnesenskii",
    "Camponotus_floridanus", "Cardiocondyla_obscurior", "Cataglyphis_hispanica",
    "Ceratina_calcarata", "Colletes_gigas", "Cyphomyrmex_costatus",
    "Dinoponera_quadriceps", "Drosophila_melanogaster", "Dufourea_novaeangliae",
    "Eufriesea_mexicana", "Formica_exsecta", "Frieseomelitta_varia",
    "Habropoda_laboriosa", "Harpegnathos_saltator", "Hylaeus_anthracinus",
    "Hylaeus_volcanicus", "Linepithema_humile", "Megachile_rotundata",
    "Megalopta_genalis", "Monomorium_pharaonis", "Nomia_melanderi",
    "Nylanderia_fulva", "Odontomachus_brunneus", "Ooceraea_biroi",
    "Osmia_bicornis_bicornis", "Osmia_lignaria", "Pogonomyrmex_barbatus",
    "Polistes_canadensis", "Polistes_dominula", "Polistes_fuscatus",
    "Polyergus_mexicanus", "Prorops_nasuta", "Pseudomyrmex_gracilis",
    "Solenopsis_invicta", "Temnothorax_curvispinosus",
    "Temnothorax_longispinosus", "Temnothorax_nylanderi", "Trachymyrmex_cornetzi",
    "Trachymyrmex_septentrionalis", "Trachymyrmex_zeteki", "Vespa_crabro",
    "Vespa_mandarinia", "Vespa_velutina", "Vespula_pensylvanica",
    "Vespula_vulgaris", "Vollenhovia_emeryi", "Wasmannia_auropunctata"
], key=len, reverse=True) # 按长度降序排序，确保最长匹配优先

# --- 筛选阈值 ---
MIN_MSA_LEN = 200        # 最小比对长度 (AA)
MAX_BIAS_RATIO = 0.10    # 成分检验失败的最大比例 (例如 0.10 = 10%)
LOW_RATE_PERCENTILE = 10  # 要移除的最慢进化速率的百分比
HIGH_RATE_PERCENTILE = 90 # 要移除的最快进化速率的百分比
MIN_AVG_BOOTSTRAP = 80.0  # 内部节点平均自举支持率的最小阈值 (%)
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
                match = re.search(r"Total tree length \(sum of branch lengths\):\s*(\d+\.?\d*)", line)
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
    解析 IQ-TREE 日志文件 (.iqtree 或 .log) 以找到成分检验失败的序列数。
    总序列数会首先尝试从 .iqtree 的 'Input data: X sequences' 行获取，
    如果失败，则尝试从对应 .log 文件的 'Alignment has X sequences' 行获取。
    失败序列数会首先尝试从 .iqtree 的 'NUMBER OF SEQUENCES FAILED: X / Y' 模式获取，
    如果失败，则尝试从对应 .log 文件的 'N sequences failed composition chi2 test' 模式获取。
    返回一个元组 (failed_count, total_count) 或 (None, None)。
    """
    total_sequences = None
    failed_count = None

    # --- 步骤 1: 确定 .log 文件路径 ---
    log_file_path = None
    if iqtree_file_path.endswith(".clipkit.iqtree"):
        base_name = iqtree_file_path[:-len(".clipkit.iqtree")]
        log_file_path = base_name + ".clipkit.log"
    else:
        print(f"警告: IQ-TREE 文件名 '{iqtree_file_path}' 格式非预期，无法安全推断 .log 文件路径。", file=sys.stderr)
        # 即使无法推断 .log，仍然尝试仅从 .iqtree 文件获取信息

    # --- 步骤 2: 尝试获取总序列数 ---
    # 首先尝试从 .iqtree 文件
    try:
        with open(iqtree_file_path, 'r') as f_iqtree_total:
            for line in f_iqtree_total:
                match_total_iq = re.search(r"Input data:\s*(\d+)\s*sequences", line)
                if match_total_iq:
                    total_sequences = int(match_total_iq.group(1))
                    break
    except FileNotFoundError:
        print(f"警告: 主要 IQ-TREE 文件 (.iqtree) 未找到: {iqtree_file_path}", file=sys.stderr)
        # 不直接返回，允许尝试从 .log 获取总序列数
    except Exception as e:
        print(f"警告: 读取主要 IQ-TREE 文件 (.iqtree) '{iqtree_file_path}' 以获取总序列数时出错: {e}", file=sys.stderr)

    # 如果在 .iqtree 中未找到总序列数，并且 .log 文件路径有效，则尝试从 .log 文件获取
    if total_sequences is None and log_file_path and os.path.exists(log_file_path):
        try:
            with open(log_file_path, 'r') as f_log_total:
                for line in f_log_total:
                    match_total_log = re.search(r"Alignment has\s*(\d+)\s*sequences", line)
                    if match_total_log:
                        total_sequences = int(match_total_log.group(1))
                        print(f"信息: 从 {log_file_path} 获取到总序列数: {total_sequences}", file=sys.stderr)
                        break
        except Exception as e:
            print(f"警告: 读取 IQ-TREE .log 文件 '{log_file_path}' 以获取总序列数时出错: {e}", file=sys.stderr)

    if total_sequences is None:
        print(f"警告: 未能在 {iqtree_file_path} 或 (如果适用) {log_file_path or 'N/A'} 中找到总序列数信息。无法继续提取偏好统计。", file=sys.stderr)
        return None, None

    # --- 步骤 3: 尝试获取成分偏好失败计数 ---
    # 首先尝试从 .iqtree 文件中的主要模式
    try:
        with open(iqtree_file_path, 'r') as f_iqtree_bias:
            for line in f_iqtree_bias:
                if "NUMBER OF SEQUENCES FAILED:" in line:
                    match_bias_iq = re.search(r'(\d+)\s*/\s*(\d+)', line)
                    if match_bias_iq:
                        failed_count = int(match_bias_iq.group(1))
                        total_from_pattern = int(match_bias_iq.group(2))
                        if total_from_pattern != total_sequences:
                            print(f"警告: {iqtree_file_path} 中成分检验的总序列数 ({total_from_pattern}) 与先前确定的总序列数 ({total_sequences}) 不符。将使用模式中的值: {total_from_pattern}。", file=sys.stderr)
                            total_sequences = total_from_pattern # 更新总数
                        return failed_count, total_sequences
                    else:
                        print(f"警告: 在 {iqtree_file_path} 中找到了 'NUMBER OF SEQUENCES FAILED:' 行但无法解析: '{line.strip()}'。", file=sys.stderr)
                        return None, None # 主要模式行找到但解析失败，不再尝试备用
            # 如果循环完成，说明在 .iqtree 文件中未找到 "NUMBER OF SEQUENCES FAILED:" 行
    except FileNotFoundError:
        # 如果 .iqtree 文件在此阶段未找到 (理论上不太可能，因为前面已尝试读取总数)
        # 允许继续尝试 .log 文件获取偏好信息
        print(f"信息: {iqtree_file_path} 在查找偏好模式时未找到，将尝试 {log_file_path or 'N/A'}。", file=sys.stderr)
        pass # 确保会尝试log文件
    except Exception as e:
        print(f"警告: 解析 IQ-TREE iqtree 文件 '{iqtree_file_path}' 寻找主要偏好模式时出错: {e}", file=sys.stderr)
        # 即使出错，也尝试 .log 文件

    # 如果在 .iqtree 中未成功获取 failed_count (未找到行或文件读取出错)，并且 .log 文件路径有效，则尝试 .log 文件中的备用模式
    if failed_count is None and log_file_path and os.path.exists(log_file_path):
        try:
            with open(log_file_path, 'r') as f_log_bias:
                for line_log in f_log_bias:
                    match_log_failed = re.search(r"(\d+)\s*sequences failed composition chi2 test", line_log)
                    if match_log_failed:
                        failed_count = int(match_log_failed.group(1))
                        print(f"信息: 在 {log_file_path} 中找到备用成分偏好格式。失败数: {failed_count}, 总序列数: {total_sequences}", file=sys.stderr)
                        return failed_count, total_sequences
            # 如果 .log 文件中也未找到备用模式
            print(f"信息: 未在 {iqtree_file_path} 中找到主要偏好模式，在 {log_file_path} 中也未找到备用偏好模式。", file=sys.stderr)
            return None, None # 明确表示两个地方都没找到

        except Exception as e:
            print(f"警告: 解析 IQ-TREE .log 文件 '{log_file_path}' 寻找备用偏好模式时出错: {e}", file=sys.stderr)
            return None, None # .log 解析出错
    elif failed_count is None: # .log 文件路径无效或不存在，且 .iqtree 中未找到
        print(f"警告: 未能在 {iqtree_file_path} 中找到 'NUMBER OF SEQUENCES FAILED' 行，且对应的 .log 文件 ({log_file_path or 'N/A'}) 不适用或未找到备用模式。", file=sys.stderr)
        return None, None

    # 此处逻辑上应已被覆盖，若到达表示前面判断有遗漏或 failed_count 非 None 但未返回
    if failed_count is not None: # 意味着 failed_count 从.iqtree的主要模式中获取但由于某种原因未返回
        print(f"警告: 状态异常，failed_count 已设置 ({failed_count}) 但未提前返回。", file=sys.stderr)
        return failed_count, total_sequences
    
    print(f"警告: 未能从 {iqtree_file_path} 或 {log_file_path or 'N/A'} 提取成分偏好统计信息。", file=sys.stderr)
    return None, None # 最终的回退

def get_species_name_from_header(header, base_species_list):
    """
    将FASTA头部映射到预定义的基本物种名列表中的一个。
    优先匹配列表中更长（更具体）的物种名。
    如果找不到匹配，则打印警告并返回原始头部。
    """
    for species_name in base_species_list:
        if header.startswith(species_name):
            return species_name
    print(f"警告: 序列头部 '{header}' 未能映射到任何预定义的BASE_SPECIES_NAMES。将使用原始头部。", file=sys.stderr)
    # 可以在此处添加flush，如果日志句柄在作用域内的话
    # if log_f_handle and hasattr(log_f_handle, 'flush'): log_f_handle.flush()
    return header # 或返回 None 并让调用者处理

def read_fasta(fasta_file_path):
    """
    读取 FASTA 文件并返回 {species_name: sequence} 字典、
    第一个序列的长度以及按出现顺序排列的映射后物种名列表。
    如果多个头部映射到同一物种名，则后出现的序列会覆盖先出现的。
    """
    sequences = {}
    ordered_species_names = [] # 存储映射后的物种名
    seq_len = 0
    header_count = 0 # 用于跟踪原始头部数量
    mapped_species_in_file = set() # 跟踪此文件中已映射的物种名

    try:
        with open(fasta_file_path, 'r') as f:
            original_header = None
            current_seq_parts = []
            first_seq_processed = False

            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    header_count += 1
                    if original_header:
                        # 处理上一个序列
                        species_name = get_species_name_from_header(original_header, BASE_SPECIES_NAMES)
                        current_seq_str = "".join(current_seq_parts)
                        
                        if species_name in sequences and species_name in mapped_species_in_file:
                            # 如果在同一个文件中，同一个物种名已经从不同的原始头部映射过来并被赋值过
                            # 这意味着一个物种在这个OG的MSA中有多个序列条目
                            print(f"警告: 在文件 {fasta_file_path} 中，物种名 '{species_name}' (来自原始头部 '{original_header}') 多次出现。将使用最后出现的序列。原有长度: {len(sequences[species_name])}, 新长度: {len(current_seq_str)}", file=sys.stderr)
                            # if log_f_handle: log_f_handle.flush() # 如果log_f_handle在此作用域
                        
                        sequences[species_name] = current_seq_str
                        if not first_seq_processed:
                            seq_len = len(current_seq_str)
                            first_seq_processed = True
                        
                        if species_name not in mapped_species_in_file:
                            ordered_species_names.append(species_name)
                            mapped_species_in_file.add(species_name)
                        elif species_name not in ordered_species_names: # 物种已映射但未在ordered_species_names，理论不应发生
                             ordered_species_names.append(species_name) 

                    original_header = line[1:] # 新的原始头部
                    current_seq_parts = []
                else:
                    current_seq_parts.append(line)
            
            # 处理文件中的最后一个序列
            if original_header:
                species_name = get_species_name_from_header(original_header, BASE_SPECIES_NAMES)
                current_seq_str = "".join(current_seq_parts)

                if species_name in sequences and species_name in mapped_species_in_file:
                     print(f"警告: 在文件 {fasta_file_path} 中，物种名 '{species_name}' (来自原始头部 '{original_header}') 多次出现 (文件末尾)。将使用最后出现的序列。原有长度: {len(sequences[species_name])}, 新长度: {len(current_seq_str)}", file=sys.stderr)
                     # if log_f_handle: log_f_handle.flush()

                sequences[species_name] = current_seq_str
                if not first_seq_processed:
                    seq_len = len(current_seq_str)
                    # first_seq_processed = True # 不需要，因为这是文件末尾

                if species_name not in mapped_species_in_file:
                    ordered_species_names.append(species_name)
                    # mapped_species_in_file.add(species_name) # 不需要，这是文件末尾
                elif species_name not in ordered_species_names:
                     ordered_species_names.append(species_name) 

        if header_count == 0 and not sequences: # 文件可能是空的或非FASTA
            print(f"警告: FASTA 文件 {fasta_file_path} 为空或不含序列。", file=sys.stderr)
            # if log_f_handle: log_f_handle.flush()
            return None, 0, []

    except FileNotFoundError:
        print(f"警告: FASTA 文件未找到: {fasta_file_path}", file=sys.stderr)
        # if log_f_handle: log_f_handle.flush()
        return None, 0, []
    except Exception as e:
        print(f"警告: 读取 FASTA 文件 '{fasta_file_path}' 时出错: {e}", file=sys.stderr)
        # if log_f_handle: log_f_handle.flush()
        return None, 0, []
    
    return sequences, seq_len, ordered_species_names

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

def get_average_bootstrap_from_treefile(treefile_path):
    """
    解析 IQ-TREE .treefile 以计算内部节点的平均自举支持率。
    返回平均自举值 (浮点数)，如果发生错误或未找到支持值，则返回 None。
    对于序列数 <=3 的树，假定通过，返回 100.0。
    """
    try:
        with open(treefile_path, 'r') as f:
            tree_string = f.readline().strip() # 假设树在第一行

        if not tree_string:
            print(f"警告: Treefile {treefile_path} 为空。", file=sys.stderr)
            return None

        # 统计Tip数量，简单树特殊处理
        num_tips = tree_string.count(',') + 1 if ',' in tree_string else tree_string.count(':') # 逗号数+1 或 冒号数
        if not tree_string.strip(";").strip("()").strip(): # 处理空树或单节点树 "();" or "A;"
             num_tips = 0 
             if ':' in tree_string: # (A:0.1); is one tip
                 num_tips = 1

        if num_tips <= 3 and num_tips > 0: # 对于1, 2, 或 3 个序列的树
            print(f"信息: 在 {treefile_path} 中序列数 ({num_tips}) <=3。假定通过自举值筛选，平均值设为 100.0。", file=sys.stderr)
            return 100.0
        elif num_tips == 0 : # 完全空的树字符串
            print(f"警告: 在 {treefile_path} 中似乎是空树或无法识别的树结构。", file=sys.stderr)
            return None

        # 正则表达式查找自举值 (整数或浮点数)，它们通常跟在内部节点的右括号后
        bootstrap_values_str = re.findall(r"\)(\d+(?:\.\d+)?)[:,\]\)]", tree_string)
        # 更宽松的匹配，仅查找右括号后的数字，可能需要根据IQ-TREE版本调整，确保只匹配内部节点
        if not bootstrap_values_str:
             bootstrap_values_str = re.findall(r"\)(\d+(?:\.\d+)?)?", tree_string)
             # 过滤掉可能由 branch length 带来的空匹配或非数字
             bootstrap_values_str = [bs for bs in bootstrap_values_str if bs and bs!='.']


        if not bootstrap_values_str:
            print(f"警告: 在 {treefile_path} (序列数: {num_tips}) 中未找到内部节点自举值。", file=sys.stderr)
            return None 

        bootstrap_values_float = []
        for val_str in bootstrap_values_str:
            try:
                bootstrap_values_float.append(float(val_str))
            except ValueError:
                print(f"警告: 在 {treefile_path} 中无法将自举值 '{val_str}' 转换为浮点数。", file=sys.stderr)
                # 可选择跳过此值或使整个OG失败
        
        if not bootstrap_values_float:
            print(f"警告: 在 {treefile_path} 中所有提取的自举值都无法解析为数字。", file=sys.stderr)
            return None

        average_bootstrap = sum(bootstrap_values_float) / len(bootstrap_values_float)
        return average_bootstrap

    except FileNotFoundError:
        print(f"警告: Treefile 未找到: {treefile_path}", file=sys.stderr)
    except Exception as e:
        print(f"警告: 解析 treefile '{treefile_path}' 时出错 (Bootstrap): {e}", file=sys.stderr)
    return None

def main(msa_dir, iqtree_results_dir, concatenated_msa_out, partition_file_out):
    """
    主函数，用于收集数据、筛选基因并生成串联 MSA 和分区文件。
    """
    # --- 日志记录设置 ---
    base_log_name = os.path.splitext(concatenated_msa_out)[0]
    log_file_name = f"{base_log_name}_pipeline.log"
    
    # 初始消息打印到控制台
    print(f"脚本开始运行。所有详细输出将被重定向到日志文件: {os.path.abspath(log_file_name)}")

    original_stdout = sys.stdout
    original_stderr = sys.stderr
    log_f_handle = None  # 初始化日志文件句柄
    
    try:
        log_f_handle = open(log_file_name, 'w', encoding='utf-8')
        sys.stdout = log_f_handle
        sys.stderr = log_f_handle

        # --- 原 main 函数的核心逻辑开始 ---
        print(f"日志记录已启动。输出到: {log_file_name}")
        if log_f_handle: log_f_handle.flush()
        print(f"扫描 MSA 目录: {msa_dir}")
        if log_f_handle: log_f_handle.flush()
        print(f"扫描 IQ-TREE 结果目录: {iqtree_results_dir}")
        if log_f_handle: log_f_handle.flush()

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
        if log_f_handle: log_f_handle.flush()
        og_stats = []
        all_species_set = set()
        
        for og_id in sorted(og_files.keys()):
            og_data = og_files[og_id]
            if 'msa' not in og_data or 'iqtree_log' not in og_data:
                print(f"跳过 {og_id}: 缺少 MSA 或 IQ-TREE 日志。", file=sys.stderr)
                if log_f_handle: log_f_handle.flush()
                continue

            print(f"  正在处理 {og_id}...")
            if log_f_handle: log_f_handle.flush()
            sequences, msa_len, current_og_species_ordered_headers = read_fasta(og_data['msa'])
            if not sequences or msa_len == 0:
                print(f"  -> 跳过 {og_id}: MSA 为空或读取失败。", file=sys.stderr)
                if log_f_handle: log_f_handle.flush()
                continue

            model = get_best_model_from_iqtree_log(og_data['iqtree_log'])
            if not model:
                print(f"  -> 跳过 {og_id}: 无法提取模型。", file=sys.stderr)
                if log_f_handle: log_f_handle.flush()
                continue

            tree_len = get_tree_length_from_iqtree_log(og_data['iqtree_log'])
            if tree_len is None:
                print(f"  -> 跳过 {og_id}: 无法提取树长。", file=sys.stderr)
                if log_f_handle: log_f_handle.flush()
                continue

            failed, total = get_compositional_bias_stats(og_data['iqtree_log'])
            if failed is None: # total will also be None
                print(f"  -> 跳过 {og_id}: 无法提取成分偏好统计。", file=sys.stderr)
                if log_f_handle: log_f_handle.flush()
                continue
            
            bias_ratio = failed / total if total > 0 else 1.0

            # 新增：获取平均自举值
            treefile_name = f"{og_id}.clipkit.treefile"
            treefile_path = os.path.join(iqtree_results_dir, treefile_name)
            avg_bootstrap = None
            if os.path.exists(treefile_path):
                avg_bootstrap = get_average_bootstrap_from_treefile(treefile_path)
            else:
                print(f"警告: {og_id} 的 Treefile 未找到于: {treefile_path}。无法计算平均自举值。", file=sys.stderr)
                if log_f_handle: log_f_handle.flush()

            og_stats.append({
                'id': og_id, 'msa': og_data['msa'], 'iqtree': og_data['iqtree_log'],
                'len': msa_len, 'tree_len': tree_len, 'bias': bias_ratio,
                'model': model,
                'avg_bootstrap': avg_bootstrap # 添加平均自举值
            })
            all_species_set.update(sequences.keys())


        print(f"成功收集了 {len(og_stats)} 个基因的数据。")
        if log_f_handle: log_f_handle.flush()
        all_species_ordered = sorted(list(all_species_set)) 
        print(f"共发现 {len(all_species_ordered)} 个物种。")
        if log_f_handle: log_f_handle.flush()

        # --- 2. 基因筛选 ---
        print("\n--- 2. 正在筛选基因 ---")
        if log_f_handle: log_f_handle.flush()
        count_before_filtering = len(og_stats) 
        
        # 筛选 1: 比对长度
        og_stats_filtered = [og for og in og_stats if og['len'] >= MIN_MSA_LEN]
        print(f"比对长度筛选 ({MIN_MSA_LEN} AA): {count_before_filtering} -> {len(og_stats_filtered)}")
        if log_f_handle: log_f_handle.flush()
        count_before_filtering = len(og_stats_filtered)

        # 筛选 2: 成分偏好
        og_stats_filtered = [og for og in og_stats_filtered if og['bias'] <= MAX_BIAS_RATIO]
        print(f"成分偏好筛选 (<={MAX_BIAS_RATIO*100:.0f}%): {count_before_filtering} -> {len(og_stats_filtered)}")
        if log_f_handle: log_f_handle.flush()
        count_before_filtering = len(og_stats_filtered)

        # 筛选 3: 进化速率 (树长)
        if len(og_stats_filtered) > 10: 
            tree_lengths = [og['tree_len'] for og in og_stats_filtered]
            low_threshold = get_percentile(tree_lengths, LOW_RATE_PERCENTILE)
            high_threshold = get_percentile(tree_lengths, HIGH_RATE_PERCENTILE)
            
            if low_threshold is not None and high_threshold is not None:
                print(f"  计算树长阈值: {low_threshold:.4f} ({LOW_RATE_PERCENTILE}%) - {high_threshold:.4f} ({HIGH_RATE_PERCENTILE}%)")
                if log_f_handle: log_f_handle.flush()
                
                og_stats_temp_rate_filter = [og for og in og_stats_filtered if low_threshold <= og['tree_len'] <= high_threshold]
                print(f"进化速率筛选: {count_before_filtering} -> {len(og_stats_temp_rate_filter)}")
                if log_f_handle: log_f_handle.flush()
                og_stats_filtered = og_stats_temp_rate_filter
            else:
                print(f"  警告: 无法计算树长阈值 (数据不足或 get_percentile 返回 None)，跳过进化速率筛选。")
                if log_f_handle: log_f_handle.flush()
        else:
            print(f"基因数量不足 (当前 {len(og_stats_filtered)} 个，需 >10)，跳过进化速率筛选。")
            if log_f_handle: log_f_handle.flush()

        # 新增：筛选 4: 平均自举值
        if MIN_AVG_BOOTSTRAP > 0: # 只有当阈值大于0时才执行此筛选
            og_stats_filtered_bootstrap = []
            for og_entry in og_stats_filtered:
                avg_bs = og_entry.get('avg_bootstrap')
                if avg_bs is not None and avg_bs >= MIN_AVG_BOOTSTRAP:
                    og_stats_filtered_bootstrap.append(og_entry)
                elif avg_bs is None:
                    print(f"信息: OG {og_entry['id']} 因无法计算平均自举值而被自举值筛选移除。", file=sys.stderr)
                    if log_f_handle: log_f_handle.flush()
                else: # avg_bs is not None but < MIN_AVG_BOOTSTRAP
                    print(f"信息: OG {og_entry['id']} 平均自举值 {avg_bs:.2f}% < {MIN_AVG_BOOTSTRAP}%，被自举值筛选移除。", file=sys.stderr)
                    if log_f_handle: log_f_handle.flush()
            
            print(f"平均自举值筛选 (>= {MIN_AVG_BOOTSTRAP}%): {count_before_filtering} -> {len(og_stats_filtered_bootstrap)}")
            if log_f_handle: log_f_handle.flush()
            og_stats_filtered = og_stats_filtered_bootstrap
            count_before_filtering = len(og_stats_filtered) # 更新计数器，以防未来有更多筛选步骤

        final_ogs = og_stats_filtered
        print(f"最终保留 {len(final_ogs)} 个基因进行串联分析。")
        if log_f_handle: log_f_handle.flush()

        if not final_ogs:
            print("错误: 筛选后没有剩余基因。请检查阈值或输入文件。", file=sys.stderr)
            if log_f_handle: log_f_handle.flush()
            return # 允许 finally 块执行

        # --- 3. 构建串联 MSA 和分区文件 ---
        print("\n--- 3. 正在构建串联 MSA 和分区文件 ---")
        if log_f_handle: log_f_handle.flush()
        concatenated_sequences = defaultdict(str)
        partitions = []
        current_msa_start_pos = 1

        for og_data in final_ogs:
            og_id = og_data['id']
            msa_path = og_data['msa']
            model = og_data['model']

            print(f"  串联 {og_id}...")
            if log_f_handle: log_f_handle.flush()
            sequences, msa_len, _ = read_fasta(msa_path) 
            if not sequences or msa_len == 0:
                print(f"  -> 严重警告: 无法重读 {og_id} 的 MSA，跳过！", file=sys.stderr)
                if log_f_handle: log_f_handle.flush()
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
        if log_f_handle: log_f_handle.flush()
        with open(concatenated_msa_out, 'w') as f_out:
            for sp_header in all_species_ordered:
                f_out.write(f">{sp_header}\n")
                seq = concatenated_sequences[sp_header]
                for i in range(0, len(seq), 60):
                    f_out.write(seq[i:i+60] + "\n")
        print("合并后的 MSA 文件写入完成。")
        if log_f_handle: log_f_handle.flush()

        print(f"\n正在写入分区文件到: {partition_file_out} ...")
        if log_f_handle: log_f_handle.flush()
        with open(partition_file_out, 'w') as f_out:
            for part_info in partitions:
                f_out.write(part_info + "\n")
        print("分区文件写入完成。")
        if log_f_handle: log_f_handle.flush()
        print(f"\n成功处理并串联了 {len(final_ogs)} 个 Orthogroups。")
        # --- 原 main 函数的核心逻辑结束 ---

    finally:
        # 恢复 stdout 和 stderr
        if sys.stdout != original_stdout and log_f_handle: # 检查是否真的发生了重定向
            sys.stdout = original_stdout
        if sys.stderr != original_stderr and log_f_handle: # 检查是否真的发生了重定向
            sys.stderr = original_stderr
        
        if log_f_handle: # 关闭日志文件句柄
            log_f_handle.close()
        
        # 最终消息打印到控制台
        print(f"脚本执行完毕。日志文件已保存到: {os.path.abspath(log_file_name)}")


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