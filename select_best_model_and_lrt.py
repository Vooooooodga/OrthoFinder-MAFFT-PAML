import os
import re
import argparse
import collections
from scipy.stats import chi2

def parse_lnL(file_path):
    """从PAML结果文件中解析 np 和 lnL 值。"""
    try:
        with open(file_path, 'r') as f:
            for line in f:
                # 适配 PAML 输出格式: lnL(ntime: X np: Y): ZZZZ.ZZZZZ
                # 或者有时是 Tree#1: lnL(ntime: X np: Y): ZZZZ.ZZZZZ
                match = re.search(r'lnL\(ntime:\s*\d+\s+np:\s*(\d+)\):\s*(-?\d+\.\d+)', line)
                if match:
                    np = int(match.group(1))
                    lnL = float(match.group(2))
                    return np, lnL
        print(f"警告: 在文件 {file_path} 中未找到有效的 'lnL' 和 'np' 行。")
        return None, None
    except FileNotFoundError:
        print(f"错误: 文件 {file_path} 未找到。")
        return None, None
    except Exception as e:
        print(f"错误: 打开或解析文件 {file_path} 时发生错误: {e}")
        return None, None

def parse_branch_w_values(file_path):
    """从PAML branch model结果文件中解析 w (dN/dS) for branches 值。"""
    try:
        with open(file_path, 'r') as f:
            for line in f:
                if "w (dN/dS) for branches:" in line:
                    w_values_str = line.split(":", 1)[1].strip()
                    return w_values_str
        return None
    except FileNotFoundError:
        return None
    except Exception:
        return None

def parse_bsa_beb_sites(file_path):
    """从PAML branch-site model结果文件中解析BEB分析中带星号的位点信息。"""
    beb_sites = []
    in_beb_section = False
    try:
        with open(file_path, 'r') as f:
            for line in f:
                if "Bayes Empirical Bayes (BEB) analysis" in line:
                    in_beb_section = True
                    continue

                if in_beb_section:
                    stripped_line = line.strip()
                    if not stripped_line or "The grid" in stripped_line or "Prob(w>1) for branches" in stripped_line:
                        in_beb_section = False
                        break
                    
                    if "*" in stripped_line:
                        # 简单的启发式检查，是否像一个位点行
                        parts = stripped_line.split()
                        if len(parts) >= 3 and parts[0].isdigit() and parts[1].isalpha() and len(parts[1]) == 1:
                            beb_sites.append(stripped_line)
        
        if beb_sites:
            return "; ".join(beb_sites)
        return None
    except FileNotFoundError:
        return None
    except Exception:
        return None

def select_best_model_from_list(file_paths, gene_id, model_type_description):
    """
    从给定的文件路径列表中选择具有最高 lnL 值的模型。
    打印选择的详情到屏幕。
    """
    best_lnL = -float('inf')
    best_np = None
    best_file_path = None

    if not file_paths:
        print(f"信息: 基因 {gene_id} ({model_type_description}): 没有找到候选模型文件。")
        return None, None, None

    for f_path in file_paths:
        np, lnL = parse_lnL(f_path)
        if np is not None and lnL is not None:
            if lnL > best_lnL:
                best_lnL = lnL
                best_np = np
                best_file_path = f_path
    
    if best_file_path:
        print(f"基因 {gene_id} ({model_type_description}): "
              f"从 {len(file_paths)} 个候选中选择 {os.path.basename(best_file_path)} "
              f"(lnL: {best_lnL:.4f}, np: {best_np})")
    else:
        print(f"警告: 基因 {gene_id} ({model_type_description}): "
              f"无法从 {len(file_paths)} 个候选文件中选出最佳模型 (可能所有文件都无法解析)。")

    return best_file_path, best_lnL, best_np

def perform_lrt(alt_lnL, alt_np, null_lnL, null_np):
    """
    执行似然比检验。
    返回 (LRT统计量, 自由度, p值)。
    """
    try:
        if alt_lnL is None or alt_np is None or null_lnL is None or null_np is None:
            print("错误: LRT的输入参数 (lnL, np) 存在None值。")
            return None, None, None

        # PAML 有时会因为数值问题，导致备择模型lnL略小于零假设模型，即使它更复杂。
        # 传统LRT要求备择模型lnL >= 零假设模型lnL。
        # if alt_lnL < null_lnL:
        #     print(f"警告: 备择模型的 lnL ({alt_lnL:.4f}) 小于零假设模型的 lnL ({null_lnL:.4f})。")
            # 根据上下文，这可能是一个问题，或者只是数值噪音。
            # 对于某些比较 (如 M1a vs M2a)，这可能意味着没有证据支持更复杂的模型。
            # 我们仍然可以计算，但p值可能不那么有意义，或者应该被保守地解释。

        df = alt_np - null_np
        if df <= 0:
            print(f"错误: 自由度 ({df}) 不是正数。备择模型的参数数量 ({alt_np}) "
                  f"必须大于零假设模型的参数数量 ({null_np})。")
            return None, df, None

        # LRT statistic: 2 * (lnL_alt - lnL_null)
        lr_stat = 2 * (alt_lnL - null_lnL)
        
        p_val = None
        if lr_stat < 0:
            # 如果 lr_stat < 0 且 df > 0, 这通常意味着备择模型拟合更差。
            # 在这种情况下，我们不能拒绝零假设，P值通常设为1。
            print(f"警告: LR统计量 ({lr_stat:.4f}) 为负 (alt_lnL: {alt_lnL:.4f}, null_lnL: {null_lnL:.4f}, alt_np: {alt_np}, null_np: {null_np})。"
                  " 这表示备择模型拟合不如零假设模型。将P值设为1。")
            p_val = 1.0
        else:
            p_val = chi2.sf(lr_stat, df) # Survival function (1 - CDF)

        return lr_stat, df, p_val
    except Exception as e:
        print(f"错误: LRT计算过程中发生错误: {e}")
        return None, None, None

def main():
    parser = argparse.ArgumentParser(
        description="对PAML codeml的输出进行处理：为每个基因选择最佳的多omega模型（branch和bsA_alt），"
                    "然后执行LRT并报告结果。",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('--paml_results_dir', required=True, 
                        help='包含所有PAML codeml输出文件的目录。')
    parser.add_argument('--output_file', required=True, 
                        help='LRT结果的CSV输出文件路径。')
    parser.add_argument('--omega_values', type=str, default="0.5,1.0,1.5,2.0",
                        help='逗号分隔的omega初始值，用于构建备择模型文件名 (例如 "0.5,1.0")。')

    args = parser.parse_args()

    paml_results_dir = os.path.expanduser(args.paml_results_dir)
    output_file = os.path.expanduser(args.output_file)

    if not os.path.isdir(paml_results_dir):
        print(f"错误: PAML结果目录 '{paml_results_dir}' 不存在或不是一个目录。")
        return

    omega_filename_suffixes = [val.replace('.', 'p') for val in args.omega_values.split(',')]

    # 正则表达式来捕获基因ID和模型类型指示符
    # 例如: OG0000001_codon.clipkit_branch_omega0p5_results.txt
    # gene_id_part = OG0000001_codon.clipkit
    # model_indicator = branch_omega0p5
    filename_pattern = re.compile(
        r"^(.*?)_((M0)|(branch_omega(?:{}))|(bsA_alt_omega(?:{}))|(bsA_null))_results\.txt$".format(
            "|".join(omega_filename_suffixes), "|".join(omega_filename_suffixes)
        )
    )
    
    # 组织文件
    gene_model_files = collections.defaultdict(lambda: {
        "M0": [], "branch_alt_candidates": [], "bsA_null": [], "bsA_alt_candidates": []
    })

    print(f"正在扫描目录: {paml_results_dir}")
    found_files_count = 0
    for fname in os.listdir(paml_results_dir):
        match = filename_pattern.match(fname)
        if match:
            found_files_count +=1
            gene_id_part = match.group(1)
            full_model_indicator = match.group(2)
            
            file_path = os.path.join(paml_results_dir, fname)

            if full_model_indicator == "M0":
                gene_model_files[gene_id_part]["M0"].append(file_path)
            elif full_model_indicator.startswith("branch_omega"):
                gene_model_files[gene_id_part]["branch_alt_candidates"].append(file_path)
            elif full_model_indicator == "bsA_null":
                gene_model_files[gene_id_part]["bsA_null"].append(file_path)
            elif full_model_indicator.startswith("bsA_alt_omega"):
                gene_model_files[gene_id_part]["bsA_alt_candidates"].append(file_path)
    
    print(f"扫描完成。从文件名中识别出 {len(gene_model_files)} 个独立基因ID，共匹配 {found_files_count} 个文件。\n")

    with open(output_file, 'w') as out_f:
        out_f.write("Gene_ID,Test_Type,Best_Alt_Model_File,Alt_lnL,Alt_np,"
                    "Null_Model_File,Null_lnL,Null_np,"
                    "LRT_Statistic,df,P_Value,Significance_0.05,"
                    "Branch_w_values,BEB_Positive_Sites\n")

        for gene_id, models in sorted(gene_model_files.items()):
            print(f"--- 正在处理基因: {gene_id} ---")
            
            # 1. Branch vs M0
            m0_files = models["M0"]
            branch_alt_candidate_files = models["branch_alt_candidates"]

            if not m0_files:
                print(f"警告: 基因 {gene_id}: 未找到M0模型文件。跳过 Branch vs M0 检验。")
            elif len(m0_files) > 1:
                print(f"警告: 基因 {gene_id}: 找到多个M0模型文件: {m0_files}。使用第一个。")
            
            if m0_files:
                m0_file_path = m0_files[0]
                m0_np, m0_lnL = parse_lnL(m0_file_path)

                if not branch_alt_candidate_files:
                    print(f"信息: 基因 {gene_id}: 未找到branch alternative候选模型文件。")
                else:
                    best_branch_alt_file, best_branch_alt_lnL, best_branch_alt_np = \
                        select_best_model_from_list(branch_alt_candidate_files, gene_id, "Branch Alternative")

                    if best_branch_alt_file and m0_np is not None and m0_lnL is not None:
                        lr_stat, df, p_val = perform_lrt(best_branch_alt_lnL, best_branch_alt_np, m0_lnL, m0_np)
                        if p_val is not None:
                            significance = '+' if p_val < 0.05 else '-'
                            branch_w_values_str = parse_branch_w_values(best_branch_alt_file) if best_branch_alt_file else ""
                            branch_w_values_csv = branch_w_values_str.replace('"', '""') if branch_w_values_str else ""
                            
                            out_f.write(
                                f"{gene_id},Branch_vs_M0,"
                                f"{os.path.basename(best_branch_alt_file)},{best_branch_alt_lnL or ''},{best_branch_alt_np or ''},"
                                f"{os.path.basename(m0_file_path)},{m0_lnL or ''},{m0_np or ''},"
                                f"{lr_stat if lr_stat is not None else ''},{df if df is not None else ''},{p_val:.6g},{significance},"
                                f"\"{branch_w_values_csv}\",\"\"\n"
                            )
                            print(f"基因 {gene_id} (Branch vs M0): LRT P-value = {p_val:.4g} ({significance})")
                        else:
                            print(f"基因 {gene_id} (Branch vs M0): LRT计算失败。")
                    elif not best_branch_alt_file:
                         print(f"信息: 基因 {gene_id}: 未能从候选者中选择最佳Branch Alternative模型。")
                    elif m0_np is None or m0_lnL is None:
                        print(f"信息: 基因 {gene_id}: M0模型 ({m0_file_path}) 解析失败。")


            # 2. Branch-Site A (Alt vs Null)
            bsa_null_files = models["bsA_null"]
            bsa_alt_candidate_files = models["bsA_alt_candidates"]

            if not bsa_null_files:
                print(f"警告: 基因 {gene_id}: 未找到Branch-Site A Null模型文件。跳过 BsA Alt vs Null 检验。")
            elif len(bsa_null_files) > 1:
                 print(f"警告: 基因 {gene_id}: 找到多个 BsA Null 模型文件: {bsa_null_files}。使用第一个。")

            if bsa_null_files:
                bsa_null_file_path = bsa_null_files[0]
                bsa_null_np, bsa_null_lnL = parse_lnL(bsa_null_file_path)
                
                if not bsa_alt_candidate_files:
                    print(f"信息: 基因 {gene_id}: 未找到Branch-Site A alternative候选模型文件。")
                else:
                    best_bsa_alt_file, best_bsa_alt_lnL, best_bsa_alt_np = \
                        select_best_model_from_list(bsa_alt_candidate_files, gene_id, "Branch-Site A Alternative")

                    if best_bsa_alt_file and bsa_null_np is not None and bsa_null_lnL is not None:
                        lr_stat, df, p_val = perform_lrt(best_bsa_alt_lnL, best_bsa_alt_np, bsa_null_lnL, bsa_null_np)
                        if p_val is not None:
                            significance = '+' if p_val < 0.05 else '-'
                            beb_sites_str = parse_bsa_beb_sites(best_bsa_alt_file) if best_bsa_alt_file else ""
                            beb_sites_csv = beb_sites_str.replace('"', '""') if beb_sites_str else ""

                            out_f.write(
                                f"{gene_id},BsA_Alt_vs_Null,"
                                f"{os.path.basename(best_bsa_alt_file)},{best_bsa_alt_lnL or ''},{best_bsa_alt_np or ''},"
                                f"{os.path.basename(bsa_null_file_path)},{bsa_null_lnL or ''},{bsa_null_np or ''},"
                                f"{lr_stat if lr_stat is not None else ''},{df if df is not None else ''},{p_val:.6g},{significance},"
                                f"\"\",\"{beb_sites_csv}\"\n"
                            )
                            print(f"基因 {gene_id} (BsA Alt vs Null): LRT P-value = {p_val:.4g} ({significance})")
                        else:
                            print(f"基因 {gene_id} (BsA Alt vs Null): LRT计算失败。")
                    elif not best_bsa_alt_file:
                        print(f"信息: 基因 {gene_id}: 未能从候选者中选择最佳BsA Alternative模型。")
                    elif bsa_null_np is None or bsa_null_lnL is None:
                        print(f"信息: 基因 {gene_id}: BsA Null模型 ({bsa_null_file_path}) 解析失败。")
            print("") # 为每个基因的输出添加空行

    print(f"\n所有分析完成。结果已保存到 {output_file}")

if __name__ == "__main__":
    main()