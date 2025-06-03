import os
import re
import argparse
from scipy.stats import chi2

def parse_lnL(file_path):
    try:
        with open(file_path, 'r') as f:
            for line in f:
                match = re.search(r'lnL\(ntime:\s*\d+\s+np:\s*(\d+)\):\s*(-?\d+\.\d+)', line)
                if match:
                    np = int(match.group(1))
                    lnL = float(match.group(2))
                    return np, lnL
        print(f"在文件 {file_path} 中未找到有效的 'lnL' 和 'np' 行。请检查文件格式。")
        return None, None
    except Exception as e:
        print(f"打开或解析文件 {file_path} 时发生错误: {e}")
        return None, None

def perform_lrt(alt_lnL, alt_np, null_lnL, null_np):
    try:
        if alt_lnL < null_lnL:
            print(f"警告: 备择模型的 lnL ({alt_lnL}) 小于零假设模型的 lnL ({null_lnL})。这可能表示模型设置或结果文件有问题。")
        
        df = alt_np - null_np
        if df <= 0:
            print(f"错误: 自由度 ({df}) 不是正数。备择模型的参数数量 ({alt_np}) 必须大于零假设模型的参数数量 ({null_np})。")
            return None

        lr_stat = 2 * (alt_lnL - null_lnL)
        if lr_stat < 0:
            print(f"警告: LR统计量 ({lr_stat}) 为负。将p值设为1。")
            return 1.0
            
        p_val = chi2.sf(lr_stat, df)
        return p_val
    except Exception as e:
        print(f"LRT计算过程中发生错误: {e}")
        return None

def main():
    parser = argparse.ArgumentParser(description="对PAML codeml的输出执行似然比检验 (LRT)。")
    parser.add_argument('--alt_dir', required=True, help='包含备择模型结果文件的目录路径。')
    parser.add_argument('--null_dir', required=True, help='包含零假设模型结果文件的目录路径。')
    parser.add_argument('--output_file', required=True, help='LRT结果的输出文件路径。')
    parser.add_argument('--alt_suffix', required=True, help="备择模型结果文件的后缀 (例如 '_alt_paml_results.txt')。")
    parser.add_argument('--null_suffix', required=True, help="零假设模型结果文件的后缀 (例如 '_null_paml_results.txt')。")
    
    args = parser.parse_args()

    alt_dir = os.path.expanduser(args.alt_dir)
    null_dir = os.path.expanduser(args.null_dir)
    output_file = os.path.expanduser(args.output_file)

    if not os.path.isdir(alt_dir):
        print(f"错误: 备择模型目录 '{alt_dir}' 不存在或不是一个目录。")
        return
    if not os.path.isdir(null_dir):
        print(f"错误: 零假设模型目录 '{null_dir}' 不存在或不是一个目录。")
        return

    try:
        alt_model_files = [f for f in os.listdir(alt_dir) if f.endswith(args.alt_suffix)]
    except FileNotFoundError:
        print(f"错误: 无法访问备择模型目录 '{alt_dir}'。")
        return
        
    if not alt_model_files:
        print(f"在目录 '{alt_dir}' 中没有找到后缀为 '{args.alt_suffix}' 的文件。")
        return

    with open(output_file, 'w') as out_f:
        out_f.write('Gene_ID\tp_val\tpositive_selection_detected\n')

        for idx, alt_filename in enumerate(alt_model_files):
            base_name = alt_filename[:-len(args.alt_suffix)] if alt_filename.endswith(args.alt_suffix) else alt_filename
            
            report_id = base_name.split('_')[0] if '_' in base_name else base_name

            alt_filepath = os.path.join(alt_dir, alt_filename)
            
            null_filename = base_name + args.null_suffix
            null_filepath = os.path.join(null_dir, null_filename)

            print(f"正在处理 {idx+1}/{len(alt_model_files)}: 基因ID '{report_id}' (备择文件: {alt_filename})...")

            if os.path.exists(null_filepath):
                alt_np, alt_lnL = parse_lnL(alt_filepath)
                null_np, null_lnL = parse_lnL(null_filepath)

                if alt_np is not None and alt_lnL is not None and null_np is not None and null_lnL is not None:
                    p_val = perform_lrt(alt_lnL, alt_np, null_lnL, null_np)
                    if p_val is not None:
                        reject_null = '+' if p_val < 0.05 else '-'
                        out_f.write(f'{report_id}\t{p_val:.6g}\t{reject_null}\n')
                        print(f"基因ID '{report_id}' 分析完成。P值: {p_val:.4g}, 正选择检测: {reject_null}")
                    else:
                        print(f"基因ID '{report_id}' 的LRT计算失败。")
                else:
                    print(f"基因ID '{report_id}' 的lnL或np数据不完整 (备择: np={alt_np}, lnL={alt_lnL}; 零假设: np={null_np}, lnL={null_lnL})。")
            else:
                print(f"未找到基因ID '{report_id}' 对应的零假设模型文件: {null_filepath}。")
        print(f"\n所有分析完成。结果已保存到 {output_file}")

if __name__ == "__main__":
    main()