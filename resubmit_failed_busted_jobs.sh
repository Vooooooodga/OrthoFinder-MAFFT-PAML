#!/bin/bash

# --- 配置 ---
# 主输出目录 (与 run_hyphy_busted.sh 中的 MAIN_OUTPUT_DIR 一致)
MAIN_OUTPUT_DIR="/home/yuhangjia/data/AlternativeSplicing/evo_rate_test_RNA_splicing_term/BUSTED"

# 存放原始子 sbatch 脚本的目录
SUB_SBATCH_SCRIPT_DIR="${MAIN_OUTPUT_DIR}/sub_sbatch_scripts"
# 存放原始子作业 SLURM 日志的目录
SLURM_LOGS_DIR="${MAIN_OUTPUT_DIR}/BUSTED_slurm_logs"

# 重新提交作业的新配置 (分区将保持原始设置)
# NEW_PARTITION="your_new_partition" # 不再修改分区
NEW_CPUS_PER_TASK=32
NEW_MEM="32G"
NEW_HYPHY_CPU_COUNT=32 # hyphy 命令内部的 CPU 参数

# 用于存放修改后用于重新提交的脚本的临时目录
RESUBMIT_SCRIPT_DIR="${MAIN_OUTPUT_DIR}/resubmit_scripts"
mkdir -p "${RESUBMIT_SCRIPT_DIR}" # 创建目录

echo "Starting to check for failed/incomplete jobs and resubmit them (will overwrite original logs, partition unchanged)..."

# 遍历 SLURM_LOGS_DIR 中的所有 .err 文件
find "${SLURM_LOGS_DIR}" -name "*_busted.err" -type f -print0 | while IFS= read -r -d $'\\0' err_file; do
    # 检查 .err 文件是否非空
    if [ -s "${err_file}" ]; then
        echo "Found non-empty .err file: ${err_file}"
        
        # 从 .err 文件名提取基因ID (假设格式为 GENEID_busted.err)
        base_err_name=$(basename "${err_file}")
        gene_id=$(echo "${base_err_name}" | sed 's/_busted\.err$//')
        
        echo "Gene ID identified as: ${gene_id}"
        
        # 构建原始子 sbatch 脚本的路径
        original_sub_script_path="${SUB_SBATCH_SCRIPT_DIR}/sub_busted_${gene_id}.sh"
        
        if [ ! -f "${original_sub_script_path}" ]; then
            echo "ERROR: Original sub-sbatch script not found for gene ${gene_id} at ${original_sub_script_path}. Skipping."
            continue
        fi

        echo "Original sub-script found: ${original_sub_script_path}"
        
        # 定义新的重提交脚本路径 (仍保存到 RESUBMIT_SCRIPT_DIR 以便追踪，但其内容会指向原始日志路径)
        resubmit_script_path="${RESUBMIT_SCRIPT_DIR}/resubmit_busted_${gene_id}.sh"
        
        # 修改原始子脚本内容以用于重新提交
        echo "Generating resubmission script (partition unchanged, to overwrite original logs): ${resubmit_script_path}"
        sed \
            -e "s/^#SBATCH --cpus-per-task=1/#SBATCH --cpus-per-task=${NEW_CPUS_PER_TASK}/" \
            -e "s/^#SBATCH --mem=4G/#SBATCH --mem=${NEW_MEM}/" \
            -e "s/hyphy CPU=1 busted/hyphy CPU=${NEW_HYPHY_CPU_COUNT} busted/" \
            "${original_sub_script_path}" > "${resubmit_script_path}" 
            # 注意移除了修改 partition, job-name, .out, .err 文件名的行

        if [ ! -s "${resubmit_script_path}" ]; then
            echo "ERROR: Failed to generate resubmission script for ${gene_id}. Skipping."
            continue
        fi
        
        chmod +x "${resubmit_script_path}"
        echo "Submitting resubmission script for ${gene_id} (partition unchanged, will overwrite original logs): ${resubmit_script_path}"
        sbatch "${resubmit_script_path}"
        echo "-----------------------------------------------------"
    else
        # echo "DEBUG: .err file is empty (or does not exist): ${err_file}" # 可选的调试信息
        true # 占位符，如果.err文件为空则什么也不做
    fi
done

echo "Resubmission process completed." 