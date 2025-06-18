#!/bin/bash
#
# 2_run_iqtree_array.sh
#
# 功能:
# 作为一个 Slurm 作业数组运行 IQ-TREE。
# 它会读取由 1_prepare_resampling_list.sh 生成的文件列表。
#
# 使用方法:
# 1. 确保 1_prepare_resampling_list.sh 已成功运行。
# 2. 提交 Slurm 作业:
#    sbatch 2_run_iqtree_array.sh
#

# --- Slurm 配置 ---
#SBATCH --job-name=iqtree_resampling_codon
#SBATCH --output=slurm_logs/iqtree_resampling_codon_%A_%a.out
#SBATCH --error=slurm_logs/iqtree_resampling_codon_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

# --- 用户可配置变量 ---

# 输出 IQ-TREE 结果的目录
OUTPUT_DIR="./gene_trees_with_codon_MFP_resampling"
# 由第一个脚本生成的、去重后的文件列表
UNIQUE_FILE_LIST="unique_sampled_files_codon_resampling.txt"

# IQ-TREE 命令
IQTree_CMD="singularity exec /usr/local/biotools/i/iqtree:2.3.6--h503566f_1 iqtree"
# 每个 IQ-TREE 任务使用的 CPU 核心数
THREADS_PER_JOB=4


# --- IQ-TREE 运行函数 (由 Slurm 作业数组中的每个任务执行) ---
run_iqtree_task() {
    # 检查文件列表是否存在
    if [ ! -f "$UNIQUE_FILE_LIST" ]; then
        echo "错误: 未找到文件列表 '$UNIQUE_FILE_LIST'。"
        echo "请先运行脚本: bash 1_prepare_resampling_list.sh"
        exit 1
    fi

    # 获取此作业数组任务应该处理的文件
    # SLURM_ARRAY_TASK_ID 是 Slurm 提供的环境变量
    fasta_file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$UNIQUE_FILE_LIST")

    if [ -z "$fasta_file" ]; then
        echo "错误: 无法从 '$UNIQUE_FILE_LIST' 获取第 ${SLURM_ARRAY_TASK_ID} 行的文件名。任务退出。"
        exit 1
    fi

    # 创建输出和日志目录 (检查以确保它们存在)
    mkdir -p "$OUTPUT_DIR"
    mkdir -p slurm_logs
    
    # 获取文件名基本部分
    base_name=$(basename "$fasta_file" .fasta)

    echo "=================================================="
    echo "Slurm Job ID: ${SLURM_JOB_ID}, Array Task ID: ${SLURM_ARRAY_TASK_ID}"
    echo "开始处理文件: $fasta_file"
    echo "输出文件前缀: $OUTPUT_DIR/$base_name"
    echo "--------------------------------------------------"

    # 运行 IQ-TREE
    $IQTree_CMD -s "$fasta_file" -B 1000 -st CODON -m MFP -T $THREADS_PER_JOB --prefix "$OUTPUT_DIR/$base_name" --quiet

    echo "任务 ${SLURM_ARRAY_TASK_ID} 完成。"
    echo "=================================================="
}

# --- 脚本主逻辑 ---

# 如果 SLURM_ARRAY_TASK_ID 存在，说明这是一个作业数组中的任务
if [[ -n "$SLURM_ARRAY_TASK_ID" ]]; then
    run_iqtree_task
else
    # 如果 SLURM_ARRAY_TASK_ID 不存在，说明这是第一次提交
    # 需要计算数组大小并使用 --array 选项重新提交
    
    if [ ! -f "$UNIQUE_FILE_LIST" ]; then
        echo "错误: 文件列表 '$UNIQUE_FILE_LIST' 不存在。"
        echo "请首先运行脚本: bash 1_prepare_resampling_list.sh"
        exit 1
    fi
    
    # 计算需要运行的任务总数
    NUM_UNIQUE_FILES=$(wc -l < "$UNIQUE_FILE_LIST")
    if [ "$NUM_UNIQUE_FILES" -eq 0 ]; then
        echo "错误: 文件列表 '$UNIQUE_FILE_LIST' 为空，没有任务可运行。"
        exit 1
    fi
    
    # Slurm 数组索引从1开始
    echo "检测到 $NUM_UNIQUE_FILES 个独立文件。将提交一个范围为 1-$NUM_UNIQUE_FILES 的作业数组。"
    
    # 使用 sbatch 提交自身，并动态设置 --array 参数
    sbatch --array=1-${NUM_UNIQUE_FILES} "$0"
    
    echo "作业数组已提交。"
fi

# 脚本结束 