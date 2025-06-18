#!/bin/bash
#
# run_iqtree_array.sh
#
# 功能:
# 作为一个 Slurm 作业数组运行 IQ-TREE。
# 它会处理 aligned_codon_clipkit_70_coverage 目录下的所有 .fasta 文件。
#
# 使用方法:
# 1. 确保输入目录 aligned_codon_clipkit_70_coverage 存在且包含 FASTA 文件。
# 2. 提交 Slurm 作业:
#    sbatch run_iqtree_for_RELAX_Chi2_test.sh
#

# --- Slurm 配置 ---
#SBATCH --job-name=iqtree_array_codon
#SBATCH --output=slurm_logs/iqtree_codon_%A_%a.out
#SBATCH --error=slurm_logs/iqtree_codon_%A_%a.err
#SBATCH --array=1-8378
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

# --- 用户可配置变量 ---

# 输入 FASTA 文件的目录
INPUT_DIR="aligned_codon_clipkit_70_coverage"
# 输出 IQ-TREE 结果的目录
OUTPUT_DIR="./gene_trees_with_codon_MFP"

# IQ-TREE 命令
IQTree_CMD="singularity exec /usr/local/biotools/i/iqtree:2.3.6--h503566f_1 iqtree"
# 每个 IQ-TREE 任务使用的 CPU 核心数
THREADS_PER_JOB=4


# --- 脚本主逻辑 (由 Slurm 作业数组中的每个任务执行) ---

# 检查输入目录是否存在
if [ ! -d "$INPUT_DIR" ]; then
    echo "错误: 输入目录 '$INPUT_DIR' 不存在。"
    exit 1
fi

# 创建输出和日志目录 (如果不存在)
mkdir -p "$OUTPUT_DIR"
mkdir -p slurm_logs

# 获取此作业数组任务应该处理的文件
# SLURM_ARRAY_TASK_ID 是 Slurm 提供的环境变量 (1-based)
# 我们创建一个包含所有目标文件的bash数组 (0-based)
file_list=(${INPUT_DIR}/OG*_codon.clipkit.fasta)
fasta_file=${file_list[$SLURM_ARRAY_TASK_ID-1]}

if [ -z "$fasta_file" ] || [ ! -f "$fasta_file" ]; then
    echo "错误: 无法找到任务 ${SLURM_ARRAY_TASK_ID} 对应的文件。计算出的文件路径是 '$fasta_file'。任务退出。"
    exit 1
fi

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

# 脚本结束 