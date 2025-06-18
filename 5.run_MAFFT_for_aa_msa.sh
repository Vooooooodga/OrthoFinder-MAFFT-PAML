#!/bin/bash
set -e # 任何命令失败时立即退出
set -o pipefail # 管道中的任何命令失败都算作失败

#SBATCH --job-name=mafft-array
#SBATCH --output=logs/mafft_array_%A_%a.out
#SBATCH --error=logs/mafft_array_%A_%a.err
#SBATCH --array=1-8378 # 请根据 translated_proteins 目录中实际的 .fa 文件数量调整此处的数字
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=4 # 每个子任务需要4个核心
#SBATCH --mem=4G          # 每个子任务需要4G内存

ulimit -s unlimited

# --- 关键检查 ---
# 获取脚本所在的目录，以确保我们总能找到位于同一目录下的文件列表
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
file_list="${SCRIPT_DIR}/mafft_input_files.list"

# 在脚本早期就检查文件列表是否存在，如果不存在则快速失败并给出提示
if [ ! -f "$file_list" ]; then
    echo "错误：在脚本目录中找不到输入文件列表 '$file_list'！" >&2
    echo "请确保 'mafft_input_files.list' 文件与您的 sbatch 脚本位于同一目录。" >&2
    echo "您可以在该目录中运行以下命令来创建它：" >&2
    echo 'find "$(pwd)/translated_proteins" -type f -name "*.fa" > mafft_input_files.list' >&2
    exit 1
fi

echo "Running on $(hostname)"
echo "Starting at $(date)"

# 输入和输出目录定义
output="./aligned_translated_proteins"
CORES_PER_MAFFT_JOB=4 # 与 --cpus-per-task 保持一致

# 确保输出目录存在
mkdir -p $output
mkdir -p logs

# 从预先生成的文件列表中读取要处理的文件
total_files=$(wc -l < "$file_list")

if [ "$total_files" -eq 0 ]; then
    echo "Input file list '$file_list' is empty or not found. Exiting."
    exit 1
fi

# SLURM_ARRAY_TASK_ID 从 1 开始
task_id=${SLURM_ARRAY_TASK_ID}

# 检查任务ID是否在文件列表范围内
if [ -z "$task_id" ]; then
    echo "错误: SLURM_ARRAY_TASK_ID 变量为空。请确认您是通过 'sbatch' 命令提交的作业数组。" >&2
    exit 1
elif [ "$task_id" -gt "$total_files" ]; then
    echo "Task ID $task_id is out of bounds. Total files: $total_files. Exiting."
    exit 0
fi

# 使用 sed 命令从文件列表中获取当前任务对应的文件路径
file_to_process=$(sed -n "${task_id}p" "$file_list")
base_name=$(basename "${file_to_process%.fa}")

# --- 日志重定向 ---
# 将此任务的 stdout 和 stderr 重定向到与输入文件同名的日志文件中
exec > "logs/${base_name}.out" 2> "logs/${base_name}.err"

echo "SLURM Job ID: ${SLURM_JOB_ID}"
echo "SLURM Array Job ID: ${SLURM_ARRAY_JOB_ID}"
echo "SLURM Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "----------------------------------------------------"
echo "Running on $(hostname)"
echo "Starting at $(date)"

echo "Processing file $task_id of $total_files: $file_to_process"
outfile="$output/${base_name}_aligned.fa"
echo "Output will be saved to $outfile"

# 运行MAFFT
singularity exec /usr/local/biotools/m/mafft:7.525--h031d066_0 mafft --auto --thread $CORES_PER_MAFFT_JOB "$file_to_process" > "$outfile"

echo "MAFFT job for $file_to_process completed."
echo "Ending at $(date)"