#!/bin/bash
#SBATCH --job-name=mafft-array
#SBATCH --output=logs/mafft_array_%A_%a.out
#SBATCH --error=logs/mafft_array_%A_%a.err
#SBATCH --array=1-8378 # 请根据 translated_proteins 目录中实际的 .fa 文件数量调整此处的数字
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=4 # 每个子任务需要4个核心
#SBATCH --mem=4G          # 每个子任务需要4G内存

ulimit -s unlimited
echo "Running on $(hostname)"
echo "Starting at $(date)"
echo "SLURM_JOB_ID: $SLURM_JOB_ID, SLURM_ARRAY_JOB_ID: $SLURM_ARRAY_JOB_ID, SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"

# 输入和输出目录定义
input="./translated_proteins"
output="./aligned_translated_proteins"
CORES_PER_MAFFT_JOB=4 # 与 --cpus-per-task 保持一致

# 确保输出目录存在
mkdir -p $output

# 获取所有输入文件的列表
files=($input/*.fa)
total_files=${#files[@]}

if [ "$total_files" -eq 0 ]; then
    echo "No .fa files found in $input. Exiting."
    exit 0
fi

# 从数组中选择当前任务要处理的文件
# SLURM_ARRAY_TASK_ID 从 1 开始，而 bash 数组索引从 0 开始
task_id=${SLURM_ARRAY_TASK_ID}
file_index=$((task_id - 1))

# 检查任务ID是否在文件列表范围内
if [ "$file_index" -ge "$total_files" ]; then
    echo "Task ID $task_id is out of bounds. Total files: $total_files. Exiting."
    exit 0
fi

file_to_process=${files[$file_index]}
base_name=$(basename "${file_to_process%.fa}")

# --- 日志重定向 ---
# 将此任务的 stdout 和 stderr 重定向到与输入文件同名的日志文件中
# 注意：这之后的 echo 和其他命令的输出都会进入新文件
exec > "logs/${base_name}.out" 2> "logs/${base_name}.err"

echo "SLURM Job ID: ${SLURM_JOB_ID}"
echo "SLURM Array Job ID: ${SLURM_ARRAY_JOB_ID}"
echo "SLURM Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "----------------------------------------------------"
echo "Running on $(hostname)"
echo "Starting at $(date)"

echo "Processing file $task_id of $total_files: $file_to_process"
outfile="$output/$(basename "${file_to_process%.fa}")_aligned.fa"
echo "Output will be saved to $outfile"

# 运行MAFFT
singularity exec /usr/local/biotools/m/mafft:7.525--h031d066_0 mafft --auto --thread $CORES_PER_MAFFT_JOB "$file_to_process" > "$outfile"

echo "MAFFT job for $file_to_process completed."
echo "Ending at $(date)"