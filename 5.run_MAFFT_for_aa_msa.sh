#!/bin/bash
#SBATCH --job-name=mafft
#SBATCH --output=mafft_%j.out
#SBATCH --error=mafft_%j.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G

ulimit -s unlimited
echo "Running on $(hostname)"
echo "Starting at $(date)"

# 输入和输出目录定义
input="./translated_proteins_deg"
output="./aligned_translated_proteins_deg"

# 并行处理参数
# 从SBATCH指令中获取总CPU核心数 (或者直接使用SLURM_CPUS_PER_TASK)
# 这里我们根据脚本中的 #SBATCH --cpus-per-task=16 来设定
TOTAL_CPUS_ALLOCATED=32
CORES_PER_MAFFT_JOB=4
MAX_CONCURRENT_JOBS=$((TOTAL_CPUS_ALLOCATED / CORES_PER_MAFFT_JOB))

echo "Total CPUs allocated: $TOTAL_CPUS_ALLOCATED"
echo "Cores per MAFFT job: $CORES_PER_MAFFT_JOB"
echo "Maximum concurrent MAFFT jobs: $MAX_CONCURRENT_JOBS"

# 确保输出目录存在
mkdir -p $output

# 检查输入文件
echo "Checking input files in $input:"
ls $input

# 为目录中的每个文件运行MAFFT
active_jobs=0
file_count=0
total_files=$(ls $input/*.fa 2>/dev/null | wc -l)

if [ "$total_files" -eq 0 ]; then
    echo "No .fa files found in $input. Exiting."
    exit 0
fi

echo "Found $total_files files to process."

for file in $input/*.fa; do
    ((file_count++))
    echo "Processing file $file_count of $total_files: $file..."
    outfile="$output/$(basename "${file%.fa}")_aligned.fa"
    echo "Output will be saved to $outfile"
    
    # 以后台方式运行MAFFT
    singularity exec /usr/local/biotools/m/mafft:7.525--h031d066_0 mafft --auto --thread $CORES_PER_MAFFT_JOB "$file" > "$outfile" &
    
    ((active_jobs++))
    echo "Launched job for $file. Active jobs: $active_jobs."

    # 如果活动的作业数达到最大并发数，则等待一个作业完成
    if [[ "$active_jobs" -ge "$MAX_CONCURRENT_JOBS" ]]; then
        echo "Reached maximum concurrent jobs ($MAX_CONCURRENT_JOBS). Waiting for a job to finish..."
        wait -n
        ((active_jobs--))
        echo "A job finished. Active jobs: $active_jobs."
    fi
done

# 等待所有剩余的后台作业完成
echo "All MAFFT jobs launched. Waiting for remaining jobs to complete..."
wait
echo "All MAFFT jobs have completed."

echo "Ending at $(date)"