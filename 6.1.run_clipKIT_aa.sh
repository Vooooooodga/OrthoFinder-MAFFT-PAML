#!/bin/bash
#SBATCH --job-name=clipkit_trim
#SBATCH --output=clipkit_trim_aa_%j.out
#SBATCH --error=clipkit_trim_aa_%j.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G

# 定义每个ClipKIT任务使用的CPU核心数
CPUS_PER_TASK=4

ulimit -s unlimited
echo "Running on $(hostname)"
echo "Starting at $(date)"

# 定义Singularity命令前缀
CLIPKIT_CMD="singularity exec /usr/local/biotools/c/clipkit:2.3.0--pyhdfd78af_0 clipkit"

# 输入和输出目录定义
input_dir="./SCOGs_msa_aa"
output_dir="./SCOGs_msa_aa_clipkit" # 修剪后的FASTA文件将存放在新目录中

# 确保输出目录存在
mkdir -p "$output_dir"

# 检查输入目录是否存在且包含文件
if [ ! -d "$input_dir" ]; then
    echo "Error: Input directory $input_dir does not exist."
    exit 1
fi

# 使用 shopt -s nullglob 使得在没有匹配文件时，通配符扩展为空，而不是保留通配符本身
shopt -s nullglob
files_to_process=("$input_dir"/*.fa)
shopt -u nullglob # (可选) 如果后续脚本部分不需要此行为，可以取消设置

# 检查是否找到了 .fa 文件
if [ ${#files_to_process[@]} -eq 0 ]; then
    echo "Error: No .fa files found in $input_dir"
    exit 1
fi

echo "Found ${#files_to_process[@]} .fa files to process in $input_dir:"

# 获取总可用CPU核心数
TOTAL_CPUS=$(nproc)
# 计算可以同时运行的最大作业数
MAX_PARALLEL_JOBS=$((TOTAL_CPUS / CPUS_PER_TASK))
echo "Total available CPUs: $TOTAL_CPUS"
echo "CPUs per task: $CPUS_PER_TASK"
echo "Maximum parallel jobs: $MAX_PARALLEL_JOBS"

# 为目录中的每个 .aln 文件运行 ClipKIT
job_count=0
for file_path in "${files_to_process[@]}"; do
    # 再次确认是文件 (通常nullglob后不需要，但作为安全检查)
    if [ -f "$file_path" ]; then
        filename=$(basename "$file_path")
        # 构建输出文件名，例如：OG0006709_codon.clipkit.fasta
        base_name_no_ext=$(basename "${filename%.fa}") # 去除 .aln 后缀
        outfile="$output_dir/${base_name_no_ext}.clipkit.fasta"

        echo "Processing $filename ..."
        echo "Input: $file_path"
        echo "Output will be saved to $outfile (FASTA format)"

        # 构建CPU亲和性掩码，这里我们简单地分配连续的CPU核心
        # 注意：这种分配方式在高负载或复杂NUMA架构下可能不是最优的
        # 但对于多数情况是一个合理的起点。
        # 我们需要确保分配的CPU核心不重复。
        # 以下的 taskset 逻辑比较复杂，暂时先不为每个任务精确分配CPU核心，
        # 而是依赖于系统的调度。如果需要更精细的控制，可以后续添加。
        # SLURM的--cpus-per-task通常会处理好CPU分配。
        # 如果不通过SLURM直接运行，并且ClipKIT本身不支持多线程，
        # 那么并行运行多个单线程的ClipKIT实例是合理的。
        # 此处我们假设ClipKIT是单线程的，并行启动多个实例。

        # 运行 ClipKIT (后台运行)
        # -m kpic-smart-gap: 保留简约性信息位点和恒定位点，并使用智能空位算法修剪
        # -co: 密码子模式，确保以完整密码子为单位修剪
        # -l: 生成日志文件 (与输出文件同名，但后缀为 .log)
        # -of fasta: 指定输出文件格式为FASTA
        (
          # 如果需要严格限制每个进程使用的CPU，可以使用taskset
          # 这里我们先不使用taskset，因为SLURM的--cpus-per-task应该已经处理了资源分配
          # 如果直接在没有SLURM的机器上运行，并且希望限制每个clipkit进程的CPU，可以取消下面的注释
          # cpu_start=$(( (job_count % MAX_PARALLEL_JOBS) * CPUS_PER_TASK ))
          # cpu_end=$(( cpu_start + CPUS_PER_TASK - 1 ))
          # taskset -c $cpu_start-$cpu_end ${CLIPKIT_CMD} "$file_path" -o "$outfile" -m kpic-smart-gap -co -l -of fasta
          ${CLIPKIT_CMD} "$file_path" -o "$outfile" -m kpic-smart-gap -l -of fasta
          echo "Finished processing $filename. Log file: ${outfile}.log"
          echo "-----------------------------------------------------"
        ) & # 将命令放入后台执行

        # 控制并行作业的数量
        job_count=$((job_count + 1))
        if (( job_count % MAX_PARALLEL_JOBS == 0 )); then
            echo "Waiting for a batch of $MAX_PARALLEL_JOBS jobs to complete..."
            wait # 等待所有后台作业完成
            echo "Batch completed."
        fi
    fi
done

# 等待剩余的后台作业完成
echo "Waiting for any remaining jobs to complete..."
wait
echo "All ClipKIT trimming jobs finished."
echo "Ending at $(date)"