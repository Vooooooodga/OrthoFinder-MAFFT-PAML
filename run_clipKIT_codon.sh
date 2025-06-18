#!/bin/bash
#SBATCH --job-name=clipkit-array
#SBATCH --array=1-8378  # 重要提示：请将此处的数字修改为您.aln文件的确切总数
#SBATCH --output=logs/clipkit_slurm-%A_%a.out  # SLURM自身的基础日志
#SBATCH --error=logs/clipkit_slurm-%A_%a.err   # SLURM自身的基础错误日志
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=1  # ClipKIT通常是单线程的，1个核心足够
#SBATCH --mem=2G

# --- 使用说明 ---
#
# 此脚本通过SLURM作业数组高效运行成千上万个ClipKIT任务。
#
# 操作仅需两步:
#
# 1.【提交前】创建任务清单：
#   在终端运行以下命令，生成一个包含所有输入文件路径的清单。
#   (请在您存放 SCOGs_msa_codon 目录的位置运行此命令)
#
#   find "$(pwd)/SCOGs_msa_codon" -type f -name "*.aln" > clipkit_input_files.list
#
# 2.【提交任务】
#   sbatch run_clipKIT_codon.sh
#
# --- 脚本正文 ---

# 任何命令失败时立即退出，增加脚本的健壮性
set -e

# --- 配置 ---
# 定义输出修剪后文件的目录
OUTPUT_DIR="./aligned_codon_clipkit_70_coverage"
# 定义存放每个任务详细日志的目录
LOG_DIR="./logs"
# 定义任务清单文件的位置。$SLURM_SUBMIT_DIR是SLURM提供的变量，指向您运行sbatch的目录。
FILE_LIST="${SLURM_SUBMIT_DIR}/clipkit_input_files.list"
# 定义Singularity命令
CLIPKIT_CMD="singularity exec /usr/local/biotools/c/clipkit:2.3.0--pyhdfd78af_0 clipkit"


# --- 环境设置 ---
# 确保输出和日志目录存在
mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

# --- 核心执行逻辑 ---
# 1. 从SLURM获取当前子任务的唯一编号
TASK_ID=${SLURM_ARRAY_TASK_ID}

# 2. 根据任务编号，从任务清单中精确取出当前任务需要处理的文件路径
INPUT_FILE=$(sed -n "${TASK_ID}p" "$FILE_LIST")

# 3. 检查是否成功获取了文件路径
if [ -z "$INPUT_FILE" ]; then
    echo "错误：未能从任务清单 '$FILE_LIST' 的第 $TASK_ID 行获取到文件名。" >&2
    echo "请检查您的 --array 范围是否正确，以及清单文件是否存在。" >&2
    exit 1
fi

# 4. 根据输入文件名，生成一个干净的文件前缀 (例如：OG0006709_codon)
BASENAME=$(basename "$INPUT_FILE" .aln)

# 5. 定义此任务专属的日志文件和输出文件
LOG_FILE="${LOG_DIR}/${BASENAME}.clipkit.log"
OUTPUT_FILE="${OUTPUT_DIR}/${BASENAME}.clipkit.fasta"

# 6. 将此任务的所有后续输出（标准输出和错误输出）重定向到它自己的日志文件中
exec > "$LOG_FILE" 2>&1

# --- 记录任务信息 ---
echo "--- ClipKIT 子任务启动 ---"
echo "主机名: $(hostname)"
echo "任务编号 (Task ID): ${TASK_ID}"
echo "处理文件: ${INPUT_FILE}"
echo "结果输出至: ${OUTPUT_FILE}"
echo "日志记录至: ${LOG_FILE}"
echo "--------------------------"

# 7. 执行ClipKIT修剪
#    -m kpic-smart-gap: 保留简约性信息位点和恒定位点，并使用智能空位算法修剪
#    -co: 密码子模式，确保以完整密码子为单位修剪
#    -l: 生成日志文件 (与输出文件同名，但后缀为 .log)
#    -of fasta: 指定输出文件格式为FASTA
${CLIPKIT_CMD} "$INPUT_FILE" -o "$OUTPUT_FILE" -m gappy -g 0.3 -co -l -of fasta

echo "--- ClipKIT 子任务完成 ---"
echo "结束时间: $(date)"