#!/bin/bash
#SBATCH --job-name=mafft-array
#SBATCH --array=1-8378  # 重要提示：请将此处的数字修改为您文件的确切总数
#SBATCH --output=logs/slurm-%A_%a.out  # SLURM自身的基础日志
#SBATCH --error=logs/slurm-%A_%a.err   # SLURM自身的基础错误日志
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G

# --- 使用说明 ---
#
# 此脚本通过SLURM作业数组高效运行成千上万个MAFFT任务。
#
# 操作仅需两步:
#
# 1.【提交前】创建任务清单：
#   在终端运行以下命令，生成一个包含所有输入文件路径的清单。
#   (请在您存放 translated_proteins 目录的位置运行此命令)
#
#   find "$(pwd)/translated_proteins" -type f -name "*.fa" > mafft_input_files.list
#
# 2.【提交任务】
#   sbatch 5.run_MAFFT_for_aa_msa.sh
#
# --- 脚本正文 ---

# 任何命令失败时立即退出，增加脚本的健壮性
set -e

# --- 配置 ---
# 定义输出比对结果的目录
OUTPUT_DIR="./aligned_translated_proteins"
# 定义存放每个任务详细日志的目录
LOG_DIR="./logs"
# 定义任务清单文件的位置。$SLURM_SUBMIT_DIR是SLURM提供的变量，指向您运行sbatch的目录。
FILE_LIST="${SLURM_SUBMIT_DIR}/mafft_input_files.list"

# --- 环境设置 ---
# 确保输出和日志目录存在
mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

# --- 核心执行逻辑 ---
# 1. 从SLURM获取当前子任务的唯一编号
TASK_ID=${SLURM_ARRAY_TASK_ID}

# 2. 根据任务编号，从任务清单中精确取出当前任务需要处理的文件路径
#    sed -n "${TASK_ID}p" 的作用是打印清单文件中的第 N 行，N就是任务编号
INPUT_FILE=$(sed -n "${TASK_ID}p" "$FILE_LIST")

# 3. 检查是否成功获取了文件路径
if [ -z "$INPUT_FILE" ]; then
    echo "错误：未能从任务清单 '$FILE_LIST' 的第 $TASK_ID 行获取到文件名。" >&2
    echo "请检查您的 --array 范围是否正确。" >&2
    exit 1
fi

# 4. 根据输入文件名，生成一个干净的文件前缀 (例如：OG0007311)
BASENAME=$(basename "$INPUT_FILE" .fa)

# 5. 定义此任务专属的日志文件和输出文件
LOG_FILE="${LOG_DIR}/${BASENAME}.log"
OUTPUT_FILE="${OUTPUT_DIR}/${BASENAME}_aligned.fa"

# 6. 将此任务的所有后续输出（标准输出和错误输出）重定向到它自己的日志文件中
exec > "$LOG_FILE" 2>&1

# --- 记录任务信息 ---
echo "--- MAFFT 子任务启动 ---"
echo "主机名: $(hostname)"
echo "任务编号 (Task ID): ${TASK_ID}"
echo "处理文件: ${INPUT_FILE}"
echo "结果输出至: ${OUTPUT_FILE}"
echo "日志记录至: ${LOG_FILE}"
echo "--------------------------"

# 7. 执行MAFFT比对
#    --thread 参数直接使用SLURM分配的CPU核心数，无需手动指定
singularity exec /usr/local/biotools/m/mafft:7.525--h031d066_0 mafft --auto --thread "$SLURM_CPUS_PER_TASK" "$INPUT_FILE" > "$OUTPUT_FILE"

echo "--- MAFFT 子任务完成 ---"
echo "结束时间: $(date)"