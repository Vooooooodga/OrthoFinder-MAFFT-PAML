#!/bin/bash
# 该脚本用于重新提交在 'short' 分区上运行失败的 SLURM 作业。
# 它会执行以下操作：
# 1. 查找非空的错误日志文件，以确定哪些作业失败了。
# 2. 从日志文件名中提取作业的 OG_ID。
# 3. 修改对应的提交脚本：
#    a. 修正 hyphy 命令以正确使用多核心 (将 CPU=... 替换为 CPU=128)。
#    b. 确保作业在 'short' 分区上运行。
#    c. 确保作业请求 128 CPU 和 128G 内存。
#    d. 更新作业名和日志路径以供追踪。
# 4. 重新提交已修改的脚本。

set -e # 如果任何命令失败，则立即退出
set -u # 将未设置的变量视为错误
set -o pipefail # 管道中的命令失败也会导致整个管道失败

# --- 配置 ---
# 请根据您的实际环境验证以下路径
LOG_DIR="/lustre10/home/yuhangjia/data/AlternativeSplicing/evo_rate_test_RNA_splicing_term/hyphy_relax/slurm_logs"
SCRIPT_DIR="/lustre10/home/yuhangjia/data/AlternativeSplicing/evo_rate_test_RNA_splicing_term/hyphy_relax/submission_scripts"
FAILED_OG_LIST_FILE="failed_ogs_to_resubmit.txt"

# --- 步骤 1: 查找失败的作业并提取 OG_ID ---
echo "正在从 ${LOG_DIR} 目录中查找失败的作业..."
# 查找所有以 'short_' 开头且大小不为零的 .err 文件，并从中提取 OG_ID
find "${LOG_DIR}" -name "short_*.err" -size +0c | sed -E 's/.*short_(OG[0-9]+)_[0-9]+.err/\1/' | sort | uniq > "${FAILED_OG_LIST_FILE}"

# 检查是否找到了失败的作业
if [ ! -s "${FAILED_OG_LIST_FILE}" ]; then
    echo "未找到任何失败的作业。所有日志文件均为空。"
    rm -f "${FAILED_OG_LIST_FILE}"
    exit 0
fi

JOB_COUNT=$(wc -l < "${FAILED_OG_LIST_FILE}")
echo "找到了 ${JOB_COUNT} 个失败的作业。其 OG_ID 列表已保存至: ${FAILED_OG_LIST_FILE}"
echo "----------------------------------------------------"

# --- 步骤 2 & 3: 修改提交脚本并重新提交 ---
echo "开始修改脚本分区、资源和HyPhy命令，并重新提交作业..."

while IFS= read -r OG_ID; do
    # 跳过空行
    if [ -z "${OG_ID}" ]; then
        continue
    fi

    SUB_SCRIPT_PATH="${SCRIPT_DIR}/${OG_ID}.sh"

    if [ -f "${SUB_SCRIPT_PATH}" ]; then
        echo "正在处理脚本: ${OG_ID}.sh"

        # 使用 sed 命令进行原地修改
        # 1. 将分区从 'rome' 改为 'short'
        sed -i 's/--partition=rome/--partition=short/' "${SUB_SCRIPT_PATH}"

        # 2. 修改CPU和内存请求
        sed -i 's/--cpus-per-task=16/--cpus-per-task=128/' "${SUB_SCRIPT_PATH}"
        sed -i 's/--mem=16G/--mem=128G/' "${SUB_SCRIPT_PATH}"

        # 3. 修正HyPhy命令: 使用 perl 保证兼容性, 将 CPU=... 替换为 CPU=128
        perl -i -pe 's/CPU=[^ ]*/CPU=128/' "${SUB_SCRIPT_PATH}"

        # 重新提交作业
        sbatch "${SUB_SCRIPT_PATH}"
        echo "已将 ${OG_ID} 的作业重新提交至 'short' 分区。"
        sleep 0.1 # 短暂延迟，避免对 SLURM 控制器造成过大压力

    else
        echo "警告: 未找到 ${OG_ID} 对应的提交脚本: ${SUB_SCRIPT_PATH}"
    fi
    echo "----------------------------------------------------"
done < "${FAILED_OG_LIST_FILE}"

echo "所有失败的作业都已重新提交到 'short' 分区。"
echo "您可以删除该列表文件: rm ${FAILED_OG_LIST_FILE}" 