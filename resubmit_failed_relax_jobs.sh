#!/bin/bash
# 该脚本用于重新提交失败的 HyPhy RELAX 作业。
# 流程:
# 1. 识别失败的作业：
#    a. 遍历 /home/yuhangjia/data/AlternativeSplicing/all_OGs_for_relax_test/subtrees_with_foreground/ 中的所有 OGxxxxxx_marked.treefile 文件作为作业全集。
#    b. 检查对应的 /home/yuhangjia/data/AlternativeSplicing/all_OGs_for_relax_test/hyphy_relax_submission_20250621_140935/relax_json_results/OGxxxxxx.RELAX.json 文件是否存在且非空。
#    c. 如果 JSON 文件不存在或为空，则将该 OG_ID 标记为失败。
# 2. 为失败的作业重新生成提交脚本并准备提交：
#    a. 根据用户提供的模板生成新的 sbatch 脚本。
#    b. 将作业提交到 'short' 分区，并配置高资源（128 CPU, 128G 内存）。
#    c. sbatch 命令默认被注释，需要用户手动取消注释以执行提交。

set -e
set -u
set -o pipefail

# --- 配置 ---
# !!! 用户请务必核对以下路径是否正确 !!!
# 1. 包含前景枝标记的树文件目录
TREE_DIR="/lustre10/home/yuhangjia/data/AlternativeSplicing/all_OGs_for_relax_test/subtrees_with_foreground"
# 2. HyPhy RELAX 的 JSON 结果输出目录
JSON_RESULT_DIR="/lustre10/home/yuhangjia/data/AlternativeSplicing/all_OGs_for_relax_test/hyphy_relax_submission_20250621_140935/relax_json_results"
# 3. 多重序列比对文件目录
MSA_DIR="/lustre10/home/yuhangjia/data/AlternativeSplicing/all_OGs_for_relax_test/aligned_codon_replace_stop_codon_70_coverage"
# 4. 为重新提交而生成的新脚本的存放目录
RESUBMIT_SCRIPT_DIR="/lustre10/home/yuhangjia/data/AlternativeSplicing/all_OGs_for_relax_test/hyphy_relax_submission_20250621_140935/submission_scripts/short_sub_scripts"
# 5. SLURM 日志文件目录
LOG_DIR="/lustre10/home/yuhangjia/data/AlternativeSplicing/all_OGs_for_relax_test/hyphy_relax_submission_20250621_140935/slurm_logs"
# 6. 存储失败任务 OG_ID 的临时文件
FAILED_OG_LIST_FILE="failed_relax_ogs_to_resubmit.txt"

# --- 步骤 1: 查找所有失败或未完成的作业 ---
echo "开始扫描所有作业，以识别失败或未完成的任务..."
# 如果存在旧的失败列表，先清空
> "${FAILED_OG_LIST_FILE}"

# 遍历所有树文件，这是所有应有任务的基准
for treefile in "${TREE_DIR}"/OG*_marked.treefile; do
    # 如果目录为空或没有匹配的文件，find可能会返回路径本身
    if [ ! -f "${treefile}" ]; then
        continue
    fi

    # 从树文件名中提取 OG_ID (例如: OG0000001_marked.treefile -> OG0000001)
    OG_ID=$(basename "${treefile}" | sed -E 's/_marked\.treefile//')
    
    JSON_FILE="${JSON_RESULT_DIR}/${OG_ID}.RELAX.json"

    # 核心判断：检查对应的 JSON 结果文件是否存在且大小不为零
    if [ -s "${JSON_FILE}" ]; then
        # 文件存在且非空，视为成功，不做任何事
        :
    else
        # 文件不存在或为空，视为失败
        echo "检测到失败/未完成的作业: ${OG_ID}"
        echo "${OG_ID}" >> "${FAILED_OG_LIST_FILE}"
    fi
done

# 检查是否找到了失败的作业
if [ ! -s "${FAILED_OG_LIST_FILE}" ]; then
    echo "----------------------------------------------------"
    echo "恭喜！未找到任何失败或未完成的作业。所有任务均已成功完成。"
    rm -f "${FAILED_OG_LIST_FILE}" # 清理空文件
    exit 0
fi

JOB_COUNT=$(wc -l < "${FAILED_OG_LIST_FILE}")
echo "----------------------------------------------------"
echo "扫描完成。共找到 ${JOB_COUNT} 个失败或未完成的作业。"
echo "其 OG_ID 列表已保存至: ${FAILED_OG_LIST_FILE}"
echo "----------------------------------------------------"

# --- 步骤 2: 为失败的作业生成新的提交脚本 ---
echo "开始为所有失败的作业生成新的 SLURM 提交脚本..."

# 确保输出目录和日志目录存在
mkdir -p "${RESUBMIT_SCRIPT_DIR}"
mkdir -p "${LOG_DIR}"

while IFS= read -r OG_ID; do
    # 跳过可能的空行
    if [ -z "${OG_ID}" ]; then
        continue
    fi

    # 定义各个文件的完整路径
    SUB_SCRIPT="${RESUBMIT_SCRIPT_DIR}/${OG_ID}.sh"
    TREE_FILE="${TREE_DIR}/${OG_ID}_marked.treefile"
    MSA_FILE="${MSA_DIR}/${OG_ID}_codon.clipkit.fasta"
    JSON_OUTPUT_FILE="${JSON_RESULT_DIR}/${OG_ID}.RELAX.json"

    # 在创建脚本之前，最终检查所需文件是否都齐全
    if [ ! -f "${MSA_FILE}" ]; then
        echo "警告: 找不到 MSA 文件，跳过为 ${OG_ID} 生成脚本。路径: ${MSA_FILE}"
        continue
    fi
    if [ ! -f "${TREE_FILE}" ]; then
        echo "警告: 找不到树文件，跳过为 ${OG_ID} 生成脚本。路径: ${TREE_FILE}"
        continue
    fi

    echo "正在为 ${OG_ID} 生成脚本: ${SUB_SCRIPT}"
    
    # 使用 heredoc 创建 sbatch 脚本
    cat << SBATCH_EOF > "${SUB_SCRIPT}"
#!/bin/bash
#SBATCH --job-name=relax_resub_${OG_ID}
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=128G
#SBATCH --output=${LOG_DIR}/resub_short_${OG_ID}_%j.out
#SBATCH --error=${LOG_DIR}/resub_short_${OG_ID}_%j.err

echo "重新运行 HyPhy RELAX for ${OG_ID} on \$(hostname)"
echo "Start time: \$(date)"

singularity exec -B /lustre10:/lustre10 "/usr/local/biotools/h/hyphy:2.5.65--he91c24d_0" hyphy CPU=128 relax \\
    --alignment "${MSA_FILE}" \\
    --tree "${TREE_FILE}" \\
    --test test \\
    --output "${JSON_OUTPUT_FILE}" \\
    --code Universal < /dev/null

EXIT_CODE=\$?
echo "Job for ${OG_ID} finished with exit code \${EXIT_CODE}."
echo "End time: \$(date)"
SBATCH_EOF

    # 给予执行权限
    chmod +x "${SUB_SCRIPT}"

    # 提交作业（已按您的要求注释掉）
    echo "-> 已为 ${OG_ID} 生成脚本。如需提交, 请在本脚本中取消下面一行的注释。"
    # sbatch "${SUB_SCRIPT}"
    
    sleep 0.05 # 短暂延迟，避免文件系统操作过快

done < "${FAILED_OG_LIST_FILE}"

echo "----------------------------------------------------"
echo "所有失败作业的重新提交脚本均已在 ${RESUBMIT_SCRIPT_DIR} 目录中生成完毕。"
echo "重要提示:"
echo "1. 请务必检查脚本顶部的 MSA_DIR 和 MSA_SUFFIX 变量是否正确。"
echo "2. 请抽查几个位于 ${RESUBMIT_SCRIPT_DIR} 的脚本，确保其内容无误。"
echo "3. 确认无误后，请取消本脚本中 'sbatch \"\${SUB_SCRIPT}\"' 行的注释，然后再次运行本脚本以批量提交所有任务。"
echo "4. 任务提交后，您可以删除临时列表文件: rm ${FAILED_OG_LIST_FILE}"
echo "----------------------------------------------------" 