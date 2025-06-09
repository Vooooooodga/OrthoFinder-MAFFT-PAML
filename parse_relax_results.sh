#!/bin/bash

# --- 使用说明 ---
# 该脚本用于解析 HyPhy RELAX 分析生成的 JSON 结果文件，
# 并将关键结果 (p-value, k-value) 汇总到一个 CSV 文件中。
#
# 用法:
# bash parse_relax_results.sh <json_directory> <output_csv_file>
#
# 参数:
#   <json_directory>: 存放 RELAX 输出的 .json 文件的目录。
#                     (在 run_hyphy_relax.sh 中是 RELAX_JSON_OUTPUT_DIR)
#   <output_csv_file>:  输出的 CSV 文件路径。
#
# 依赖:
#   - jq: 一个轻量级的命令行 JSON 解析器。请确保已经安装。
#     (例如, 在 Ubuntu/Debian 上: sudo apt-get install jq)
# ---

# --- 输入验证 ---
if [ "$#" -ne 2 ]; then
    echo "错误: 需要提供两个参数。"
    echo "用法: bash $0 <json_directory> <output_csv_file>"
    exit 1
fi

JSON_DIR=$1
OUTPUT_CSV=$2

if [ ! -d "${JSON_DIR}" ]; then
    echo "错误: 输入目录 '${JSON_DIR}' 不存在。"
    exit 1
fi

# 检查 jq 是否安装
if ! command -v jq &> /dev/null; then
    echo "错误: 核心依赖 'jq' 未找到。请安装 jq 后再运行此脚本。"
    echo "例如, 在 Ubuntu/Debian 上: sudo apt-get install jq"
    echo "在 CentOS/RHEL 上: sudo yum install jq"
    echo "在 macOS (使用 Homebrew) 上: brew install jq"
    exit 1
fi

# --- 主逻辑 ---
echo "开始解析 HyPhy RELAX 结果..."
echo "输入目录: ${JSON_DIR}"
echo "输出文件: ${OUTPUT_CSV}"

# 创建或清空输出文件并写入 CSV 头部
echo "Gene,p-value,k" > "${OUTPUT_CSV}"

# 遍历目录中所有的 .RELAX.json 文件
# 使用 find 命令可以更好地处理没有匹配文件的情况
find "${JSON_DIR}" -maxdepth 1 -name "*.RELAX.json" -print0 | while IFS= read -r -d $'\0' json_file; do
    # 提取基因ID (从文件名中去除 .RELAX.json 后缀)
    base_name=$(basename "${json_file}")
    gene_id=${base_name%.RELAX.json}

    echo "正在处理: ${base_name}"

    # 使用 jq 解析 p-value 和 k 参数 (relaxation or intensification parameter)
    # -r 选项移除输出字符串的双引号
    p_value=$(jq -r '."test results"."p-value"' "${json_file}")
    k_value=$(jq -r '."test results"."relaxation or intensification parameter"' "${json_file}")

    # 检查 jq 是否成功解析
    if [ "$p_value" == "null" ] || [ "$k_value" == "null" ]; then
        echo "警告: 未能在文件 '${base_name}' 中找到 'p-value' 或 'k' 值。可能是个空文件或格式错误。跳过此文件。"
        # 可以选择写入错误日志或在主CSV中标记为错误
        echo "${gene_id},NA,NA" >> "${OUTPUT_CSV}"
    else
        # 将结果追加到 CSV 文件
        echo "${gene_id},${p_value},${k_value}" >> "${OUTPUT_CSV}"
    fi
done

echo "-----------------------------------------------------"
echo "解析完成！结果已保存到 ${OUTPUT_CSV}"

# 检查输出文件是否只包含头部 (即没有找到任何json文件)
if [ "$(wc -l < "${OUTPUT_CSV}")" -le 1 ]; then
    echo "警告: 没有在目录 '${JSON_DIR}' 中找到任何 '*.RELAX.json' 文件，或者所有文件都无法解析。"
fi 