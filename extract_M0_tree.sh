#!/bin/bash

# -----------------------------------------------------------------------------
# 脚本名称: extract_and_remark_m0_trees_optimized.sh
# 功能描述:
#   1. 从 PAML M0 模型输出的 .txt 文件中提取物种树（带支长和物种名）。
#   2. 将所有提取的树收集到一个临时目录。
#   3. 调用 mark_foreground.py 脚本一次，对临时目录中所有的树文件
#      根据社会性数据进行前景标记。
#
# 如何使用:
#   1. 修改下面的 "--- 配置 ---" 部分，设置正确的路径和参数。
#   2. 给脚本执行权限: chmod +x extract_and_remark_m0_trees_optimized.sh
#   3. 运行脚本: ./extract_and_remark_m0_trees_optimized.sh
# -----------------------------------------------------------------------------

# --- 配置 ---
# !!! 用户：请仔细检查并根据您的实际情况修改以下路径和参数 !!!

# 项目基础路径 (例如: /home/user/my_project)
PROJECT_BASE_PATH="/Users/jiayuhang/Desktop/doc/OrthoFinder-MAFFT-PAML" # <--- 修改这里

# 包含 PAML M0 模型输出文件的目录 (通常是 branch_model_results 子目录)
# 文件名应为 *_M0_paml_results.txt
M0_RESULTS_DIR="${PROJECT_BASE_PATH}/branch_model_results" # <--- 修改这里 (假设M0结果在此)

# mark_foreground.py 脚本的完整路径
MARK_PY_SCRIPT="${PROJECT_BASE_PATH}/mark_foreground.py" # <--- 修改这里

# 物种社会性信息文件 (CSV 或 TSV) 的完整路径
SOCIALITY_FILE="${PROJECT_BASE_PATH}/species_sociality_data.csv" # <--- 修改这里

# 最终保存标记后树文件的目录
# 所有标记后的树文件将直接保存在此目录中
FINAL_MARKED_TREES_OUTPUT_DIR="${PROJECT_BASE_PATH}/gene_trees_from_M0_remarked_all"

# 用于存放所有提取出来的原始树文件的临时目录 (脚本会自动创建和清理)
TEMP_EXTRACTED_TREES_DIR="${PROJECT_BASE_PATH}/temp_all_extracted_m0_trees"

# mark_foreground.py 脚本所需的参数
TARGET_SOCIALITY="Advanced Eusocial" # <--- 根据需要修改 (例如: "Solitary", "Primitive Eusocial")
PAML_MARKER="#1"                     # <--- 根据需要修改 (PAML 前景标记)

# --- 配置结束 ---

echo "脚本开始：提取 M0 树并重新标记 (优化版)。"
echo "  项目基础路径: $PROJECT_BASE_PATH"
echo "  M0 结果目录: $M0_RESULTS_DIR"
echo "  mark_foreground.py 脚本: $MARK_PY_SCRIPT"
echo "  社会性文件: $SOCIALITY_FILE"
echo "  最终输出目录: $FINAL_MARKED_TREES_OUTPUT_DIR"
echo "  临时提取树目录: $TEMP_EXTRACTED_TREES_DIR"
echo "  目标社会性: $TARGET_SOCIALITY"
echo "  PAML 标记: $PAML_MARKER"
echo "-----------------------------------------------------"

# --- 预检查 ---
if [ ! -d "$PROJECT_BASE_PATH" ]; then
  echo "错误: 项目基础路径 '$PROJECT_BASE_PATH' 未找到。"
  exit 1
fi
if [ ! -d "$M0_RESULTS_DIR" ]; then
  echo "错误: M0 结果目录 '$M0_RESULTS_DIR' 未找到。"
  exit 1
fi
if [ ! -f "$MARK_PY_SCRIPT" ]; then
  echo "错误: mark_foreground.py 脚本 '$MARK_PY_SCRIPT' 未找到。"
  exit 1
fi
if [ ! -f "$SOCIALITY_FILE" ]; then
  echo "错误: 社会性文件 '$SOCIALITY_FILE' 未找到。"
  exit 1
fi
if ! command -v python3 &> /dev/null; then
    echo "错误: python3 命令未找到。请确保 Python 3 已安装并位于 PATH中。"
    exit 1
fi
if ! command -v awk &> /dev/null; then
    echo "错误: awk 命令未找到。请确保 awk 已安装并位于 PATH中。"
    exit 1
fi

# --- 主要逻辑 ---
# 创建或清空临时提取树目录和最终输出目录
mkdir -p "$FINAL_MARKED_TREES_OUTPUT_DIR"
rm -rf "$TEMP_EXTRACTED_TREES_DIR" # 清理旧的临时文件（如果存在）
mkdir -p "$TEMP_EXTRACTED_TREES_DIR"

# Awk 脚本，用于从 PAML M0 输出中提取最后一个带非数字（即名称）叶标签的 Newick 树
read -r -d '' AWK_SCRIPT_EXTRACT_TREE << 'EOFawk'
BEGIN {
    tree_buffer = "";
    last_valid_tree = "";
    in_tree = 0;
}
{
    current_line_content = $0;
    if (!in_tree && current_line_content ~ /^\s*\(/) {
        in_tree = 1;
        sub(/^[ \t]+/, "", current_line_content);
        tree_buffer = current_line_content;
    } else if (in_tree) {
        gsub(/^[ \t]+|[ \t]+$/, "", current_line_content);
        tree_buffer = tree_buffer current_line_content;
    }
    if (in_tree && tree_buffer ~ /\);\s*$/) {
        clean_tree_buffer = tree_buffer;
        gsub(/[ \t\r\n]+$/, "", clean_tree_buffer);
        gsub(/^[ \t\r\n]+/, "", clean_tree_buffer);
        if (clean_tree_buffer ~ /^\([A-Za-z_]/ ) {
             temp_label_check = clean_tree_buffer;
             sub(/^\([ \t]*/, "", temp_label_check);
             sub(/:.*/, "", temp_label_check);
             if (temp_label_check ~ /[A-Za-z_]/) {
                last_valid_tree = clean_tree_buffer;
             }
        }
        in_tree = 0;
        tree_buffer = "";
    }
}
END {
    if (last_valid_tree != "") {
        print last_valid_tree;
    }
}
EOFawk

echo "步骤 1: 提取所有 M0 树到临时目录..."
extracted_tree_count=0
# 查找并处理每个 PAML M0 结果文件
find "$M0_RESULTS_DIR" -type f -name "*_M0_paml_results.txt" -print0 | while IFS= read -r -d $'\0' m0_result_file; do
    gene_base_name=$(basename "$m0_result_file" "_M0_paml_results.txt")
    echo "  正在从 $m0_result_file 提取树 (基因: $gene_base_name)..."

    extracted_tree_string=$(awk "$AWK_SCRIPT_EXTRACT_TREE" "$m0_result_file")

    if [ -z "$extracted_tree_string" ]; then
        echo "    警告: 未能从 '$m0_result_file' 中提取到带名称的 Newick 树。跳过。"
        continue
    fi

    # 定义提取的树文件在临时目录中的路径
    # 文件名将是类似 OG0001155_codon.clipkit_from_M0.treefile
    extracted_tree_filepath="${TEMP_EXTRACTED_TREES_DIR}/${gene_base_name}_from_M0.treefile"
    echo "$extracted_tree_string" > "$extracted_tree_filepath"
    echo "    提取的树已保存至: $extracted_tree_filepath"
    extracted_tree_count=$((extracted_tree_count + 1))
done

echo "步骤 1 完成: 共提取了 $extracted_tree_count 个树文件到 $TEMP_EXTRACTED_TREES_DIR"
echo "-----------------------------------------------------"

if [ "$extracted_tree_count" -eq 0 ]; then
    echo "错误: 未能从 $M0_RESULTS_DIR 中的任何文件提取树。无法继续进行标记。"
    rm -rf "$TEMP_EXTRACTED_TREES_DIR" # 清理空的临时目录
    exit 1
fi

echo "步骤 2: 运行 mark_foreground.py 对所有提取的树进行标记..."
python3 "$MARK_PY_SCRIPT" \
    --tree_folder "$TEMP_EXTRACTED_TREES_DIR" \
    --sociality_file "$SOCIALITY_FILE" \
    --output_folder "$FINAL_MARKED_TREES_OUTPUT_DIR" \
    --target_sociality "$TARGET_SOCIALITY" \
    --paml_marker "$PAML_MARKER"

if [ $? -eq 0 ]; then
    echo "步骤 2 完成: mark_foreground.py 成功执行。"
    echo "标记后的树文件已保存到: $FINAL_MARKED_TREES_OUTPUT_DIR"
else
    echo "错误: mark_foreground.py 执行失败。请检查相关日志（如果有）。"
    echo "最终输出可能不完整或不存在于: $FINAL_MARKED_TREES_OUTPUT_DIR"
    # 保留临时目录以便调试
    echo "提示: 提取的原始树文件保留在 $TEMP_EXTRACTED_TREES_DIR 以便进行调试。"
    exit 1 # 发生错误时退出
fi
echo "-----------------------------------------------------"

echo "步骤 3: 清理临时文件..."
rm -rf "$TEMP_EXTRACTED_TREES_DIR"
echo "  已清理临时目录: $TEMP_EXTRACTED_TREES_DIR"
echo "-----------------------------------------------------"

echo "所有处理已完成。"
echo "最终的标记树 (源自M0输出) 保存在: $FINAL_MARKED_TREES_OUTPUT_DIR"