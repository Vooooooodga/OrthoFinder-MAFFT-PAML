#!/bin/bash

# -----------------------------------------------------------------------------
# 脚本名称: extract_and_remark_m0_trees.sh
# 功能描述:
#   1. 从 PAML M0 模型输出的 .txt 文件中提取物种树（带支长和物种名）。
#   2. 使用 mark_foreground.py 脚本根据社会性数据重新标记提取的树。
#
# 如何使用:
#   1. 修改下面的 "--- 配置 ---" 部分，设置正确的路径和参数。
#   2. 给脚本执行权限: chmod +x extract_and_remark_m0_trees.sh
#   3. 运行脚本: ./extract_and_remark_m0_trees.sh
# -----------------------------------------------------------------------------

# --- 配置 ---
# !!! 用户：请仔细检查并根据您的实际情况修改以下路径和参数 !!!

# 项目基础路径 (例如: /home/user/my_project)
# 这应该是您运行 run_codeml_branch_models.sh 的项目根目录或类似路径。
PROJECT_BASE_PATH="/home/yuhangjia/data/AlternativeSplicing/9_overlap_DAS_evo_rate_analysis" # <--- 修改这里

# 包含 PAML M0 模型输出文件的目录 (通常是 branch_model 子目录)
# 文件名应为 *_M0_paml_results.txt
M0_RESULTS_DIR="${PROJECT_BASE_PATH}/branch_model_results" # <--- 修改这里 (假设M0结果在此)

# mark_foreground.py 脚本的完整路径
MARK_PY_SCRIPT="${PROJECT_BASE_PATH}/mark_foreground.py" # <--- 修改这里

# 物种社会性信息文件 (CSV 或 TSV) 的完整路径
SOCIALITY_FILE="${PROJECT_BASE_PATH}/species_sociality_data.csv" # <--- 修改这里

# 最终保存标记后树文件的基础目录
# 每个基因的标记树将保存在此目录下的一个同名子目录中
FINAL_MARKED_TREES_OUTPUT_DIR="${PROJECT_BASE_PATH}/gene_trees_from_M0_remarked"

# 用于处理的临时目录的基础路径 (脚本会自动创建和清理)
TEMP_PROCESSING_DIR_BASE="${PROJECT_BASE_PATH}/temp_m0_tree_processing_intermediate"

# mark_foreground.py 脚本所需的参数
TARGET_SOCIALITY="Advanced Eusocial" # <--- 根据需要修改 (例如: "Solitary", "Primitive Eusocial")
PAML_MARKER="#1"                     # <--- 根据需要修改 (PAML 前景标记)

# --- 配置结束 ---

echo "脚本开始：提取 M0 树并重新标记。"
echo "  项目基础路径: $PROJECT_BASE_PATH"
echo "  M0 结果目录: $M0_RESULTS_DIR"
echo "  mark_foreground.py 脚本: $MARK_PY_SCRIPT"
echo "  社会性文件: $SOCIALITY_FILE"
echo "  最终输出目录: $FINAL_MARKED_TREES_OUTPUT_DIR"
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
mkdir -p "$FINAL_MARKED_TREES_OUTPUT_DIR"
mkdir -p "$TEMP_PROCESSING_DIR_BASE"

# Awk 脚本，用于从 PAML M0 输出中提取最后一个带非数字（即名称）叶标签的 Newick 树
# 它处理可能跨越多行的树，并旨在捕获 PAML 输出中带有物种/基因名称的树
read -r -d '' AWK_SCRIPT_EXTRACT_TREE << 'EOFawk'
BEGIN {
    tree_buffer = "";
    last_valid_tree = "";
    in_tree = 0;
    # FS="\n"; # Process line by line
}
{
    current_line_content = $0;

    # 如果不在树内，并且当前行以可选空格后跟 '(' 开头，则开始一个新树
    if (!in_tree && current_line_content ~ /^\s*\(/) {
        in_tree = 1;
        # 清理前导空格
        sub(/^[ \t]+/, "", current_line_content);
        tree_buffer = current_line_content;
    } else if (in_tree) {
        # 如果在树内，则将当前行（去除首尾空格）追加到缓冲区
        # PAML输出的Newick树，如果跨行，后续行可能也有缩进。
        # 我们直接连接，因为Newick格式本身对空白不敏感（除了名称内）
        # 为了安全，先去除行首尾的空白
        gsub(/^[ \t]+|[ \t]+$/, "", current_line_content);
        tree_buffer = tree_buffer current_line_content;
    }

    # 如果在树内，并且缓冲区以可选空格后跟 ');' 结尾，则树结束
    if (in_tree && tree_buffer ~ /\);\s*$/) {
        clean_tree_buffer = tree_buffer;
        # 进一步清理整个缓冲区的尾部空白和换行符，以进行准确匹配
        gsub(/[ \t\r\n]+$/, "", clean_tree_buffer);
        gsub(/^[ \t\r\n]+/, "", clean_tree_buffer); # 也清理头部

        # 检查这是否是一个有效的带名称的树
        # 树应以 '(' 开始，后跟一个字母或下划线（表示名称的一部分）
        if (clean_tree_buffer ~ /^\([A-Za-z_]/ ) {
             # 进一步检查：确保第一个标签部分（第一个冒号之前）包含字母/下划线
             # 这有助于与 (1:0.1, 2:0.2); 这样的数字标签树区分开
             temp_label_check = clean_tree_buffer;
             sub(/^\([ \t]*/, "", temp_label_check); # 移除前导 '(' 和任何空格
             sub(/:.*/, "", temp_label_check);      # 移除从第一个 ':' 开始的所有内容

             if (temp_label_check ~ /[A-Za-z_]/) { # 如果提取的标签部分包含字母/下划线
                last_valid_tree = clean_tree_buffer;
             }
        }
        in_tree = 0;        # 重置状态
        tree_buffer = "";   # 重置缓冲区
    }
}
END {
    if (last_valid_tree != "") {
        print last_valid_tree;
    }
}
EOFawk

# 查找并处理每个 PAML M0 结果文件
# 使用 -print0 和 read -d $'\0' 来安全处理包含空格或特殊字符的文件名
find "$M0_RESULTS_DIR" -type f -name "*_M0_paml_results.txt" -print0 | while IFS= read -r -d $'\0' m0_result_file; do
    # 从 M0 结果文件名中提取基因/直系同源群的基础名称
    # 例如: 从 OG0001155_codon.clipkit_M0_paml_results.txt -> OG0001155_codon.clipkit
    gene_base_name=$(basename "$m0_result_file" "_M0_paml_results.txt")

    echo "正在处理 M0 结果文件: $m0_result_file (基因基础名: $gene_base_name)"

    # 使用 awk 脚本提取树字符串
    extracted_tree_string=$(awk "$AWK_SCRIPT_EXTRACT_TREE" "$m0_result_file")

    if [ -z "$extracted_tree_string" ]; then
        echo "  警告: 未能从 '$m0_result_file' 中提取到带名称的 Newick 树。跳过此文件。"
        continue
    fi
    echo "  成功提取树字符串。"
    # echo "  调试信息: 提取的树: $extracted_tree_string" # 取消注释以进行调试

    # 为当前基因的提取树文件准备一个临时目录
    # 此目录将包含一个 .treefile 文件，供 mark_foreground.py 处理
    current_gene_temp_tree_input_dir="${TEMP_PROCESSING_DIR_BASE}/${gene_base_name}_temp_tree_for_marking"
    mkdir -p "$current_gene_temp_tree_input_dir"

    # 定义提取的树文件的路径
    extracted_tree_filepath="${current_gene_temp_tree_input_dir}/${gene_base_name}_from_M0.treefile"
    echo "$extracted_tree_string" > "$extracted_tree_filepath"
    echo "  提取的树已保存至: $extracted_tree_filepath"

    # 定义当前基因标记后树的最终输出目录
    # mark_foreground.py 会将其输出 (例如: ${gene_base_name}_from_M0_marked.treefile) 放入此目录
    current_gene_final_output_dir="${FINAL_MARKED_TREES_OUTPUT_DIR}/${gene_base_name}_marked_M0_tree"
    mkdir -p "$current_gene_final_output_dir"

    echo "  为 $gene_base_name 运行 mark_foreground.py..."
    python3 "$MARK_PY_SCRIPT" \
        --tree_folder "$current_gene_temp_tree_input_dir" \
        --sociality_file "$SOCIALITY_FILE" \
        --output_folder "$current_gene_final_output_dir" \
        --target_sociality "$TARGET_SOCIALITY" \
        --paml_marker "$PAML_MARKER"

    if [ $? -eq 0 ]; then
        echo "  成功为 $gene_base_name 生成标记树。输出位于: $current_gene_final_output_dir"
    else
        echo "  错误: mark_foreground.py 未能成功处理 $gene_base_name。请检查相关日志（如果有）。"
        # 在此可以考虑添加更复杂的错误处理或日志记录机制
    fi

    # 清理当前基因的临时目录
    rm -rf "$current_gene_temp_tree_input_dir"
    echo "  已清理临时目录: $current_gene_temp_tree_input_dir"
    echo "-----------------------------------------------------"
done

# 清理基础临时处理目录 (此时它应该是空的)
# 2>/dev/null 用于抑制 rmdir 在目录非空时的错误消息
rmdir "$TEMP_PROCESSING_DIR_BASE" 2>/dev/null || echo "  提示: 基础临时目录 '$TEMP_PROCESSING_DIR_BASE' 非空或无法移除 (通常情况下这不影响结果)。"

echo "所有处理已完成。"
echo "最终的标记树 (源自M0输出) 保存在位于以下路径的相应子目录中: $FINAL_MARKED_TREES_OUTPUT_DIR"