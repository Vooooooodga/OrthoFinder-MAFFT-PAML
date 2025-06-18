#!/bin/bash
#
# 1_prepare_resampling_list.sh
#
# 功能:
# 1. 从源目录过滤掉在排除列表中的文件。
# 2. 对过滤后的文件列表进行多次重抽样。
# 3. 合并所有抽样结果，去重，生成一个最终的、唯一的文件列表。
#
# 使用方法:
# bash 1_prepare_resampling_list.sh
#

# --- 用户可配置变量 ---

# 输入包含 FASTA 文件的目录
INPUT_DIR="./aligned_codon_clipkit_70_coverage"
# 需要排除的 OG 列表文件
EXCLUSION_FILE="test_OGs.txt"
# 存放抽样列表的目录
SAMPLING_LIST_DIR="./sampling_lists_codon_resampling"
# 去重后最终待处理的文件列表
UNIQUE_FILE_LIST="unique_sampled_files_codon_resampling.txt"

# 抽样参数
NUM_SAMPLES=147
NUM_REPEATS=100

# --- 脚本主逻辑 ---
echo "=================================================="
echo "          开始准备文件列表 (准备阶段)           "
echo "=================================================="

# 检查输入目录
if [ ! -d "$INPUT_DIR" ]; then
    echo "错误: 输入目录 '$INPUT_DIR' 不存在。"
    exit 1
fi
# 检查排除文件
if [ ! -f "$EXCLUSION_FILE" ]; then
    echo "错误: 排除文件 '$EXCLUSION_FILE' 不存在。"
    exit 1
fi

# 创建输出目录
mkdir -p "$SAMPLING_LIST_DIR"

# 1. 获取所有符合条件的 fasta 文件
echo "步骤 1: 从 '$INPUT_DIR' 获取所有 .clipkit.fasta 文件..."
find "$INPUT_DIR" -type f -name "*.clipkit.fasta" > all_files.tmp

# 2. 根据 EXCLUSION_FILE 过滤文件
echo "步骤 2: 根据 '$EXCLUSION_FILE' 过滤文件..."
grep -v -f "$EXCLUSION_FILE" all_files.tmp > filtered_files.tmp
FILTERED_COUNT=$(wc -l < filtered_files.tmp)
echo "过滤后剩下 $FILTERED_COUNT 个文件可用于抽样。"

if [ "$FILTERED_COUNT" -lt "$NUM_SAMPLES" ]; then
    echo "错误：过滤后的文件数量 ($FILTERED_COUNT) 少于单次抽样所需的数量 ($NUM_SAMPLES)。"
    rm all_files.tmp filtered_files.tmp
    exit 1
fi

# 3. 执行100次重抽样
echo "步骤 3: 执行 ${NUM_REPEATS} 次重抽样, 每次抽取 ${NUM_SAMPLES} 个文件..."
for i in $(seq 1 $NUM_REPEATS); do
    shuf -n $NUM_SAMPLES filtered_files.tmp > "$SAMPLING_LIST_DIR/sample_${i}.txt"
done
echo "抽样列表已生成在 '$SAMPLING_LIST_DIR' 目录中。"

# 4. 合并抽样结果并去重
echo "步骤 4: 合并所有抽样列表并去重..."
cat "$SAMPLING_LIST_DIR"/sample_*.txt | sort | uniq > "$UNIQUE_FILE_LIST"

UNIQUE_COUNT=$(wc -l < "$UNIQUE_FILE_LIST")
echo "合并去重后，总共有 $UNIQUE_COUNT 个独立文件需要处理。"
echo "最终文件列表保存在: $UNIQUE_FILE_LIST"

# 清理临时文件
rm all_files.tmp filtered_files.tmp

echo "=================================================="
echo "准备阶段完成！"
echo "现在，您可以运行第二个脚本来提交 Slurm 作业数组。"
echo "命令: sbatch 2_run_iqtree_array.sh"
echo "=================================================="

# 脚本结束 