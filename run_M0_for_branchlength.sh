#!/bin/bash

#SBATCH --job-name=paml_M0_brlen_est # 作业名称: PAML M0 分支长度估计
#SBATCH --output=paml_M0_brlen_est_%j.out # 标准输出文件
#SBATCH --error=paml_M0_brlen_est_%j.err  # 标准错误文件
#SBATCH -N 1                             # 请求1个节点
#SBATCH -n 1                             # 在该节点上运行1个任务
#SBATCH --cpus-per-task=16               # 为该任务请求CPU核心数
#SBATCH --mem=32G                        # 请求内存

# 定义项目基础路径 (请根据您的实际情况修改)
PROJECT_BASE_PATH="/home/yuhangjia/data/AlternativeSplicing/9_overlap_DAS_evo_rate_analysis" 

# 硬编码路径变量
SEQ_ALN_DIR="$PROJECT_BASE_PATH/DAS_aligned_codon_clipkit_for_paml" # 序列比对文件目录
TREE_DIR_PATH="$PROJECT_BASE_PATH/gene_trees_with_foreground" # 包含原始树文件的目录 (M0会重新估计分支长度，但需要拓扑)
OUTPUT_DIR_NAME="$PROJECT_BASE_PATH/paml_M0_branch_length_estimation" # 输出目录名
M0_BRANCHLENGTH_CTL_TEMPLATE="./M0_for_branchlength_estimation.ctl" # 新的M0 CTL模板 (假设在当前脚本同目录下)

# 获取输出目录的绝对路径
ABS_OUTPUT_DIR=$(realpath "$OUTPUT_DIR_NAME")

echo "配置信息:"
echo "  项目基础路径: $PROJECT_BASE_PATH"
echo "  序列比对目录: $SEQ_ALN_DIR"
echo "  树文件目录: $TREE_DIR_PATH"
echo "  输出目录 (绝对): $ABS_OUTPUT_DIR"
echo "  M0分支长度估计 CTL模板: $M0_BRANCHLENGTH_CTL_TEMPLATE"

# 检查输入文件/目录是否存在
if [ ! -d "$SEQ_ALN_DIR" ]; then
  echo "错误: 序列比对目录 '$SEQ_ALN_DIR' 未找到。"
  exit 1
fi
if [ ! -d "$TREE_DIR_PATH" ]; then
  echo "错误: 树文件目录 '$TREE_DIR_PATH' 未找到或不是一个目录。"
  exit 1
fi
if [ ! -f "$M0_BRANCHLENGTH_CTL_TEMPLATE" ]; then
  echo "错误: M0分支长度估计CTL模板 '$M0_BRANCHLENGTH_CTL_TEMPLATE' 未找到。"
  exit 1
fi

# 创建主输出目录 (如果不存在)
mkdir -p "$ABS_OUTPUT_DIR"
echo "确保主输出目录存在: $ABS_OUTPUT_DIR"

# 加载CTL模板内容
m0_brlen_ctl_template_content=$(cat "$M0_BRANCHLENGTH_CTL_TEMPLATE")

# 定义CTL模板中的通用占位符 (与 M0_for_branchlength_estimation.ctl 中的占位符一致)
SEQFILE_PLACEHOLDER="<你的序列比对文件名.phy>"
TREEFILE_PLACEHOLDER_IN_CTL="/home/kosukesano/tools/for_paml/IQTREE_6sp/data/new_tree_IQTREE_ultrametric.nwk" # 确保此占位符与CTL文件实际使用的一致
M0_OUTFILE_PLACEHOLDER="<gene_name>_M0_for_branchlength_estimation_results.txt" # 确保此占位符与CTL文件实际使用的一致


# 设置并行任务数
job_count=0
max_jobs=1 # 默认至少1个
if [ -n "$SLURM_CPUS_PER_TASK" ] && [ "$SLURM_CPUS_PER_TASK" -gt 1 ]; then
    max_jobs=$((SLURM_CPUS_PER_TASK ))
elif type nproc &>/dev/null && [ "$(nproc)" -gt 1 ]; then
    max_jobs=$(($(nproc) -1))
fi
if [ "$max_jobs" -lt 1 ]; then max_jobs=1; fi
echo "最大并行PAML任务数: $max_jobs"

# 处理目录中的每个序列比对文件
for seq_file_path in "$SEQ_ALN_DIR"/*_codon.clipkit.fasta; do
  if [ -f "$seq_file_path" ]; then
    seq_file_basename_full=$(basename "$seq_file_path")
    # 假设基因名可以从序列文件名中提取，例如 OG0000001_codon.clipkit.fasta -> OG0000001
    gene_name="${seq_file_basename_full%_codon.clipkit.fasta}" 

    # 假设树文件名与序列文件名前缀一致，但后缀不同
    # 例如，如果序列文件是 OG0000001_codon.clipkit.fasta，对应的树文件可能是 OG0000001_codon.clipkit.fasta_from_M0_marked.treefile
    # 或者，您可能需要调整下面的逻辑以匹配您的树文件命名约定
    current_tree_file="$TREE_DIR_PATH/${seq_file_basename_full}_from_M0_marked.treefile" # 这是基于您原脚本的命名
    # 如果您的树文件命名更简单，如 ${gene_name}.treefile，请相应修改
    # current_tree_file="$TREE_DIR_PATH/${gene_name}.treefile" 

    if [ ! -f "$current_tree_file" ]; then
      echo "警告: 序列文件 '$seq_file_path' 对应的树文件 '$current_tree_file' 未找到。跳过此基因 '$gene_name'。"
      continue
    fi

    abs_seq_file_path=$(realpath "$seq_file_path")
    abs_tree_file_path=$(realpath "$current_tree_file")

    # --- M0 模型 (用于分支长度估计) ---
    m0_workdir="$ABS_OUTPUT_DIR/${gene_name}_M0_brlen_work"
    mkdir -p "$m0_workdir"
    m0_ctl_filename_basename="${gene_name}_M0_brlen.ctl"
    m0_ctl_path_in_workdir="$m0_workdir/$m0_ctl_filename_basename"
    m0_paml_outfile_path_abs="$ABS_OUTPUT_DIR/${gene_name}_M0_brlen_results.txt" # 与CTL中的outfile对应
    m0_codeml_log_path_abs="$ABS_OUTPUT_DIR/${gene_name}_M0_brlen.codeml.log"

    current_m0_ctl_content=$(echo "$m0_brlen_ctl_template_content" | \
        sed "s|$SEQFILE_PLACEHOLDER|$abs_seq_file_path|g" | \
        sed "s#$TREEFILE_PLACEHOLDER_IN_CTL#$abs_tree_file_path#g" | \
        sed "s|$M0_OUTFILE_PLACEHOLDER|$m0_paml_outfile_path_abs|g")
    echo "$current_m0_ctl_content" > "$m0_ctl_path_in_workdir"

    echo "正在为 $gene_name 启动PAML (M0 用于分支长度估计)... 工作目录: $m0_workdir ; 树文件: $abs_tree_file_path"
    # 注意：这里的 singularity 命令假设您的环境已配置好。
    # 您可能需要调整 -B 挂载路径以匹配您的系统。
    (cd "$m0_workdir" && singularity exec -e -B /lustre10:/lustre10 /usr/local/biotools/p/paml:4.9--h779adbc_6 codeml "$m0_ctl_filename_basename" > "$m0_codeml_log_path_abs" 2>&1) &
    
    job_count=$((job_count + 1))
    if [ "$job_count" -ge "$max_jobs" ]; then
        echo "达到最大并行任务数 ($max_jobs), 等待一个任务完成..."
        wait -n
        job_count=$((job_count - 1))
    fi
    
    echo "已为 $gene_name 提交PAML M0 (分支长度估计) 任务。CTL文件在 $m0_workdir 目录中。"
  else
    echo "警告: '$seq_file_path' 不是一个文件, 跳过。"
  fi
done

echo "所有PAML M0 (分支长度估计) 任务已提交。等待剩余任务完成..."
wait # 等待所有后台任务执行完毕
echo "所有PAML M0 (分支长度估计) 分析已完成。结果保存在: $ABS_OUTPUT_DIR (主要结果和日志), 各自的 _work 子目录包含辅助文件。"

# 脚本结束 