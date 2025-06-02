#!/bin/bash

#SBATCH --job-name=codeml_branch # 中文注释：作业名称，branch models
#SBATCH --output=codeml_branch_%j.out # 中文注释：标准输出文件
#SBATCH --error=codeml_branch_%j.err  # 中文注释：标准错误文件
#SBATCH -N 1                             # 中文注释：请求1个节点
#SBATCH -n 1                             # 中文注释：在该节点上运行1个任务 (整个bash脚本视为1个任务)
#SBATCH --cpus-per-task=16               # 中文注释：为该任务请求CPU核心数，用于并行运行codeml
#SBATCH --mem=64G                        # 中文注释：请求内存

# BASE_NAMES 数组 (如果需要，可以用于标记树等，当前脚本不直接使用它进行PAML运行)
BASE_NAMES=(
"Acromyrmex_echinatior"
"Apis_cerana"
"Apis_dorsata"
"Apis_florea"
"Apis_laboriosa"
"Apis_mellifera"
"Atta_cephalotes"
"Atta_colombica"
"Bombus_affinis"
"Bombus_bifarius"
"Bombus_fervidus"
"Bombus_flavifrons"
"Bombus_huntii"
"Bombus_impatiens"
"Bombus_pascuorum"
"Bombus_pyrosoma"
"Bombus_terrestris"
"Bombus_vancouverensis_nearcticus"
"Bombus_vosnesenskii"
"Camponotus_floridanus"
"Cardiocondyla_obscurior"
"Cataglyphis_hispanica"
"Ceratina_calcarata"
"Colletes_gigas"
"Cyphomyrmex_costatus"
"Dinoponera_quadriceps"
"Drosophila_melanogaster"
"Dufourea_novaeangliae"
"Eufriesea_mexicana"
"Formica_exsecta"
"Frieseomelitta_varia"
"Habropoda_laboriosa"
"Harpegnathos_saltator"
"Hylaeus_anthracinus"
"Hylaeus_volcanicus"
"Linepithema_humile"
"Megachile_rotundata"
"Megalopta_genalis"
"Monomorium_pharaonis"
"Nomia_melanderi"
"Nylanderia_fulva"
"Odontomachus_brunneus"
"Ooceraea_biroi"
"Osmia_bicornis_bicornis"
"Osmia_lignaria"
"Pogonomyrmex_barbatus"
"Polistes_canadensis"
"Polistes_dominula"
"Polistes_fuscatus"
"Polyergus_mexicanus"
"Prorops_nasuta"
"Pseudomyrmex_gracilis"
"Solenopsis_invicta"
"Temnothorax_curvispinosus"
"Temnothorax_longispinosus"
"Temnothorax_nylanderi"
"Trachymyrmex_cornetzi"
"Trachymyrmex_septentrionalis"
"Trachymyrmex_zeteki"
"Vespa_crabro"
"Vespa_mandarinia"
"Vespa_velutina"
"Vespula_pensylvanica"
"Vespula_vulgaris"
"Vollenhovia_emeryi"
"Wasmannia_auropunctata"
)

# 硬编码路径变量
SEQ_ALN_DIR="DAS_aligned_codon_clipkit"
TREE_DIR_PATH="gene_trees_with_foreground" # 包含 *_codon.clipkit_marked.treefile 文件
OUTPUT_DIR="branch_model"      # 新的输出目录名
M0_CTL_TEMPLATE="branch_M0.ctl"          # M0 模型CTL模板
TWO_RATIO_CTL_TEMPLATE="branch_2ratio.ctl" # 2-ratio 模型CTL模板

echo "配置信息 (硬编码):"
echo "  序列比对目录: $SEQ_ALN_DIR"
echo "  树文件目录: $TREE_DIR_PATH"
echo "  输出目录: $OUTPUT_DIR"
echo "  M0 模型CTL模板: $M0_CTL_TEMPLATE"
echo "  2-ratio 模型CTL模板: $TWO_RATIO_CTL_TEMPLATE"

# 检查输入文件/目录是否存在
if [ ! -d "$SEQ_ALN_DIR" ]; then
  echo "错误: 序列比对目录 '$SEQ_ALN_DIR' 未找到。"
  exit 1
fi
if [ ! -d "$TREE_DIR_PATH" ]; then
  echo "错误: 树文件目录 '$TREE_DIR_PATH' 未找到或不是一个目录。"
  exit 1
fi
if [ ! -f "$M0_CTL_TEMPLATE" ]; then
  echo "错误: M0模型CTL模板 '$M0_CTL_TEMPLATE' 未找到。请确保它位于脚本执行的当前目录。"
  exit 1
fi
if [ ! -f "$TWO_RATIO_CTL_TEMPLATE" ]; then
  echo "错误: 2-ratio模型CTL模板 '$TWO_RATIO_CTL_TEMPLATE' 未找到。请确保它位于脚本执行的当前目录。"
  exit 1
fi

# 创建输出目录 (如果不存在)
mkdir -p "$OUTPUT_DIR"
echo "确保输出目录存在: $OUTPUT_DIR"

# 加载CTL模板内容
m0_ctl_template_content=$(cat "$M0_CTL_TEMPLATE")
two_ratio_ctl_template_content=$(cat "$TWO_RATIO_CTL_TEMPLATE")

# 定义CTL模板中的通用占位符
SEQFILE_PLACEHOLDER="<你的序列比对文件名.phy>"
TREEFILE_PLACEHOLDER_IN_CTL="/home/kosukesano/tools/for_paml/IQTREE_6sp/data/new_tree_IQTREE_ultrametric.nwk"

# 定义M0模型特定的占位符
M0_OUTFILE_PLACEHOLDER="<geneX_AE_branch_M0_results.txt>"

# 定义2-ratio模型特定的占位符
TWO_RATIO_OUTFILE_PLACEHOLDER="<geneX_AE_branch_2ratio_results.txt>"

# 设置并行任务数
job_count=0
max_jobs=1 # 默认至少1个
if [ -n "$SLURM_CPUS_PER_TASK" ] && [ "$SLURM_CPUS_PER_TASK" -gt 1 ]; then
    max_jobs=$((SLURM_CPUS_PER_TASK -1)) # 留一个核心给主脚本或其他开销
elif type nproc &>/dev/null && [ "$(nproc)" -gt 1 ]; then
    max_jobs=$(($(nproc) -1)) # 非SLURM环境下的备用方案
fi
if [ "$max_jobs" -lt 1 ]; then max_jobs=1; fi
echo "最大并行PAML任务数: $max_jobs"

# 处理目录中的每个序列比对文件
for seq_file_path in "$SEQ_ALN_DIR"/*_codon.clipkit.fasta; do
  if [ -f "$seq_file_path" ]; then
    seq_file_basename_full=$(basename "$seq_file_path")
    og_base_name="${seq_file_basename_full%.fasta}" #例如 OG0001155_codon.clipkit
    gene_name="$og_base_name"

    # 构建对应的树文件路径 (期望 *_codon.clipkit_marked.treefile)
    current_tree_file="$TREE_DIR_PATH/${og_base_name}_marked.treefile"

    if [ ! -f "$current_tree_file" ]; then
      echo "警告: 序列文件 '$seq_file_path' 对应的树文件 '$current_tree_file' 未找到。跳过此序列。"
      continue # 跳过当前序列文件
    fi

    abs_seq_file_path=$(realpath "$seq_file_path")
    abs_tree_file_path=$(realpath "$current_tree_file")

    # --- M0 模型 (one-ratio) ---
    m0_ctl_filename="${gene_name}_M0.ctl"
    m0_ctl_path="$OUTPUT_DIR/$m0_ctl_filename"
    m0_paml_outfile_path="$OUTPUT_DIR/${gene_name}_M0_paml_results.txt"
    m0_codeml_log_path="$OUTPUT_DIR/${gene_name}_M0.codeml.log"

    current_m0_ctl_content=$(echo "$m0_ctl_template_content" | \
        sed "s|$SEQFILE_PLACEHOLDER|$abs_seq_file_path|g" | \
        sed "s|$TREEFILE_PLACEHOLDER_IN_CTL|$abs_tree_file_path|g" | \
        sed "s|$M0_OUTFILE_PLACEHOLDER|$m0_paml_outfile_path|g")
    echo "$current_m0_ctl_content" > "$m0_ctl_path"

    echo "正在为 $gene_name 启动PAML (M0 模型)..."
    singularity exec -e -B /lustre10:/lustre10 /usr/local/biotools/p/paml:4.9--h779adbc_6 codeml "$m0_ctl_path" > "$m0_codeml_log_path" 2>&1 &
    job_count=$((job_count + 1))
    if [ "$job_count" -ge "$max_jobs" ]; then
        echo "达到最大并行任务数 ($max_jobs), 等待一个任务完成..."
        wait -n
        job_count=$((job_count - 1))
    fi

    # --- 2-ratio 模型 ---
    two_ratio_ctl_filename="${gene_name}_2ratio.ctl"
    two_ratio_ctl_path="$OUTPUT_DIR/$two_ratio_ctl_filename"
    two_ratio_paml_outfile_path="$OUTPUT_DIR/${gene_name}_2ratio_paml_results.txt"
    two_ratio_codeml_log_path="$OUTPUT_DIR/${gene_name}_2ratio.codeml.log"

    current_two_ratio_ctl_content=$(echo "$two_ratio_ctl_template_content" | \
        sed "s|$SEQFILE_PLACEHOLDER|$abs_seq_file_path|g" | \
        sed "s|$TREEFILE_PLACEHOLDER_IN_CTL|$abs_tree_file_path|g" | \
        sed "s|$TWO_RATIO_OUTFILE_PLACEHOLDER|$two_ratio_paml_outfile_path|g")
    echo "$current_two_ratio_ctl_content" > "$two_ratio_ctl_path"

    echo "正在为 $gene_name 启动PAML (2-ratio 模型)..."
    singularity exec -e -B /lustre10:/lustre10 /usr/local/biotools/p/paml:4.9--h779adbc_6 codeml "$two_ratio_ctl_path" > "$two_ratio_codeml_log_path" 2>&1 &
    job_count=$((job_count + 1))
    if [ "$job_count" -ge "$max_jobs" ]; then
        echo "达到最大并行任务数 ($max_jobs), 等待一个任务完成..."
        wait -n
        job_count=$((job_count - 1))
    fi

    echo "已为 $gene_name 提交PAML任务 (M0 模型 和 2-ratio 模型)。CTL文件: $m0_ctl_path, $two_ratio_ctl_path"
  else
    echo "警告: '$seq_file_path' 不是一个文件, 跳过。"
  fi
done

echo "所有PAML任务已提交。等待剩余任务完成..."
wait # 等待所有后台任务执行完毕
echo "所有PAML (branch models) 分析已完成。结果保存在: $OUTPUT_DIR" 