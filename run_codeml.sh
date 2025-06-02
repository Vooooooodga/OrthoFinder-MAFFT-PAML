#!/bin/bash

#SBATCH --job-name=codeml_batch # 中文注释：作业名称，修改为批处理codeml
#SBATCH --output=codeml_batch_%j.out # 中文注释：标准输出文件
#SBATCH --error=codeml_batch_%j.err  # 中文注释：标准错误文件
#SBATCH -N 1                             # 中文注释：请求1个节点
#SBATCH -n 1                             # 中文注释：在该节点上运行1个任务 (整个bash脚本视为1个任务)
#SBATCH --cpus-per-task=32               # 中文注释：为该任务请求CPU核心数，用于并行运行codeml
#SBATCH --mem=128G                        # 中文注释：请求内存

# 定义项目基础路径
PROJECT_BASE_PATH="/home/yuhangjia/data/AlternativeSplicing/9_overlap_DAS_evo_rate_analysis"

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
"Bombus_vancouverensis_nearcticus" # 包含亚种名
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
"Osmia_bicornis_bicornis" # 包含亚种名
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

# 硬编码路径变量 (使用绝对路径)
SEQ_ALN_DIR="$PROJECT_BASE_PATH/DAS_aligned_codon_clipkit"
TREE_DIR_PATH="$PROJECT_BASE_PATH/gene_trees_with_foreground"
OUTPUT_DIR_NAME="$PROJECT_BASE_PATH/bsA" # 输出目录的名称
ALT_CTL_TEMPLATE="$PROJECT_BASE_PATH/bsA_alt.ctl"   # 备择模型CTL模板
NULL_CTL_TEMPLATE="$PROJECT_BASE_PATH/bsA_null.ctl" # 零假设模型CTL模板

# 获取输出目录的绝对路径
ABS_OUTPUT_DIR=$(realpath "$OUTPUT_DIR_NAME")

echo "配置信息 (硬编码):"
echo "  项目基础路径: $PROJECT_BASE_PATH"
echo "  序列比对目录: $SEQ_ALN_DIR"
echo "  树文件目录: $TREE_DIR_PATH"
echo "  输出目录 (绝对): $ABS_OUTPUT_DIR"
echo "  备择模型CTL模板: $ALT_CTL_TEMPLATE"
echo "  零假设模型CTL模板: $NULL_CTL_TEMPLATE"

# 检查输入文件/目录是否存在
if [ ! -d "$SEQ_ALN_DIR" ]; then
  echo "错误: 序列比对目录 '$SEQ_ALN_DIR' 未找到。"
  exit 1
fi
if [ ! -d "$TREE_DIR_PATH" ]; then
  echo "错误: 树文件目录 '$TREE_DIR_PATH' 未找到或不是一个目录。"
  exit 1
fi
if [ ! -f "$ALT_CTL_TEMPLATE" ]; then
  echo "错误: 备择模型CTL模板 '$ALT_CTL_TEMPLATE' 未找到。"
  exit 1
fi
if [ ! -f "$NULL_CTL_TEMPLATE" ]; then
  echo "错误: 零假设模型CTL模板 '$NULL_CTL_TEMPLATE' 未找到。"
  exit 1
fi

# 创建主输出目录 (如果不存在)
mkdir -p "$ABS_OUTPUT_DIR"
echo "确保主输出目录存在: $ABS_OUTPUT_DIR"

# 加载CTL模板内容
alt_ctl_template_content=$(cat "$ALT_CTL_TEMPLATE")
null_ctl_template_content=$(cat "$NULL_CTL_TEMPLATE")

# 定义CTL模板中的占位符 (根据您提供的bsA_*.ctl文件)
SEQFILE_PLACEHOLDER="<你的序列比对文件名.phy>"
TREEFILE_PLACEHOLDER_IN_CTL="/home/kosukesano/tools/for_paml/IQTREE_6sp/data/new_tree_IQTREE_ultrametric.nwk"
ALT_OUTFILE_PLACEHOLDER="<geneX_AE_bsA_alt_results.txt>"
NULL_OUTFILE_PLACEHOLDER="<geneX_AE_bsA_null_results.txt>"

# 设置并行任务数
job_count=0
max_jobs=1 # 默认至少1个
if [ -n "$SLURM_CPUS_PER_TASK" ] && [ "$SLURM_CPUS_PER_TASK" -gt 1 ]; then
    max_jobs=$((SLURM_CPUS_PER_TASK - 1))
elif type nproc &>/dev/null && [ "$(nproc)" -gt 1 ]; then
    max_jobs=$(($(nproc) - 1))
fi
if [ "$max_jobs" -lt 1 ]; then max_jobs=1; fi
echo "最大并行PAML任务数: $max_jobs"

# 处理目录中的每个序列比对文件
for seq_file_path in "$SEQ_ALN_DIR"/*_codon.clipkit.fasta; do
  if [ -f "$seq_file_path" ]; then
    seq_file_basename_full=$(basename "$seq_file_path")
    og_base_name="${seq_file_basename_full%.fasta}"
    gene_name="$og_base_name"

    current_tree_file="$TREE_DIR_PATH/${og_base_name}_marked.treefile"

    if [ ! -f "$current_tree_file" ]; then
      echo "警告: 序列文件 '$seq_file_path' 对应的树文件 '$current_tree_file' 未找到。跳过此序列。"
      continue
    fi

    abs_seq_file_path=$(realpath "$seq_file_path")
    abs_tree_file_path=$(realpath "$current_tree_file")

    # --- 备择模型 (Alternative Model) ---
    alt_workdir="$ABS_OUTPUT_DIR/${gene_name}_alt_work"
    mkdir -p "$alt_workdir"
    alt_ctl_filename_basename="${gene_name}_alt.ctl"
    alt_ctl_path_in_workdir="$alt_workdir/$alt_ctl_filename_basename"

    alt_paml_outfile_path_abs="$ABS_OUTPUT_DIR/${gene_name}_alt_paml_results.txt"
    alt_codeml_log_path_abs="$ABS_OUTPUT_DIR/${gene_name}_alt.codeml.log"

    current_alt_ctl_content=$(echo "$alt_ctl_template_content" | \
        sed "s|$SEQFILE_PLACEHOLDER|$abs_seq_file_path|g" | \
        sed "s|$TREEFILE_PLACEHOLDER_IN_CTL|$abs_tree_file_path|g" | \
        sed "s|$ALT_OUTFILE_PLACEHOLDER|$alt_paml_outfile_path_abs|g")
    echo "$current_alt_ctl_content" > "$alt_ctl_path_in_workdir"

    echo "正在为 $gene_name 启动PAML (备择模型)... 工作目录: $alt_workdir"
    (cd "$alt_workdir" && singularity exec -e -B /lustre10:/lustre10 /usr/local/biotools/p/paml:4.9--h779adbc_6 codeml "$alt_ctl_filename_basename" > "$alt_codeml_log_path_abs" 2>&1) &
    job_count=$((job_count + 1))
    if [ "$job_count" -ge "$max_jobs" ]; then
        wait -n
        job_count=$((job_count - 1))
    fi

    # --- 零假设模型 (Null Model) ---
    null_workdir="$ABS_OUTPUT_DIR/${gene_name}_null_work"
    mkdir -p "$null_workdir"
    null_ctl_filename_basename="${gene_name}_null.ctl"
    null_ctl_path_in_workdir="$null_workdir/$null_ctl_filename_basename"

    null_paml_outfile_path_abs="$ABS_OUTPUT_DIR/${gene_name}_null_paml_results.txt"
    null_codeml_log_path_abs="$ABS_OUTPUT_DIR/${gene_name}_null.codeml.log"

    current_null_ctl_content=$(echo "$null_ctl_template_content" | \
        sed "s|$SEQFILE_PLACEHOLDER|$abs_seq_file_path|g" | \
        sed "s|$TREEFILE_PLACEHOLDER_IN_CTL|$abs_tree_file_path|g" | \
        sed "s|$NULL_OUTFILE_PLACEHOLDER|$null_paml_outfile_path_abs|g")
    echo "$current_null_ctl_content" > "$null_ctl_path_in_workdir"

    echo "正在为 $gene_name 启动PAML (零假设模型)... 工作目录: $null_workdir"
    (cd "$null_workdir" && singularity exec -e -B /lustre10:/lustre10 /usr/local/biotools/p/paml:4.9--h779adbc_6 codeml "$null_ctl_filename_basename" > "$null_codeml_log_path_abs" 2>&1) &
    job_count=$((job_count + 1))
    if [ "$job_count" -ge "$max_jobs" ]; then
        wait -n
        job_count=$((job_count - 1))
    fi
    echo "已为 $gene_name 提交PAML任务 (备择模型和零假设模型)。CTL文件在各自的 _work 目录中。"
  else
    echo "警告: '$seq_file_path' 不是一个文件, 跳过。"
  fi
done

echo "所有PAML任务已提交。等待剩余任务完成..."
wait
echo "所有PAML (branch-site models) 分析已完成。结果保存在: $ABS_OUTPUT_DIR (主要结果和日志), 各自的 _work 子目录包含辅助文件。"