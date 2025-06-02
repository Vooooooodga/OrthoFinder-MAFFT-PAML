#!/bin/bash

#SBATCH --job-name=codeml_batch # 中文注释：作业名称，修改为批处理codeml
#SBATCH --output=codeml_batch_%j.out # 中文注释：标准输出文件
#SBATCH --error=codeml_batch_%j.err  # 中文注释：标准错误文件
#SBATCH -N 1                             # 中文注释：请求1个节点
#SBATCH -n 1                             # 中文注释：在该节点上运行1个任务 (整个bash脚本视为1个任务)
#SBATCH --cpus-per-task=32               # 中文注释：为该任务请求CPU核心数，用于并行运行codeml
#SBATCH --mem=128G                        # 中文注释：请求内存

# 函数：显示用法信息
usage() {
  echo "用法: $0 -s <seq_aln_dir> -t <tree_file_path> -o <output_dir> -a <alt_ctl_template> -n <null_ctl_template>"
  echo "  -s: 序列比对文件目录 (包含 *_maffted.fna 文件)"
  echo "  -t: Newick 格式树文件的完整路径"
  echo "  -o: 输出目录 (用于存放生成的CTL文件和PAML结果)"
  echo "  -a: 备择模型CTL模板文件路径 (例如 bsA_alt.ctl)"
  echo "  -n: 零假设模型CTL模板文件路径 (例如 bsA_null.ctl)"
  exit 1
}

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

# 初始化变量
SEQ_ALN_DIR=""
TREE_FILE_PATH=""
OUTPUT_DIR=""
ALT_CTL_TEMPLATE=""
NULL_CTL_TEMPLATE=""

# 解析命令行参数
while getopts ":s:t:o:a:n:" opt; do
  case $opt in
    s) SEQ_ALN_DIR="$OPTARG" ;;
    t) TREE_FILE_PATH="$OPTARG" ;;
    o) OUTPUT_DIR="$OPTARG" ;;
    a) ALT_CTL_TEMPLATE="$OPTARG" ;;
    n) NULL_CTL_TEMPLATE="$OPTARG" ;;
    \?) echo "无效选项: -$OPTARG" >&2; usage ;;
    :) echo "选项 -$OPTARG 需要一个参数." >&2; usage ;;
  esac
done

# 检查是否提供了所有必需参数
if [ -z "$SEQ_ALN_DIR" ] || [ -z "$TREE_FILE_PATH" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$ALT_CTL_TEMPLATE" ] || [ -z "$NULL_CTL_TEMPLATE" ]; then
  echo "错误: 缺少一个或多个必需参数。"
  usage
fi

# 检查输入文件/目录是否存在
if [ ! -d "$SEQ_ALN_DIR" ]; then
  echo "错误: 序列比对目录 '$SEQ_ALN_DIR' 未找到。"
  exit 1
fi
if [ ! -f "$TREE_FILE_PATH" ]; then
  echo "错误: 树文件 '$TREE_FILE_PATH' 未找到。"
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

# 创建输出目录 (如果不存在)
mkdir -p "$OUTPUT_DIR"
echo "输出目录: $OUTPUT_DIR"

# 加载CTL模板内容
alt_ctl_template_content=$(cat "$ALT_CTL_TEMPLATE")
null_ctl_template_content=$(cat "$NULL_CTL_TEMPLATE")

# 定义CTL模板中的占位符 (根据您提供的bsA_*.ctl文件)
SEQFILE_PLACEHOLDER="<你的序列比对文件名.phy>"
# 注意：TREEFILE_PLACEHOLDER_IN_CTL 是模板中硬编码的树文件路径，需要被替换
TREEFILE_PLACEHOLDER_IN_CTL="/home/kosukesano/tools/for_paml/IQTREE_6sp/data/new_tree_IQTREE_ultrametric.nwk"
ALT_OUTFILE_PLACEHOLDER="<geneX_AE_bsA_alt_results.txt>"
NULL_OUTFILE_PLACEHOLDER="<geneX_AE_bsA_null_results.txt>"

# 设置并行任务数
job_count=0
max_jobs=1 # 默认至少1个
if [ -n "$SLURM_CPUS_PER_TASK" ] && [ "$SLURM_CPUS_PER_TASK" -gt 1 ]; then
    max_jobs=$((SLURM_CPUS_PER_TASK - 1)) # 留一个核心给主脚本或其他开销
elif type nproc &>/dev/null && [ "$(nproc)" -gt 1 ]; then
    max_jobs=$(($(nproc) - 1)) # 非SLURM环境下的备用方案
fi
if [ "$max_jobs" -lt 1 ]; then max_jobs=1; fi
echo "最大并行PAML任务数: $max_jobs"

# 处理目录中的每个序列比对文件
for seq_file_path in "$SEQ_ALN_DIR"/*_maffted.fna; do
  if [ -f "$seq_file_path" ]; then
    seq_file_basename=$(basename "$seq_file_path")
    gene_name="${seq_file_basename%_maffted.fna}" # 获取基因名, 例如 OG0000000

    # 获取文件的绝对路径，确保codeml能找到它们
    # realpath可能需要coreutils包。如果系统没有，请确保路径正确或使用其他方法获取绝对路径。
    abs_seq_file_path=$(realpath "$seq_file_path")
    abs_tree_file_path=$(realpath "$TREE_FILE_PATH")

    # --- 备择模型 (Alternative Model) ---
    alt_ctl_filename="${gene_name}_alt.ctl"
    alt_ctl_path="$OUTPUT_DIR/$alt_ctl_filename"
    # PAML输出文件的路径 (将在CTL文件中指定)
    alt_paml_outfile_path="$OUTPUT_DIR/${gene_name}_alt_paml_results.txt"
    # codeml命令本身的标准输出/错误重定向文件
    alt_codeml_log_path="$OUTPUT_DIR/${gene_name}_alt.codeml.log"

    # 生成备择模型的CTL文件内容
    # 使用'|'作为sed的分隔符，以避免路径中的'/'导致问题
    current_alt_ctl_content=$(echo "$alt_ctl_template_content" | \
        sed "s|$SEQFILE_PLACEHOLDER|$abs_seq_file_path|g" | \
        sed "s|$TREEFILE_PLACEHOLDER_IN_CTL|$abs_tree_file_path|g" | \
        sed "s|$ALT_OUTFILE_PLACEHOLDER|$alt_paml_outfile_path|g")
    echo "$current_alt_ctl_content" > "$alt_ctl_path"

    echo "正在为 $gene_name 启动PAML (备择模型)..."
    # 后台运行codeml，并重定向其标准输出和错误到日志文件
    singularity exec -e /usr/local/biotools/p/paml:4.9--h779adbc_6 codeml "$alt_ctl_path" > "$alt_codeml_log_path" 2>&1 &
    job_count=$((job_count + 1))
    if [ "$job_count" -ge "$max_jobs" ]; then
        echo "达到最大并行任务数 ($max_jobs), 等待一个任务完成..."
        wait -n # 等待任意一个后台任务结束
        job_count=$((job_count - 1))
    fi

    # --- 零假设模型 (Null Model) ---
    null_ctl_filename="${gene_name}_null.ctl"
    null_ctl_path="$OUTPUT_DIR/$null_ctl_filename"
    null_paml_outfile_path="$OUTPUT_DIR/${gene_name}_null_paml_results.txt"
    null_codeml_log_path="$OUTPUT_DIR/${gene_name}_null.codeml.log"

    # 生成零假设模型的CTL文件内容
    current_null_ctl_content=$(echo "$null_ctl_template_content" | \
        sed "s|$SEQFILE_PLACEHOLDER|$abs_seq_file_path|g" | \
        sed "s|$TREEFILE_PLACEHOLDER_IN_CTL|$abs_tree_file_path|g" | \
        sed "s|$NULL_OUTFILE_PLACEHOLDER|$null_paml_outfile_path|g")
    echo "$current_null_ctl_content" > "$null_ctl_path"

    echo "正在为 $gene_name 启动PAML (零假设模型)..."
    singularity exec -e /usr/local/biotools/p/paml:4.9--h779adbc_6 codeml "$null_ctl_path" > "$null_codeml_log_path" 2>&1 &
    job_count=$((job_count + 1))
    if [ "$job_count" -ge "$max_jobs" ]; then
        echo "达到最大并行任务数 ($max_jobs), 等待一个任务完成..."
        wait -n
        job_count=$((job_count - 1))
    fi
    echo "已为 $gene_name 提交PAML任务 (备择模型和零假设模型)。CTL文件: $alt_ctl_path, $null_ctl_path"
  else
    echo "警告: '$seq_file_path' 不是一个文件, 跳过。"
  fi
done

echo "所有PAML任务已提交。等待剩余任务完成..."
wait # 等待所有后台任务执行完毕
echo "所有PAML分析已完成。结果保存在: $OUTPUT_DIR"