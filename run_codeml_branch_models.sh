#!/bin/bash

#SBATCH --job-name=paml_M0_Br_BsA_multiOmega # 中文注释：作业名称，M0, Branch (multi-omega), Branch-site A (multi-omega Alt & Null)
#SBATCH --output=paml_M0_Br_BsA_multiOmega_%j.out # 中文注释：标准输出文件
#SBATCH --error=paml_M0_Br_BsA_multiOmega_%j.err  # 中文注释：标准错误文件
#SBATCH -N 1                             # 中文注释：请求1个节点
#SBATCH -n 1                             # 中文注释：在该节点上运行1个任务 (整个bash脚本视为1个任务)
#SBATCH --cpus-per-task=16               # 中文注释：为该任务请求CPU核心数，用于并行运行codeml
#SBATCH --mem=64G                        # 中文注释：请求内存

# 定义项目基础路径
PROJECT_BASE_PATH="/home/yuhangjia/data/AlternativeSplicing/9_overlap_DAS_evo_rate_analysis"

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

# 硬编码路径变量 (使用绝对路径)
SEQ_ALN_DIR="$PROJECT_BASE_PATH/DAS_aligned_codon_clipkit"
TREE_DIR_PATH="$PROJECT_BASE_PATH/gene_trees_from_M0_remarked" # 使用M0提取并标记的树
OUTPUT_DIR_NAME="$PROJECT_BASE_PATH/paml_M0_Branch_BranchSiteA_multiOmega_on_M0trees" # 新的输出目录名

M0_CTL_TEMPLATE="$PROJECT_BASE_PATH/branch_M0.ctl" # M0 模型CTL模板

# 分支模型 (原2-ratio) CTL模板 - 不同初始omega
BRANCH_CTL_OMEGA05_TEMPLATE="$PROJECT_BASE_PATH/branch_2ratio_omega0.5.ctl"
BRANCH_CTL_OMEGA10_TEMPLATE="$PROJECT_BASE_PATH/branch_2ratio_omega1.0.ctl"
BRANCH_CTL_OMEGA15_TEMPLATE="$PROJECT_BASE_PATH/branch_2ratio_omega1.5.ctl"
BRANCH_CTL_OMEGA20_TEMPLATE="$PROJECT_BASE_PATH/branch_2ratio_omega2.0.ctl"

# 分支位点模型A - 备择模型 CTL模板 - 不同初始omega
BSA_ALT_CTL_OMEGA05_TEMPLATE="$PROJECT_BASE_PATH/bsA_alt_omega0.5.ctl"
BSA_ALT_CTL_OMEGA10_TEMPLATE="$PROJECT_BASE_PATH/bsA_alt_omega1.0.ctl"
BSA_ALT_CTL_OMEGA15_TEMPLATE="$PROJECT_BASE_PATH/bsA_alt_omega1.5.ctl"
BSA_ALT_CTL_OMEGA20_TEMPLATE="$PROJECT_BASE_PATH/bsA_alt_omega2.0.ctl"

BSA_NULL_CTL_TEMPLATE="$PROJECT_BASE_PATH/bsA_null.ctl" # 分支位点模型A - 空模型 CTL模板

# 获取输出目录的绝对路径
ABS_OUTPUT_DIR=$(realpath "$OUTPUT_DIR_NAME")

echo "配置信息 (硬编码):"
echo "  项目基础路径: $PROJECT_BASE_PATH"
echo "  序列比对目录: $SEQ_ALN_DIR"
echo "  树文件目录 (使用M0提取并标记的树): $TREE_DIR_PATH"
echo "  输出目录 (绝对): $ABS_OUTPUT_DIR"
echo "  M0 模型 CTL模板: $M0_CTL_TEMPLATE"
echo "  分支模型 CTL模板 (omega 0.5): $BRANCH_CTL_OMEGA05_TEMPLATE"
echo "  分支模型 CTL模板 (omega 1.0): $BRANCH_CTL_OMEGA10_TEMPLATE"
echo "  分支模型 CTL模板 (omega 1.5): $BRANCH_CTL_OMEGA15_TEMPLATE"
echo "  分支模型 CTL模板 (omega 2.0): $BRANCH_CTL_OMEGA20_TEMPLATE"
echo "  分支位点模型A - 备择 CTL模板 (omega 0.5): $BSA_ALT_CTL_OMEGA05_TEMPLATE"
echo "  分支位点模型A - 备择 CTL模板 (omega 1.0): $BSA_ALT_CTL_OMEGA10_TEMPLATE"
echo "  分支位点模型A - 备择 CTL模板 (omega 1.5): $BSA_ALT_CTL_OMEGA15_TEMPLATE"
echo "  分支位点模型A - 备择 CTL模板 (omega 2.0): $BSA_ALT_CTL_OMEGA20_TEMPLATE"
echo "  分支位点模型A - 空 CTL模板: $BSA_NULL_CTL_TEMPLATE"

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
  echo "错误: M0模型CTL模板 '$M0_CTL_TEMPLATE' 未找到。"
  exit 1
fi
if [ ! -f "$BRANCH_CTL_OMEGA05_TEMPLATE" ]; then
  echo "错误: 分支模型 (omega 0.5) CTL模板 '$BRANCH_CTL_OMEGA05_TEMPLATE' 未找到。"
  exit 1
fi
if [ ! -f "$BRANCH_CTL_OMEGA10_TEMPLATE" ]; then
  echo "错误: 分支模型 (omega 1.0) CTL模板 '$BRANCH_CTL_OMEGA10_TEMPLATE' 未找到。"
  exit 1
fi
if [ ! -f "$BRANCH_CTL_OMEGA15_TEMPLATE" ]; then
  echo "错误: 分支模型 (omega 1.5) CTL模板 '$BRANCH_CTL_OMEGA15_TEMPLATE' 未找到。"
  exit 1
fi
if [ ! -f "$BRANCH_CTL_OMEGA20_TEMPLATE" ]; then
  echo "错误: 分支模型 (omega 2.0) CTL模板 '$BRANCH_CTL_OMEGA20_TEMPLATE' 未找到。"
  exit 1
fi
if [ ! -f "$BSA_ALT_CTL_OMEGA05_TEMPLATE" ]; then
  echo "错误: 分支位点模型A - 备择 (omega 0.5) CTL模板 '$BSA_ALT_CTL_OMEGA05_TEMPLATE' 未找到。"
  exit 1
fi
if [ ! -f "$BSA_ALT_CTL_OMEGA10_TEMPLATE" ]; then
  echo "错误: 分支位点模型A - 备择 (omega 1.0) CTL模板 '$BSA_ALT_CTL_OMEGA10_TEMPLATE' 未找到。"
  exit 1
fi
if [ ! -f "$BSA_ALT_CTL_OMEGA15_TEMPLATE" ]; then
  echo "错误: 分支位点模型A - 备择 (omega 1.5) CTL模板 '$BSA_ALT_CTL_OMEGA15_TEMPLATE' 未找到。"
  exit 1
fi
if [ ! -f "$BSA_ALT_CTL_OMEGA20_TEMPLATE" ]; then
  echo "错误: 分支位点模型A - 备择 (omega 2.0) CTL模板 '$BSA_ALT_CTL_OMEGA20_TEMPLATE' 未找到。"
  exit 1
fi
if [ ! -f "$BSA_NULL_CTL_TEMPLATE" ]; then
  echo "错误: 分支位点模型A - 空 CTL模板 '$BSA_NULL_CTL_TEMPLATE' 未找到。"
  exit 1
fi

# 创建主输出目录 (如果不存在)
mkdir -p "$ABS_OUTPUT_DIR"
echo "确保主输出目录存在: $ABS_OUTPUT_DIR"

# 加载CTL模板内容
m0_ctl_template_content=$(cat "$M0_CTL_TEMPLATE")
branch_ctl_omega05_content=$(cat "$BRANCH_CTL_OMEGA05_TEMPLATE")
branch_ctl_omega10_content=$(cat "$BRANCH_CTL_OMEGA10_TEMPLATE")
branch_ctl_omega15_content=$(cat "$BRANCH_CTL_OMEGA15_TEMPLATE")
branch_ctl_omega20_content=$(cat "$BRANCH_CTL_OMEGA20_TEMPLATE")
bsa_alt_ctl_omega05_content=$(cat "$BSA_ALT_CTL_OMEGA05_TEMPLATE")
bsa_alt_ctl_omega10_content=$(cat "$BSA_ALT_CTL_OMEGA10_TEMPLATE")
bsa_alt_ctl_omega15_content=$(cat "$BSA_ALT_CTL_OMEGA15_TEMPLATE")
bsa_alt_ctl_omega20_content=$(cat "$BSA_ALT_CTL_OMEGA20_TEMPLATE")
bsa_null_ctl_template_content=$(cat "$BSA_NULL_CTL_TEMPLATE")

# 定义CTL模板中的通用占位符
SEQFILE_PLACEHOLDER="<你的序列比对文件名.phy>"
TREEFILE_PLACEHOLDER_IN_CTL="/home/kosukesano/tools/for_paml/IQTREE_6sp/data/new_tree_IQTREE_ultrametric.nwk" # 确保此占位符与CTL文件实际使用的一致

# 定义各模型输出文件占位符 (确保这些与您CTL模板中的占位符一致)
M0_OUTFILE_PLACEHOLDER="<geneX_AE_branch_M0_results.txt>"
BRANCH_OUTFILE_PLACEHOLDER_BASE="<geneX_AE_branch_2ratio_results.txt>" # 基础占位符，实际会在脚本中添加omega后缀
BSA_ALT_OUTFILE_PLACEHOLDER_BASE="<geneX_AE_bsA_alt_results.txt>"     # 基础占位符
BSA_NULL_OUTFILE_PLACEHOLDER="<geneX_AE_bsA_null_results.txt>"


# 设置并行任务数
job_count=0
max_jobs=1 # 默认至少1个
if [ -n "$SLURM_CPUS_PER_TASK" ] && [ "$SLURM_CPUS_PER_TASK" -gt 1 ]; then
    max_jobs=$((SLURM_CPUS_PER_TASK -1))
elif type nproc &>/dev/null && [ "$(nproc)" -gt 1 ]; then
    max_jobs=$(($(nproc) -1))
fi
if [ "$max_jobs" -lt 1 ]; then max_jobs=1; fi
echo "最大并行PAML任务数: $max_jobs"

initial_omegas_branch=("0.5" "1.0" "1.5" "2.0")
initial_omegas_bsa_alt=("0.5" "1.0" "1.5" "2.0")

# 处理目录中的每个序列比对文件
for seq_file_path in "$SEQ_ALN_DIR"/*_codon.clipkit.fasta; do
  if [ -f "$seq_file_path" ]; then
    seq_file_basename_full=$(basename "$seq_file_path")
    og_base_name="${seq_file_basename_full%.fasta}"
    gene_name="$og_base_name"

    current_tree_file="$TREE_DIR_PATH/${og_base_name}_from_M0_marked.treefile"

    if [ ! -f "$current_tree_file" ]; then
      echo "警告: 序列文件 '$seq_file_path' 对应的树文件 '$current_tree_file' 未找到。跳过此序列。"
      continue
    fi

    abs_seq_file_path=$(realpath "$seq_file_path")
    abs_tree_file_path=$(realpath "$current_tree_file")

    # --- M0 模型 ---
    m0_workdir="$ABS_OUTPUT_DIR/${gene_name}_M0_work"
    mkdir -p "$m0_workdir"
    m0_ctl_filename_basename="${gene_name}_M0.ctl"
    m0_ctl_path_in_workdir="$m0_workdir/$m0_ctl_filename_basename"
    m0_paml_outfile_path_abs="$ABS_OUTPUT_DIR/${gene_name}_M0_results.txt"
    m0_codeml_log_path_abs="$ABS_OUTPUT_DIR/${gene_name}_M0.codeml.log"

    current_m0_ctl_content=$(echo "$m0_ctl_template_content" | \
        sed "s|$SEQFILE_PLACEHOLDER|$abs_seq_file_path|g" | \
        sed "s#$TREEFILE_PLACEHOLDER_IN_CTL#$abs_tree_file_path#g" | \
        sed "s|$M0_OUTFILE_PLACEHOLDER|$m0_paml_outfile_path_abs|g")
    echo "$current_m0_ctl_content" > "$m0_ctl_path_in_workdir"

    echo "正在为 $gene_name 启动PAML (M0 模型)... 工作目录: $m0_workdir ; 树文件: $abs_tree_file_path"
    (cd "$m0_workdir" && singularity exec -e -B /lustre10:/lustre10 /usr/local/biotools/p/paml:4.9--h779adbc_6 codeml "$m0_ctl_filename_basename" > "$m0_codeml_log_path_abs" 2>&1) &
    job_count=$((job_count + 1))
    if [ "$job_count" -ge "$max_jobs" ]; then
        echo "达到最大并行任务数 ($max_jobs), 等待一个任务完成..."
        wait -n
        job_count=$((job_count - 1))
    fi

    # --- 分支模型 (不同初始omega) ---
    for omega_val in "${initial_omegas_branch[@]}"; do
        omega_suffix_fs=$(echo "$omega_val" | sed 's/\./p/') # for filenames, e.g., 0.5 -> 0p5
        current_branch_ctl_template_var_name="branch_ctl_omega${omega_suffix_fs}_content"
        current_branch_ctl_content="${!current_branch_ctl_template_var_name}" # Indirect expansion

        branch_workdir="$ABS_OUTPUT_DIR/${gene_name}_branch_omega${omega_suffix_fs}_work"
        mkdir -p "$branch_workdir"
        branch_ctl_filename_basename="${gene_name}_branch_omega${omega_suffix_fs}.ctl"
        branch_ctl_path_in_workdir="$branch_workdir/$branch_ctl_filename_basename"
        branch_paml_outfile_path_abs="$ABS_OUTPUT_DIR/${gene_name}_branch_omega${omega_suffix_fs}_results.txt"
        branch_codeml_log_path_abs="$ABS_OUTPUT_DIR/${gene_name}_branch_omega${omega_suffix_fs}.codeml.log"
        
        # 替换CTL中的输出文件名占位符。注意：这里假设原始CTL模板中的输出文件名占位符是BRANCH_OUTFILE_PLACEHOLDER_BASE
        # 并且不包含omega值。我们在这里将实际的、包含omega的输出文件名替换掉那个基础占位符。
        temp_ctl_content_branch=$(echo "$current_branch_ctl_content" | \
            sed "s|$SEQFILE_PLACEHOLDER|$abs_seq_file_path|g" | \
            sed "s#$TREEFILE_PLACEHOLDER_IN_CTL#$abs_tree_file_path#g" | \
            sed "s|$BRANCH_OUTFILE_PLACEHOLDER_BASE|$branch_paml_outfile_path_abs|g")
        echo "$temp_ctl_content_branch" > "$branch_ctl_path_in_workdir"

        echo "正在为 $gene_name 启动PAML (分支模型, omega_init=$omega_val)... 工作目录: $branch_workdir"
        (cd "$branch_workdir" && singularity exec -e -B /lustre10:/lustre10 /usr/local/biotools/p/paml:4.9--h779adbc_6 codeml "$branch_ctl_filename_basename" > "$branch_codeml_log_path_abs" 2>&1) &
        job_count=$((job_count + 1))
        if [ "$job_count" -ge "$max_jobs" ]; then
            echo "达到最大并行任务数 ($max_jobs), 等待一个任务完成..."
            wait -n
            job_count=$((job_count - 1))
        fi
    done

    # --- 分支位点模型 A - 备择模型 (不同初始omega) ---
    for omega_val in "${initial_omegas_bsa_alt[@]}"; do
        omega_suffix_fs=$(echo "$omega_val" | sed 's/\./p/') # e.g., 0.5 -> 0p5
        current_bsa_alt_ctl_template_var_name="bsa_alt_ctl_omega${omega_suffix_fs}_content"
        current_bsa_alt_ctl_content="${!current_bsa_alt_ctl_template_var_name}" # Indirect expansion

        bsa_alt_workdir="$ABS_OUTPUT_DIR/${gene_name}_bsA_alt_omega${omega_suffix_fs}_work"
        mkdir -p "$bsa_alt_workdir"
        bsa_alt_ctl_filename_basename="${gene_name}_bsA_alt_omega${omega_suffix_fs}.ctl"
        bsa_alt_ctl_path_in_workdir="$bsa_alt_workdir/$bsa_alt_ctl_filename_basename"
        bsa_alt_paml_outfile_path_abs="$ABS_OUTPUT_DIR/${gene_name}_bsA_alt_omega${omega_suffix_fs}_results.txt"
        bsa_alt_codeml_log_path_abs="$ABS_OUTPUT_DIR/${gene_name}_bsA_alt_omega${omega_suffix_fs}.codeml.log"

        temp_ctl_content_bsa_alt=$(echo "$current_bsa_alt_ctl_content" | \
            sed "s|$SEQFILE_PLACEHOLDER|$abs_seq_file_path|g" | \
            sed "s#$TREEFILE_PLACEHOLDER_IN_CTL#$abs_tree_file_path#g" | \
            sed "s|$BSA_ALT_OUTFILE_PLACEHOLDER_BASE|$bsa_alt_paml_outfile_path_abs|g")
        echo "$temp_ctl_content_bsa_alt" > "$bsa_alt_ctl_path_in_workdir"

        echo "正在为 $gene_name 启动PAML (分支位点模型 A - 备择, omega_init=$omega_val)... 工作目录: $bsa_alt_workdir"
        (cd "$bsa_alt_workdir" && singularity exec -e -B /lustre10:/lustre10 /usr/local/biotools/p/paml:4.9--h779adbc_6 codeml "$bsa_alt_ctl_filename_basename" > "$bsa_alt_codeml_log_path_abs" 2>&1) &
        job_count=$((job_count + 1))
        if [ "$job_count" -ge "$max_jobs" ]; then
            echo "达到最大并行任务数 ($max_jobs), 等待一个任务完成..."
            wait -n
            job_count=$((job_count - 1))
        fi
    done

    # --- 分支位点模型 A - 空模型 ---
    bsa_null_workdir="$ABS_OUTPUT_DIR/${gene_name}_bsA_null_work"
    mkdir -p "$bsa_null_workdir"
    bsa_null_ctl_filename_basename="${gene_name}_bsA_null.ctl"
    bsa_null_ctl_path_in_workdir="$bsa_null_workdir/$bsa_null_ctl_filename_basename"
    bsa_null_paml_outfile_path_abs="$ABS_OUTPUT_DIR/${gene_name}_bsA_null_results.txt"
    bsa_null_codeml_log_path_abs="$ABS_OUTPUT_DIR/${gene_name}_bsA_null.codeml.log"

    current_bsa_null_ctl_content=$(echo "$bsa_null_ctl_template_content" | \
        sed "s|$SEQFILE_PLACEHOLDER|$abs_seq_file_path|g" | \
        sed "s#$TREEFILE_PLACEHOLDER_IN_CTL#$abs_tree_file_path#g" | \
        sed "s|$BSA_NULL_OUTFILE_PLACEHOLDER|$bsa_null_paml_outfile_path_abs|g")
    echo "$current_bsa_null_ctl_content" > "$bsa_null_ctl_path_in_workdir"

    echo "正在为 $gene_name 启动PAML (分支位点模型 A - 空)... 工作目录: $bsa_null_workdir"
    (cd "$bsa_null_workdir" && singularity exec -e -B /lustre10:/lustre10 /usr/local/biotools/p/paml:4.9--h779adbc_6 codeml "$bsa_null_ctl_filename_basename" > "$bsa_null_codeml_log_path_abs" 2>&1) &
    job_count=$((job_count + 1))
    if [ "$job_count" -ge "$max_jobs" ]; then
        echo "达到最大并行任务数 ($max_jobs), 等待一个任务完成..."
        wait -n
        job_count=$((job_count - 1))
    fi

    echo "已为 $gene_name 提交所有PAML任务 (M0, Branch-multiOmega, BsA-multiOmegaAlt, BsA-Null)。CTL文件在各自的 _work 目录中。"
  else
    echo "警告: '$seq_file_path' 不是一个文件, 跳过。"
  fi
done

echo "所有PAML任务已提交。等待剩余任务完成..."
wait # 等待所有后台任务执行完毕
echo "所有PAML (M0, Branch-multiOmega, BranchSiteA-multiOmegaAlt & Null) 分析已完成。结果保存在: $ABS_OUTPUT_DIR (主要结果和日志), 各自的 _work 子目录包含辅助文件。" 