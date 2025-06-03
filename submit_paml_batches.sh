#!/bin/bash

# --- Configuration for the main submission script ---
PROJECT_BASE_PATH="/home/yuhangjia/data/AlternativeSplicing/evo_rate_test_RNA_splicing_term" # 需要与原脚本一致或作为参数传入
SEQ_ALN_DIR_PARAM="${PROJECT_BASE_PATH}/GO_0008380_msa_codon_clipkit"
TREE_DIR_PATH_PARAM="${PROJECT_BASE_PATH}/gene_trees_from_M0_remarked"
MAIN_OUTPUT_DIR_PARAM="${PROJECT_BASE_PATH}/paml_M0_Branch_BranchSiteA_batched_output" # 新的主输出目录

# CTL模板路径 (与原脚本一致)
M0_CTL_TEMPLATE_PARAM="${PROJECT_BASE_PATH}/branch_M0.ctl"
BRANCH_CTL_OMEGA0P5_TEMPLATE_PARAM="${PROJECT_BASE_PATH}/branch_2ratio_omega0.5.ctl"
BRANCH_CTL_OMEGA1P0_TEMPLATE_PARAM="${PROJECT_BASE_PATH}/branch_2ratio_omega1.0.ctl"
BRANCH_CTL_OMEGA1P5_TEMPLATE_PARAM="${PROJECT_BASE_PATH}/branch_2ratio_omega1.5.ctl"
BRANCH_CTL_OMEGA2P0_TEMPLATE_PARAM="${PROJECT_BASE_PATH}/branch_2ratio_omega2.0.ctl"
BSA_ALT_CTL_OMEGA0P5_TEMPLATE_PARAM="${PROJECT_BASE_PATH}/bsA_alt_omega0.5.ctl"
BSA_ALT_CTL_OMEGA1P0_TEMPLATE_PARAM="${PROJECT_BASE_PATH}/bsA_alt_omega1.0.ctl"
BSA_ALT_CTL_OMEGA1P5_TEMPLATE_PARAM="${PROJECT_BASE_PATH}/bsA_alt_omega1.5.ctl"
BSA_ALT_CTL_OMEGA2P0_TEMPLATE_PARAM="${PROJECT_BASE_PATH}/bsA_alt_omega2.0.ctl"
BSA_NULL_CTL_TEMPLATE_PARAM="${PROJECT_BASE_PATH}/bsA_null.ctl"

GENES_PER_BATCH=12 # 用户指定：每个 SLURM 作业处理多少个基因
SLURM_CPUS_PER_SUB_JOB=128 # 用户指定：每个子 SLURM 作业请求多少 CPU 核心 (用于内部 PAML 并行)
SLURM_MEM_PER_SUB_JOB="128G" # 用户指定：每个子 SLURM 作业请求多少内存
SLURM_TIME_LIMIT="01:00:00" # 用户指定：每个子 SLURM 作业的时间限制
SUBMISSION_COMMAND="sbatch" # 作业提交命令

BATCH_SCRIPT_DIR="${MAIN_OUTPUT_DIR_PARAM}/batch_scripts" # 存储生成的批处理脚本
SINGULARITY_IMAGE_PATH="/usr/local/biotools/p/paml:4.9--h779adbc_6" # PAML Singularity 镜像路径
SINGULARITY_BIND_PATHS="-B /lustre10:/lustre10" # Singularity 绑定路径, 根据需要修改, e.g. "-B /path1:/path1 -B /path2:/path2"

# --- End Configuration ---

echo "主提交脚本开始执行..."
echo "  项目基础路径: ${PROJECT_BASE_PATH}"
echo "  序列比对目录: ${SEQ_ALN_DIR_PARAM}"
echo "  树文件目录: ${TREE_DIR_PATH_PARAM}"
echo "  主输出目录: ${MAIN_OUTPUT_DIR_PARAM}"
echo "  每个批次处理基因数: ${GENES_PER_BATCH}"
echo "  每个子作业CPU核心数: ${SLURM_CPUS_PER_SUB_JOB}"
echo "  每个子作业内存: ${SLURM_MEM_PER_SUB_JOB}"
echo "  每个子作业时间限制: ${SLURM_TIME_LIMIT}"
echo "  PAML Singularity 镜像: ${SINGULARITY_IMAGE_PATH}"
echo "  Singularity 绑定: ${SINGULARITY_BIND_PATHS}"

# 检查输入文件/目录是否存在 (CTL模板等)
check_path_exists() {
    local path_to_check="$1"
    local type_flag="$2" # "-f" for file, "-d" for directory
    local entity_type="路径"
    if [ "$type_flag" == "-f" ]; then entity_type="文件"; fi
    if [ "$type_flag" == "-d" ]; then entity_type="目录"; fi

    if [ ! "$type_flag" "$path_to_check" ]; then
        echo "错误: ${entity_type} '${path_to_check}' 未找到。"
        exit 1
    fi
}

check_path_exists "${SEQ_ALN_DIR_PARAM}" "-d"
check_path_exists "${TREE_DIR_PATH_PARAM}" "-d"
check_path_exists "${M0_CTL_TEMPLATE_PARAM}" "-f"
check_path_exists "${BRANCH_CTL_OMEGA0P5_TEMPLATE_PARAM}" "-f"
check_path_exists "${BRANCH_CTL_OMEGA1P0_TEMPLATE_PARAM}" "-f"
check_path_exists "${BRANCH_CTL_OMEGA1P5_TEMPLATE_PARAM}" "-f"
check_path_exists "${BRANCH_CTL_OMEGA2P0_TEMPLATE_PARAM}" "-f"
check_path_exists "${BSA_ALT_CTL_OMEGA0P5_TEMPLATE_PARAM}" "-f"
check_path_exists "${BSA_ALT_CTL_OMEGA1P0_TEMPLATE_PARAM}" "-f"
check_path_exists "${BSA_ALT_CTL_OMEGA1P5_TEMPLATE_PARAM}" "-f"
check_path_exists "${BSA_ALT_CTL_OMEGA2P0_TEMPLATE_PARAM}" "-f"
check_path_exists "${BSA_NULL_CTL_TEMPLATE_PARAM}" "-f"


# 创建主输出目录和批处理脚本存储目录
mkdir -p "${MAIN_OUTPUT_DIR_PARAM}"
mkdir -p "${BATCH_SCRIPT_DIR}"
echo "主输出目录已创建: ${MAIN_OUTPUT_DIR_PARAM}"
echo "批处理脚本将存储在: ${BATCH_SCRIPT_DIR}"

# 获取所有序列比对文件路径
mapfile -t all_seq_files < <(find "${SEQ_ALN_DIR_PARAM}" -name "*_codon.clipkit.fasta" -type f | sort)

num_total_genes=${#all_seq_files[@]}
if [ "$num_total_genes" -eq 0 ]; then
    echo "错误: 在 ${SEQ_ALN_DIR_PARAM} 中没有找到序列文件 (*_codon.clipkit.fasta)。"
    exit 1
fi
echo "总共找到 ${num_total_genes} 个基因进行处理。"

num_batches=$(( (num_total_genes + GENES_PER_BATCH - 1) / GENES_PER_BATCH ))
echo "这些基因将被分成 ${num_batches} 个批次进行提交。"

for (( batch_num=0; batch_num<num_batches; batch_num++ )); do
    start_index=$(( batch_num * GENES_PER_BATCH ))
    
    current_batch_gene_files_abs_paths=()
    for (( i=0; i<GENES_PER_BATCH; i++ )); do
        gene_index=$(( start_index + i ))
        if [ "$gene_index" -lt "$num_total_genes" ]; then
            # 获取文件的绝对路径
            abs_path_gene_file=$(realpath "${all_seq_files[$gene_index]}")
            current_batch_gene_files_abs_paths+=("${abs_path_gene_file}")
        else
            break 
        fi
    done

    if [ ${#current_batch_gene_files_abs_paths[@]} -eq 0 ]; then
        echo "批次 $((batch_num + 1)) 为空，跳过。"
        continue
    fi
    
    batch_id=$((batch_num + 1))
    batch_job_name="paml_batch_${batch_id}"
    batch_script_path="${BATCH_SCRIPT_DIR}/${batch_job_name}.sh"
    # 主输出目录对所有子作业是共享的，使用绝对路径
    batch_output_dir_abs=$(realpath "${MAIN_OUTPUT_DIR_PARAM}")

    echo "正在为批次 ${batch_id} 生成脚本: ${batch_script_path}"
    echo "  该批次包含 ${#current_batch_gene_files_abs_paths[@]} 个基因。"

    # 生成子作业脚本
    # 使用 cat 和 EOF 来创建脚本文件，确保变量正确展开或转义
    # 注意：在EOF块内部，需要转义那些我们不希望立即展开的变量（例如 $SLURM_JOB_ID）
    # 而那些我们希望从主脚本代入的变量（例如 ${batch_job_name}）则不需要转义。
    cat > "${batch_script_path}" <<EOF
#!/bin/bash
#SBATCH --job-name=${batch_job_name}
#SBATCH --output=${batch_output_dir_abs}/${batch_job_name}_%j.out
#SBATCH --error=${batch_output_dir_abs}/${batch_job_name}_%j.err
#SBATCH -N 1
#SBATCH -n 1 
#SBATCH --cpus-per-task=${SLURM_CPUS_PER_SUB_JOB}
#SBATCH --mem=${SLURM_MEM_PER_SUB_JOB}
#SBATCH --time=${SLURM_TIME_LIMIT}
#SBATCH --partition=short

echo "PAML批处理作业 ${batch_job_name} (Job ID: \$SLURM_JOB_ID) 开始于 \$(date)"
echo "  处理的基因数: ${#current_batch_gene_files_abs_paths[@]}"
echo "  输出目录: ${batch_output_dir_abs}"
echo "  CPU核心数: \$SLURM_CPUS_PER_TASK"
echo "  内存: \${SLURM_MEM_PER_TASK:-${SLURM_MEM_PER_SUB_JOB}}" # Use SLURM_MEM_PER_TASK if set, else fallback
echo "  时间限制: ${SLURM_TIME_LIMIT}"
echo "  PAML Singularity 镜像: ${SINGULARITY_IMAGE_PATH}"
echo "  Singularity 绑定: ${SINGULARITY_BIND_PATHS}"


# --- 配置 (从主脚本继承或硬编码，因为它们对于所有子作业是相同的) ---
# 注意: 这些路径应该是子作业脚本执行时可访问的绝对路径
TREE_DIR_PATH="${TREE_DIR_PATH_PARAM}" # 绝对路径
ABS_OUTPUT_DIR="${batch_output_dir_abs}" # 绝对路径

M0_CTL_TEMPLATE="${M0_CTL_TEMPLATE_PARAM}"
BRANCH_CTL_OMEGA0P5_TEMPLATE="${BRANCH_CTL_OMEGA0P5_TEMPLATE_PARAM}"
BRANCH_CTL_OMEGA1P0_TEMPLATE="${BRANCH_CTL_OMEGA1P0_TEMPLATE_PARAM}"
BRANCH_CTL_OMEGA1P5_TEMPLATE="${BRANCH_CTL_OMEGA1P5_TEMPLATE_PARAM}"
BRANCH_CTL_OMEGA2P0_TEMPLATE="${BRANCH_CTL_OMEGA2P0_TEMPLATE_PARAM}"
BSA_ALT_CTL_OMEGA0P5_TEMPLATE="${BSA_ALT_CTL_OMEGA0P5_TEMPLATE_PARAM}"
BSA_ALT_CTL_OMEGA1P0_TEMPLATE="${BSA_ALT_CTL_OMEGA1P0_TEMPLATE_PARAM}"
BSA_ALT_CTL_OMEGA1P5_TEMPLATE="${BSA_ALT_CTL_OMEGA1P5_TEMPLATE_PARAM}"
BSA_ALT_CTL_OMEGA2P0_TEMPLATE="${BSA_ALT_CTL_OMEGA2P0_TEMPLATE_PARAM}"
BSA_NULL_CTL_TEMPLATE="${BSA_NULL_CTL_TEMPLATE_PARAM}"

# 确保输出目录存在 (通常主脚本已创建，但以防万一)
mkdir -p "\${ABS_OUTPUT_DIR}"

# 加载CTL模板内容 (与原脚本相同)
m0_ctl_template_content=\$(cat "\${M0_CTL_TEMPLATE}")
branch_ctl_omega0p5_content=\$(cat "\${BRANCH_CTL_OMEGA0P5_TEMPLATE}")
branch_ctl_omega1p0_content=\$(cat "\${BRANCH_CTL_OMEGA1P0_TEMPLATE}")
branch_ctl_omega1p5_content=\$(cat "\${BRANCH_CTL_OMEGA1P5_TEMPLATE}")
branch_ctl_omega2p0_content=\$(cat "\${BRANCH_CTL_OMEGA2P0_TEMPLATE}")
bsa_alt_ctl_omega0p5_content=\$(cat "\${BSA_ALT_CTL_OMEGA0P5_TEMPLATE}")
bsa_alt_ctl_omega1p0_content=\$(cat "\${BSA_ALT_CTL_OMEGA1P0_TEMPLATE}")
bsa_alt_ctl_omega1p5_content=\$(cat "\${BSA_ALT_CTL_OMEGA1P5_TEMPLATE}")
bsa_alt_ctl_omega2p0_content=\$(cat "\${BSA_ALT_CTL_OMEGA2P0_TEMPLATE}")
bsa_null_ctl_template_content=\$(cat "\${BSA_NULL_CTL_TEMPLATE}")

# 定义CTL模板中的通用占位符 (与原脚本相同)
SEQFILE_PLACEHOLDER="<你的序列比对文件名.phy>"
TREEFILE_PLACEHOLDER_IN_CTL="/home/kosukesano/tools/for_paml/IQTREE_6sp/data/new_tree_IQTREE_ultrametric.nwk" # 确保此占位符与CTL文件实际使用的一致
M0_OUTFILE_PLACEHOLDER="<geneX_AE_branch_M0_results.txt>"
BRANCH_OUTFILE_PLACEHOLDER_BASE="<geneX_AE_branch_2ratio_results.txt>" 
BSA_ALT_OUTFILE_PLACEHOLDER_BASE="<geneX_AE_bsA_alt_results.txt>"     
BSA_NULL_OUTFILE_PLACEHOLDER="<geneX_AE_bsA_null_results.txt>"

# PAML 并行设置 (与原脚本相同)
job_count=0
max_jobs=1 
if [ -n "\$SLURM_CPUS_PER_TASK" ] && [ "\$SLURM_CPUS_PER_TASK" -gt 0 ]; then # Check for > 0
    max_jobs=\$SLURM_CPUS_PER_TASK # 使用所有分配的核心
elif type nproc &>/dev/null && [ "\$(nproc)" -gt 0 ]; then # Check for > 0
    max_jobs=\$(nproc)
fi
if [ "\$max_jobs" -lt 1 ]; then max_jobs=1; fi
echo "  子作业 PAML 并行数 (max_jobs): \${max_jobs}"

initial_omegas_branch=("0.5" "1.0" "1.5" "2.0")
initial_omegas_bsa_alt=("0.5" "1.0" "1.5" "2.0")

# 基因文件列表 (直接嵌入 - 路径已是绝对路径)
batch_gene_files_abs_paths=(
EOF

    # 将当前批次的绝对路径写入子脚本
    for gene_f_path in "${current_batch_gene_files_abs_paths[@]}"; do
        # 引号是为了处理路径中可能存在的空格等特殊字符
        echo "    \"${gene_f_path}\"" >> "${batch_script_path}"
    done

    # 继续生成子作业脚本的剩余部分
    cat >> "${batch_script_path}" <<EOF
)

echo "  将要处理的第一个基因文件示例 (绝对路径): \${batch_gene_files_abs_paths[0]}"

# 处理此批次中的每个序列比对文件
for abs_seq_file_path in "\${batch_gene_files_abs_paths[@]}"; do
  if [ -f "\${abs_seq_file_path}" ]; then
    seq_file_basename_full=\$(basename "\${abs_seq_file_path}")
    og_base_name="\${seq_file_basename_full%.fasta}" # e.g., OGxxxx_codon.clipkit
    gene_name="\${og_base_name}" # gene_name = OGxxxx_codon.clipkit

    # 构建树文件路径 (使用绝对路径的 TREE_DIR_PATH)
    # og_base_name 是 xxx_codon.clipkit, 树文件名是 xxx_codon.clipkit_from_M0_marked.treefile
    current_tree_file="\${TREE_DIR_PATH}/\${og_base_name}_from_M0_marked.treefile"
    
    if [ ! -f "\${current_tree_file}" ]; then
      echo "警告 (基因 \${gene_name}): 对应的树文件 '\${current_tree_file}' 未找到。跳过此基因。"
      continue
    fi
    abs_tree_file_path=\$(realpath "\${current_tree_file}") # 树文件也用绝对路径

    echo "" #空行分隔不同基因处理日志
    echo "--- [\$(date +'%Y-%m-%d %H:%M:%S')] 开始处理基因: \${gene_name} (源: \${abs_seq_file_path}) ---"
    
    # --- M0 模型 ---
    m0_workdir="\${ABS_OUTPUT_DIR}/\${gene_name}_M0_work"
    mkdir -p "\${m0_workdir}"
    m0_ctl_filename_basename="\${gene_name}_M0.ctl"
    m0_ctl_path_in_workdir="\${m0_workdir}/\${m0_ctl_filename_basename}"
    m0_paml_outfile_path_abs="\${ABS_OUTPUT_DIR}/\${gene_name}_M0_results.txt" # 绝对路径
    m0_codeml_log_path_abs="\${ABS_OUTPUT_DIR}/\${gene_name}_M0.codeml.log"   # 绝对路径

    current_m0_ctl_content=\$(echo "\${m0_ctl_template_content}" | \\
        sed "s|\${SEQFILE_PLACEHOLDER}|\${abs_seq_file_path}|g" | \\
        sed "s#\${TREEFILE_PLACEHOLDER_IN_CTL}#\${abs_tree_file_path}#g" | \\
        sed "s|\${M0_OUTFILE_PLACEHOLDER}|\${m0_paml_outfile_path_abs}|g")
    echo "\${current_m0_ctl_content}" > "\${m0_ctl_path_in_workdir}"

    echo "  [\${gene_name} M0] 启动PAML... CTL: \${m0_ctl_path_in_workdir} 输出到: \${m0_paml_outfile_path_abs}"
    (cd "\${m0_workdir}" && singularity exec -e ${SINGULARITY_BIND_PATHS} ${SINGULARITY_IMAGE_PATH} codeml "\${m0_ctl_filename_basename}" > "\${m0_codeml_log_path_abs}" 2>&1) &
    job_count=\$((job_count + 1))
    if [ "\${job_count}" -ge "\${max_jobs}" ]; then
        echo "    M0: 达到最大并行 (\${max_jobs}), 等待一个PAML任务完成..."
        wait -n
        job_count=\$((job_count - 1))
    fi

    # --- 分支模型 (不同初始omega) ---
    for omega_val in "\${initial_omegas_branch[@]}"; do
        omega_suffix_fs=\$(echo "\${omega_val}" | sed 's/\\./p/') 
        current_branch_ctl_template_var_name="branch_ctl_omega\${omega_suffix_fs}_content"
        # 使用 bash 间接变量展开
        printf -v current_branch_ctl_content "%s" "\${!current_branch_ctl_template_var_name}"


        branch_workdir="\${ABS_OUTPUT_DIR}/\${gene_name}_branch_omega\${omega_suffix_fs}_work"
        mkdir -p "\${branch_workdir}"
        branch_ctl_filename_basename="\${gene_name}_branch_omega\${omega_suffix_fs}.ctl"
        branch_ctl_path_in_workdir="\${branch_workdir}/\${branch_ctl_filename_basename}"
        branch_paml_outfile_path_abs="\${ABS_OUTPUT_DIR}/\${gene_name}_branch_omega\${omega_suffix_fs}_results.txt"
        branch_codeml_log_path_abs="\${ABS_OUTPUT_DIR}/\${gene_name}_branch_omega\${omega_suffix_fs}.codeml.log"
        
        temp_ctl_content_branch=\$(echo "\${current_branch_ctl_content}" | \\
            sed "s|\${SEQFILE_PLACEHOLDER}|\${abs_seq_file_path}|g" | \\
            sed "s#\${TREEFILE_PLACEHOLDER_IN_CTL}#\${abs_tree_file_path}#g" | \\
            sed "s|\${BRANCH_OUTFILE_PLACEHOLDER_BASE}|\${branch_paml_outfile_path_abs}|g")
        echo "\${temp_ctl_content_branch}" > "\${branch_ctl_path_in_workdir}"

        echo "  [\${gene_name} Branch omega_init=\${omega_val}] 启动PAML... 输出到: \${branch_paml_outfile_path_abs}"
        (cd "\${branch_workdir}" && singularity exec -e ${SINGULARITY_BIND_PATHS} ${SINGULARITY_IMAGE_PATH} codeml "\${branch_ctl_filename_basename}" > "\${branch_codeml_log_path_abs}" 2>&1) &
        job_count=\$((job_count + 1))
        if [ "\${job_count}" -ge "\${max_jobs}" ]; then
            echo "    Branch \${omega_val}: 达到最大并行 (\${max_jobs}), 等待一个PAML任务完成..."
            wait -n
            job_count=\$((job_count - 1))
        fi
    done

    # --- 分支位点模型 A - 备择模型 (不同初始omega) ---
    for omega_val in "\${initial_omegas_bsa_alt[@]}"; do
        omega_suffix_fs=\$(echo "\${omega_val}" | sed 's/\\./p/')
        current_bsa_alt_ctl_template_var_name="bsa_alt_ctl_omega\${omega_suffix_fs}_content"
        printf -v current_bsa_alt_ctl_content "%s" "\${!current_bsa_alt_ctl_template_var_name}"


        bsa_alt_workdir="\${ABS_OUTPUT_DIR}/\${gene_name}_bsA_alt_omega\${omega_suffix_fs}_work"
        mkdir -p "\${bsa_alt_workdir}"
        bsa_alt_ctl_filename_basename="\${gene_name}_bsA_alt_omega\${omega_suffix_fs}.ctl"
        bsa_alt_ctl_path_in_workdir="\${bsa_alt_workdir}/\${bsa_alt_ctl_filename_basename}"
        bsa_alt_paml_outfile_path_abs="\${ABS_OUTPUT_DIR}/\${gene_name}_bsA_alt_omega\${omega_suffix_fs}_results.txt"
        bsa_alt_codeml_log_path_abs="\${ABS_OUTPUT_DIR}/\${gene_name}_bsA_alt_omega\${omega_suffix_fs}.codeml.log"

        temp_ctl_content_bsa_alt=\$(echo "\${current_bsa_alt_ctl_content}" | \\
            sed "s|\${SEQFILE_PLACEHOLDER}|\${abs_seq_file_path}|g" | \\
            sed "s#\${TREEFILE_PLACEHOLDER_IN_CTL}#\${abs_tree_file_path}#g" | \\
            sed "s|\${BSA_ALT_OUTFILE_PLACEHOLDER_BASE}|\${bsa_alt_paml_outfile_path_abs}|g")
        echo "\${temp_ctl_content_bsa_alt}" > "\${bsa_alt_ctl_path_in_workdir}"

        echo "  [\${gene_name} BsA-Alt omega_init=\${omega_val}] 启动PAML... 输出到: \${bsa_alt_paml_outfile_path_abs}"
        (cd "\${bsa_alt_workdir}" && singularity exec -e ${SINGULARITY_BIND_PATHS} ${SINGULARITY_IMAGE_PATH} codeml "\${bsa_alt_ctl_filename_basename}" > "\${bsa_alt_codeml_log_path_abs}" 2>&1) &
        job_count=\$((job_count + 1))
        if [ "\${job_count}" -ge "\${max_jobs}" ]; then
            echo "    BsA-Alt \${omega_val}: 达到最大并行 (\${max_jobs}), 等待一个PAML任务完成..."
            wait -n
            job_count=\$((job_count - 1))
        fi
    done

    # --- 分支位点模型 A - 空模型 ---
    bsa_null_workdir="\${ABS_OUTPUT_DIR}/\${gene_name}_bsA_null_work"
    mkdir -p "\${bsa_null_workdir}"
    bsa_null_ctl_filename_basename="\${gene_name}_bsA_null.ctl"
    bsa_null_ctl_path_in_workdir="\${bsa_null_workdir}/\${bsa_null_ctl_filename_basename}"
    bsa_null_paml_outfile_path_abs="\${ABS_OUTPUT_DIR}/\${gene_name}_bsA_null_results.txt"
    bsa_null_codeml_log_path_abs="\${ABS_OUTPUT_DIR}/\${gene_name}_bsA_null.codeml.log"

    current_bsa_null_ctl_content=\$(echo "\${bsa_null_ctl_template_content}" | \\
        sed "s|\${SEQFILE_PLACEHOLDER}|\${abs_seq_file_path}|g" | \\
        sed "s#\${TREEFILE_PLACEHOLDER_IN_CTL}#\${abs_tree_file_path}#g" | \\
        sed "s|\${BSA_NULL_OUTFILE_PLACEHOLDER}|\${bsa_null_paml_outfile_path_abs}|g")
    echo "\${current_bsa_null_ctl_content}" > "\${bsa_null_ctl_path_in_workdir}"

    echo "  [\${gene_name} BsA-Null] 启动PAML... 输出到: \${bsa_null_paml_outfile_path_abs}"
    (cd "\${bsa_null_workdir}" && singularity exec -e ${SINGULARITY_BIND_PATHS} ${SINGULARITY_IMAGE_PATH} codeml "\${bsa_null_ctl_filename_basename}" > "\${bsa_null_codeml_log_path_abs}" 2>&1) &
    job_count=\$((job_count + 1))
    if [ "\${job_count}" -ge "\${max_jobs}" ]; then
        echo "    BsA-Null: 达到最大并行 (\${max_jobs}), 等待一个PAML任务完成..."
        wait -n
        job_count=\$((job_count - 1))
    fi
    echo "--- [\$(date +'%Y-%m-%d %H:%M:%S')] 完成基因 \${gene_name} 的所有PAML任务提交 ---"
  else
    echo "警告 (子作业 \$SLURM_JOB_ID): 基因文件 '\${abs_seq_file_path}' 不是一个文件或未找到, 跳过。"
  fi
done

echo ""
echo "批处理作业 ${batch_job_name}: 所有基因的PAML任务已提交。等待所有剩余的PAML后台进程完成..."
wait
echo "批处理作业 ${batch_job_name} (Job ID: \$SLURM_JOB_ID) 所有PAML进程执行完毕于 \$(date)。"
EOF
# End of cat >> for the rest of the sub-job script

    chmod +x "${batch_script_path}"
    echo "  子作业脚本已创建并设为可执行: ${batch_script_path}"
    
    # 提交子作业脚本
    echo "  正在提交批次 ${batch_id} (${batch_script_path})..."
    # Store sbatch output, if any
    submission_output=$(${SUBMISSION_COMMAND} "${batch_script_path}" 2>&1)
    submission_status=$?
    
    if [ ${submission_status} -eq 0 ]; then
        echo "    成功提交: ${submission_output}"
    else
        echo "    错误: 提交失败 (状态码: ${submission_status})。输出:"
        echo "${submission_output}"
    fi
    
    # 可选：在提交下一个批次前稍作停顿，以避免瞬间给SLURM控制器过大压力
    # sleep 2
done

echo ""
echo "所有批次的作业已提交到SLURM。"
echo "请监控 '${MAIN_OUTPUT_DIR_PARAM}' 目录下的输出和日志文件。"
echo "生成的批处理脚本位于 '${BATCH_SCRIPT_DIR}'。"
echo "主提交脚本执行完毕。" 