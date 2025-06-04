#!/bin/bash

#SBATCH --job-name=hyphy_relax_PARENT # 主脚本的作业名
#SBATCH --partition=your_partition      # 主脚本运行的分区 (可能与子作业不同)
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1             # 主脚本本身不需要很多CPU
#SBATCH --mem=2G                    # 主脚本本身不需要很多内存
#SBATCH --time=02:00:00               # 主脚本提交作业的时间，根据基因数量调整
#SBATCH --output=hyphy_relax_parent_%j.out
#SBATCH --error=hyphy_relax_parent_%j.err

# 加载 singularity 模块 (如果需要)
# module load singularity # 取消注释并根据您的集群环境修改

# --- 全局定义 ---
# TODO: 用户需要根据实际情况修改以下路径
MSA_DIR="/home/yuhangjia/data/AlternativeSplicing/hyphy_relax"
TREE_DIR="/home/yuhangjia/data/AlternativeSplicing/hyphy_relax"
MAIN_OUTPUT_DIR_RELAX="/home/yuhangjia/data/AlternativeSplicing/hyphy_relax/results" # 主输出目录 for RELAX
HYPHY_IMAGE="/usr/local/biotools/h/hyphy:2.5.65--he91c24d_0" # HyPhy singularity 镜像路径
SINGULARITY_BIND_OPTS="-B /lustre10:/lustre10" # 根据需要修改 Singularity 绑定路径

# 新的目录，用于存放修改后的树文件 (由父脚本统一处理)
MODIFIED_TREE_DIR_RELAX="${MAIN_OUTPUT_DIR_RELAX}/hyphy_modified_trees_relax"
# 新的目录，用于存放生成的子 sbatch 脚本
SUB_SBATCH_SCRIPT_DIR_RELAX="${MAIN_OUTPUT_DIR_RELAX}/sub_sbatch_scripts_relax"
# 子作业的输出将直接写入 MAIN_OUTPUT_DIR_RELAX/RELAX_json (json本身) 和 MAIN_OUTPUT_DIR_RELAX/RELAX_slurm_logs
RELAX_JSON_OUTPUT_DIR="${MAIN_OUTPUT_DIR_RELAX}/RELAX_json"
SLURM_LOGS_DIR_RELAX="${MAIN_OUTPUT_DIR_RELAX}/RELAX_slurm_logs"

# --- 创建所需目录 ---
echo "Creating directories for RELAX analysis..."
mkdir -p "${MODIFIED_TREE_DIR_RELAX}"
mkdir -p "${SUB_SBATCH_SCRIPT_DIR_RELAX}"
mkdir -p "${RELAX_JSON_OUTPUT_DIR}"
mkdir -p "${SLURM_LOGS_DIR_RELAX}"

echo "Starting to generate and submit HyPhy RELAX sub-jobs..."

# --- 主循环：遍历MSA文件，预处理树，生成并提交子脚本 ---
# DEBUG: List all files matching the pattern in MSA_DIR
echo "DEBUG: Files found in ${MSA_DIR} matching *_codon.clipkit.fasta:"
ls -1 "${MSA_DIR}"/*_codon.clipkit.fasta 2>/dev/null | tee /dev/stderr | wc -l
echo "DEBUG: ---- End of file list ----"

for msa_file_abs_path in "${MSA_DIR}"/*_codon.clipkit.fasta; do
    echo "DEBUG: Current msa_file_abs_path variable is: [${msa_file_abs_path}]"
    if [ -f "${msa_file_abs_path}" ]; then
        base_name=$(basename "${msa_file_abs_path}")
        gene_id=${base_name%_codon.clipkit.fasta}

        original_tree_file_abs_path="${TREE_DIR}/${gene_id}_from_M0_marked.treefile"
        # 输出的JSON文件路径 (子脚本将使用此路径)
        output_json_abs_path_relax="${RELAX_JSON_OUTPUT_DIR}/${gene_id}.RELAX.json"
        # 修改后的树文件路径 (子脚本将使用此路径)
        modified_tree_file_abs_path_relax="${MODIFIED_TREE_DIR_RELAX}/${gene_id}_from_M0_marked_relax.treefile"
        
        # 子 sbatch 脚本的路径
        sub_script_path_relax="${SUB_SBATCH_SCRIPT_DIR_RELAX}/sub_relax_${gene_id}.sh"

        # 1. 检查原始 tree 文件是否存在
        if [ ! -f "${original_tree_file_abs_path}" ]; then
            echo "WARNING: Original tree file not found for ${gene_id}: ${original_tree_file_abs_path}. Skipping gene."
            continue
        fi

        # 2. 预处理树文件 (由父脚本完成)
        #    HyPhy RELAX 需要树文件中标记 {test} 分支。
        #    此脚本将原始树文件中的 #1 标记 (代表前景分支) 替换为 {test}。
        #    假设未在树中用 {test} 标记的其他所有分支将被 HyPhy RELAX 自动视为参考 (reference) 分支。
        echo "Preprocessing tree file for ${gene_id} (RELAX): Replacing #1 with {test} (assuming other branches are reference)"
        sed 's/#1/{test}/g' "${original_tree_file_abs_path}" > "${modified_tree_file_abs_path_relax}"
        
        if [ ! -s "${modified_tree_file_abs_path_relax}" ]; then
            echo "ERROR: Failed to preprocess tree file for ${gene_id} (RELAX) or modified tree is empty. Original: ${original_tree_file_abs_path}. Check for the #1 tag. Skipping gene."
            rm -f "${modified_tree_file_abs_path_relax}" # 清理空的修改树文件
            continue
        fi

        # 3. 生成子 sbatch 脚本内容
        echo "Generating sbatch script for ${gene_id} (RELAX) at ${sub_script_path_relax}"
        cat << EOF > "${sub_script_path_relax}"
#!/bin/bash
#SBATCH --job-name=relax_${gene_id}
#SBATCH --partition=short # TODO: 根据需要调整分区
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16 # RELAX 通常也是单CPU
#SBATCH --mem=16G          # TODO: 根据需要调整内存
#SBATCH --output=${SLURM_LOGS_DIR_RELAX}/${gene_id}_relax.out
#SBATCH --error=${SLURM_LOGS_DIR_RELAX}/${gene_id}_relax.err

# --- Begin sub-script for ${gene_id} (RELAX) ---
echo "Running HyPhy RELAX for ${gene_id}"
echo "MSA: ${msa_file_abs_path}"
echo "Tree (modified for RELAX): ${modified_tree_file_abs_path_relax}"
echo "Output JSON (RELAX): ${output_json_abs_path_relax}"

# Make sure necessary parent directories for output exist
mkdir -p "$(dirname "${output_json_abs_path_relax}")"

singularity exec ${SINGULARITY_BIND_OPTS} "${HYPHY_IMAGE}" hyphy CPU=16 relax \\
    --alignment "${msa_file_abs_path}" \\
    --tree "${modified_tree_file_abs_path_relax}" \\
    --output "${output_json_abs_path_relax}" \\
    --code Universal < /dev/null

exit_code=$?
echo "HyPhy RELAX for ${gene_id} finished with exit code \${exit_code}."
# --- End sub-script for ${gene_id} (RELAX) ---
EOF

        # 4. 使子脚本可执行
        chmod +x "${sub_script_path_relax}"

        # 5. 提交子脚本
        echo "Submitting sbatch script for ${gene_id} (RELAX): ${sub_script_path_relax}"
        sbatch "${sub_script_path_relax}"
        echo "-----------------------------------------------------"

    else
        echo "Skipping non-file item: ${msa_file_abs_path}"
    fi
done

echo "All HyPhy RELAX sub-jobs have been generated and submitted." 