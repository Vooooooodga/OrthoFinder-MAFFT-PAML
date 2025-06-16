#!/bin/bash

#SBATCH --job-name=hyphy_busted-ph_PARENT # 主脚本的作业名
#SBATCH --partition=your_partition      # 主脚本运行的分区 (可能与子作业不同)
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1             # 主脚本本身不需要很多CPU
#SBATCH --mem=2G                    # 主脚本本身不需要很多内存
#SBATCH --time=02:00:00               # 主脚本提交作业的时间，根据基因数量调整
#SBATCH --output=hyphy_busted-ph_parent_%j.out
#SBATCH --error=hyphy_busted-ph_parent_%j.err

# 加载 singularity 模块 (如果需要)
# module load singularity # 取消注释并根据您的集群环境修改

# --- 全局定义 ---
MSA_DIR="/home/yuhangjia/data/AlternativeSplicing/evo_rate_test_RNA_splicing_term/GO_0008380_msa_codon_clipkit_for_paml"
TREE_DIR="/home/yuhangjia/data/AlternativeSplicing/evo_rate_test_RNA_splicing_term/gene_trees_from_M0_remarked_v2"
MAIN_OUTPUT_DIR="/home/yuhangjia/data/AlternativeSplicing/evo_rate_test_RNA_splicing_term/BUSTED-PH" # 主输出目录
HYPHY_IMAGE="/usr/local/biotools/h/hyphy:2.5.65--he91c24d_0"
SINGULARITY_BIND_OPTS="-B /lustre10:/lustre10" # 根据需要修改

# 新的目录，用于存放修改后的树文件 (由父脚本统一处理)
MODIFIED_TREE_DIR="${MAIN_OUTPUT_DIR}/hyphy_modified_trees"
# 新的目录，用于存放生成的子 sbatch 脚本
SUB_SBATCH_SCRIPT_DIR="${MAIN_OUTPUT_DIR}/sub_sbatch_scripts"
# 子作业的输出将直接写入 MAIN_OUTPUT_DIR/BUSTED-PH_json (json本身) 和 MAIN_OUTPUT_DIR/BUSTED-PH_slurm_logs
BUSTED_PH_JSON_OUTPUT_DIR="${MAIN_OUTPUT_DIR}/BUSTED-PH_json"
SLURM_LOGS_DIR="${MAIN_OUTPUT_DIR}/BUSTED-PH_slurm_logs"

# --- 创建所需目录 ---
echo "Creating directories..."
mkdir -p "${MODIFIED_TREE_DIR}"
mkdir -p "${SUB_SBATCH_SCRIPT_DIR}"
mkdir -p "${BUSTED_PH_JSON_OUTPUT_DIR}"
mkdir -p "${SLURM_LOGS_DIR}"

echo "Starting to generate and submit HyPhy BUSTED-PH sub-jobs..."

# --- 主循环：遍历MSA文件，预处理树，生成并提交子脚本 ---
# DEBUG: List all files matching the pattern in MSA_DIR
echo "DEBUG: Files found in ${MSA_DIR} matching *_codon.clipkit.fasta:"
ls -1 "${MSA_DIR}"/*_codon.clipkit.fasta 2>/dev/null | tee /dev/stderr | wc -l
echo "DEBUG: ---- End of file list ----"

for msa_file_abs_path in "${MSA_DIR}"/*_codon.clipkit.fasta; do
    echo "DEBUG: Current msa_file_abs_path variable is: [${msa_file_abs_path}]"
    if [ -f "${msa_file_abs_path}" ]; then
        base_name=$(basename "${msa_file_abs_path}")
        gene_id=$(echo "$base_name" | sed 's/_codon\.clipkit\.fasta$//')

        original_tree_file_abs_path="${TREE_DIR}/${gene_id}_from_M0_marked.treefile"
        # 输出的JSON文件路径 (子脚本将使用此路径)
        output_json_abs_path="${BUSTED_PH_JSON_OUTPUT_DIR}/${gene_id}.BUSTED-PH.json"
        # 修改后的树文件路径 (子脚本将使用此路径)
        modified_tree_file_abs_path="${MODIFIED_TREE_DIR}/${gene_id}_from_M0_marked.treefile"
        
        # 子 sbatch 脚本的路径
        sub_script_path="${SUB_SBATCH_SCRIPT_DIR}/sub_busted-ph_${gene_id}.sh"

        # 1. 检查原始 tree 文件是否存在
        if [ ! -f "${original_tree_file_abs_path}" ]; then
            echo "WARNING: Original tree file not found for ${gene_id}: ${original_tree_file_abs_path}. Skipping gene."
            continue
        fi

        # 2. 预处理树文件 (由父脚本完成)
        echo "Preprocessing tree file for ${gene_id}: Replacing #1 with {foreground}"
        sed 's/#1/{foreground}/g' "${original_tree_file_abs_path}" > "${modified_tree_file_abs_path}"
        if [ ! -s "${modified_tree_file_abs_path}" ]; then
            echo "ERROR: Failed to preprocess tree file for ${gene_id} or modified tree is empty. Original: ${original_tree_file_abs_path}. Skipping gene."
            rm -f "${modified_tree_file_abs_path}"
            continue
        fi

        # 3. 生成子 sbatch 脚本内容
        echo "Generating sbatch script for ${gene_id} at ${sub_script_path}"
        cat << EOF > "${sub_script_path}"
#!/bin/bash
#SBATCH --job-name=busted-ph_${gene_id}
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=${SLURM_LOGS_DIR}/${gene_id}_busted-ph.out
#SBATCH --error=${SLURM_LOGS_DIR}/${gene_id}_busted-ph.err

# --- Begin sub-script for ${gene_id} ---
echo "Running HyPhy BUSTED-PH for ${gene_id}"
echo "MSA: ${msa_file_abs_path}"
echo "Tree: ${modified_tree_file_abs_path}"
echo "Output JSON: ${output_json_abs_path}"

# Make sure necessary parent directories for output exist (though parent script should create BUSTED_PH_JSON_OUTPUT_DIR)
mkdir -p "$(dirname "${output_json_abs_path}")"

singularity exec ${SINGULARITY_BIND_OPTS} "${HYPHY_IMAGE}" hyphy CPU=1 busted BUSTED-PH.bf \\
    --alignment "${msa_file_abs_path}" \\
    --tree "${modified_tree_file_abs_path}" \\
    --output "${output_json_abs_path}" \\
    --code Universal \\
    --branches foreground < /dev/null

exit_code=$?
echo "HyPhy BUSTED-PH for ${gene_id} finished with exit code ${exit_code}."
# --- End sub-script for ${gene_id} ---
EOF

        # 4. 使子脚本可执行
        chmod +x "${sub_script_path}"

        # 5. 提交子脚本
        echo "Submitting sbatch script for ${gene_id}: ${sub_script_path}"
        sbatch "${sub_script_path}"
        echo "-----------------------------------------------------"

    else
        echo "Skipping non-file item: ${msa_file_abs_path}"
    fi
done

echo "All HyPhy BUSTED-PH sub-jobs have been generated and submitted." 