#!/bin/bash

#SBATCH --job-name=hyphy_aBSREL_PARENT # 主脚本的作业名
#SBATCH --partition=your_partition      # 主脚本运行的分区 (可能与子作业不同)
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1             # 主脚本本身不需要很多CPU
#SBATCH --mem=2G                    # 主脚本本身不需要很多内存
#SBATCH --time=02:00:00               # 主脚本提交作业的时间，根据基因数量调整
#SBATCH --output=hyphy_aBSREL_parent_%j.out
#SBATCH --error=hyphy_aBSREL_parent_%j.err

# 加载 singularity 模块 (如果需要)
# module load singularity # 取消注释并根据您的集群环境修改

# --- 全局定义 ---
MSA_DIR="/home/yuhangjia/data/AlternativeSplicing/9_overlap_DAS_evo_rate_analysis/DAS_aligned_codon_clipkit_for_paml"
TREE_DIR="/home/yuhangjia/data/AlternativeSplicing/9_overlap_DAS_evo_rate_analysis/gene_trees_from_M0_remarked_v2"
MAIN_OUTPUT_DIR="/home/yuhangjia/data/AlternativeSplicing/9_overlap_DAS_evo_rate_analysis/aBSREL" # 主输出目录
HYPHY_IMAGE="/usr/local/biotools/h/hyphy:2.5.65--he91c24d_0"
SINGULARITY_BIND_OPTS="-B /lustre10:/lustre10" # 根据需要修改

# 新的目录，用于存放生成的子 sbatch 脚本
SUB_SBATCH_SCRIPT_DIR="${MAIN_OUTPUT_DIR}/sub_sbatch_scripts"
# 子作业的输出将直接写入 MAIN_OUTPUT_DIR/aBSREL_json (json本身) 和 MAIN_OUTPUT_DIR/aBSREL_slurm_logs
ABSREL_JSON_OUTPUT_DIR="${MAIN_OUTPUT_DIR}/aBSREL_json"
SLURM_LOGS_DIR="${MAIN_OUTPUT_DIR}/aBSREL_slurm_logs"

# --- 创建所需目录 ---
echo "Creating directories..."
mkdir -p "${SUB_SBATCH_SCRIPT_DIR}"
mkdir -p "${ABSREL_JSON_OUTPUT_DIR}"
mkdir -p "${SLURM_LOGS_DIR}"

echo "Starting to generate and submit HyPhy aBSREL sub-jobs..."

# --- 主循环：遍历MSA文件，生成并提交子脚本 ---
for msa_file_abs_path in "${MSA_DIR}"/*_codon.clipkit.fasta; do
    if [ -f "${msa_file_abs_path}" ]; then
        base_name=$(basename "${msa_file_abs_path}")
        gene_id=$(echo "$base_name" | sed 's/_codon\.clipkit\.fasta$//')

        original_tree_file_abs_path="${TREE_DIR}/${gene_id}_from_M0_marked.treefile"
        # 输出的JSON文件路径 (子脚本将使用此路径)
        output_json_abs_path="${ABSREL_JSON_OUTPUT_DIR}/${gene_id}.aBSREL.json"
        
        # 子 sbatch 脚本的路径
        sub_script_path="${SUB_SBATCH_SCRIPT_DIR}/sub_absrel_${gene_id}.sh"

        # 1. 检查原始 tree 文件是否存在
        if [ ! -f "${original_tree_file_abs_path}" ]; then
            echo "WARNING: Original tree file not found for ${gene_id}: ${original_tree_file_abs_path}. Skipping gene."
            continue
        fi

        # 2. 生成子 sbatch 脚本内容
        echo "Generating sbatch script for ${gene_id} at ${sub_script_path}"
        cat << EOF > "${sub_script_path}"
#!/bin/bash
#SBATCH --job-name=absrel_${gene_id}
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --output=${SLURM_LOGS_DIR}/${gene_id}_absrel.out
#SBATCH --error=${SLURM_LOGS_DIR}/${gene_id}_absrel.err

# --- Begin sub-script for ${gene_id} ---
echo "Running HyPhy aBSREL for ${gene_id}"
echo "MSA: ${msa_file_abs_path}"
echo "Tree: ${original_tree_file_abs_path}"
echo "Output JSON: ${output_json_abs_path}"

# Make sure necessary parent directories for output exist
mkdir -p "$(dirname "${output_json_abs_path}")"

singularity exec ${SINGULARITY_BIND_OPTS} "${HYPHY_IMAGE}" hyphy CPU=32 aBSREL \\
    --alignment "${msa_file_abs_path}" \\
    --tree "${original_tree_file_abs_path}" \\
    --branches All \\
    --output "${output_json_abs_path}" \\
    --code Universal < /dev/null

exit_code=$?
echo "HyPhy aBSREL for ${gene_id} finished with exit code ${exit_code}."
# --- End sub-script for ${gene_id} ---
EOF

        # 3. 使子脚本可执行
        chmod +x "${sub_script_path}"

        # 4. 提交子脚本
        echo "Submitting sbatch script for ${gene_id}: ${sub_script_path}"
        sbatch "${sub_script_path}"
        echo "-----------------------------------------------------"

    else
        echo "Skipping non-file item: ${msa_file_abs_path}"
    fi
done

echo "All HyPhy aBSREL sub-jobs have been generated and submitted." 