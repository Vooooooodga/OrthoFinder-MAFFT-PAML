#!/bin/bash

#SBATCH --job-name=hyphy_meme_PARENT # 主脚本的作业名
#SBATCH --partition=your_partition      # 主脚本运行的分区 (可能与子作业不同)
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1             # 主脚本本身不需要很多CPU
#SBATCH --mem=2G                    # 主脚本本身不需要很多内存
#SBATCH --time=02:00:00               # 主脚本提交作业的时间，根据基因数量调整
#SBATCH --output=hyphy_meme_parent_%j.out
#SBATCH --error=hyphy_meme_parent_%j.err

# 加载 singularity 模块 (如果需要)
# module load singularity # 取消注释并根据您的集群环境修改

# --- 全局定义 ---
# TODO: 用户需要根据实际情况修改以下路径
MSA_DIR="/home/yuhangjia/data/AlternativeSplicing/9_overlap_DAS_evo_rate_analysis/DAS_aligned_codon_clipkit_for_paml"  
TREE_DIR="/home/yuhangjia/data/AlternativeSplicing/9_overlap_DAS_evo_rate_analysis/gene_trees_from_M0_remarked_v2" 
MAIN_OUTPUT_DIR_MEME="/home/yuhangjia/data/AlternativeSplicing/9_overlap_DAS_evo_rate_analysis/MEME" # 主输出目录 for MEME
HYPHY_IMAGE="/usr/local/biotools/h/hyphy:2.5.65--he91c24d_0" # HyPhy singularity 镜像路径
SINGULARITY_BIND_OPTS="-B /lustre10:/lustre10" # 根据需要修改 Singularity 绑定路径

# 新的目录，用于存放生成的子 sbatch 脚本
SUB_SBATCH_SCRIPT_DIR_MEME="${MAIN_OUTPUT_DIR_MEME}/sub_sbatch_scripts_meme"
# 新的目录，用于存放为MEME去除了标记的树文件
UNMARKED_TREE_DIR_MEME="${MAIN_OUTPUT_DIR_MEME}/hyphy_unmarked_trees_meme"
# 子作业的输出将直接写入 MAIN_OUTPUT_DIR_MEME/MEME_json (json本身) 和 MAIN_OUTPUT_DIR_MEME/MEME_slurm_logs
MEME_JSON_OUTPUT_DIR="${MAIN_OUTPUT_DIR_MEME}/MEME_json"
SLURM_LOGS_DIR_MEME="${MAIN_OUTPUT_DIR_MEME}/MEME_slurm_logs"

# --- 创建所需目录 ---
echo "Creating directories for MEME analysis..."
mkdir -p "${SUB_SBATCH_SCRIPT_DIR_MEME}"
mkdir -p "${UNMARKED_TREE_DIR_MEME}"
mkdir -p "${MEME_JSON_OUTPUT_DIR}"
mkdir -p "${SLURM_LOGS_DIR_MEME}"

echo "Starting to generate and submit HyPhy MEME sub-jobs..."

# --- 主循环：遍历MSA文件，生成并提交子脚本 ---
# DEBUG: List all files matching the pattern in MSA_DIR
echo "DEBUG: Files found in ${MSA_DIR} matching *_codon.clipkit.fasta:"
ls -1 "${MSA_DIR}"/*_codon.clipkit.fasta 2>/dev/null | tee /dev/stderr | wc -l
echo "DEBUG: ---- End of file list ----"

for msa_file_abs_path in "${MSA_DIR}"/*_codon.clipkit.fasta; do
    echo "DEBUG: Current msa_file_abs_path variable is: [${msa_file_abs_path}]"
    if [ -f "${msa_file_abs_path}" ]; then
        base_name=$(basename "${msa_file_abs_path}")
        gene_id=${base_name%_codon.clipkit.fasta}

        # 原始树文件路径 (MEME将直接使用此树，无需修改)
        original_tree_file_abs_path="${TREE_DIR}/${gene_id}_from_M0_marked.treefile"
        # 为MEME准备的去除了标记的树文件路径
        unmarked_tree_file_abs_path_meme="${UNMARKED_TREE_DIR_MEME}/${gene_id}_from_M0_unmarked_for_meme.treefile"
        # 输出的JSON文件路径 (子脚本将使用此路径)
        output_json_abs_path_meme="${MEME_JSON_OUTPUT_DIR}/${gene_id}.MEME.json"
        
        # 子 sbatch 脚本的路径
        sub_script_path_meme="${SUB_SBATCH_SCRIPT_DIR_MEME}/sub_meme_${gene_id}.sh"

        # 1. 检查原始 tree 文件是否存在
        if [ ! -f "${original_tree_file_abs_path}" ]; then
            echo "WARNING: Original tree file not found for ${gene_id}: ${original_tree_file_abs_path}. Skipping gene."
            continue
        fi

        # 2. 为 MEME 准备输入树文件.
        #    MEME分析的输入树不应包含前景分支标记 (如 #1).
        #    因此, 我们将从原始树文件中移除这些标记, 并生成一个新的、无标记的树文件.
        echo "Preparing unmarked tree file for ${gene_id} (MEME)..."
        echo "Original (marked) tree: ${original_tree_file_abs_path}"
        # 移除类似 #1, #2 等标记
        sed 's/#[0-9]\+//g' "${original_tree_file_abs_path}" > "${unmarked_tree_file_abs_path_meme}"

        if [ ! -s "${unmarked_tree_file_abs_path_meme}" ]; then
            echo "ERROR: Failed to create unmarked tree file for ${gene_id} (MEME) or the file is empty. Original: ${original_tree_file_abs_path}. Skipping gene."
            rm -f "${unmarked_tree_file_abs_path_meme}" # 清理空的树文件
            continue
        fi
        echo "Unmarked tree for MEME: ${unmarked_tree_file_abs_path_meme}"

        # 3. 生成子 sbatch 脚本内容
        echo "Generating sbatch script for ${gene_id} (MEME) at ${sub_script_path_meme}"
        cat << EOF > "${sub_script_path_meme}"
#!/bin/bash
#SBATCH --job-name=meme_${gene_id}
#SBATCH --partition=short # TODO: 根据需要调整分区
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32 # MEME 通常也是单CPU或可多CPU，根据hyphy版本和MEME实现调整
#SBATCH --mem=32G          # TODO: 根据需要调整内存
#SBATCH --output=${SLURM_LOGS_DIR_MEME}/${gene_id}_meme.out
#SBATCH --error=${SLURM_LOGS_DIR_MEME}/${gene_id}_meme.err

# --- Begin sub-script for ${gene_id} (MEME) ---
echo "Running HyPhy MEME for ${gene_id}"
echo "MSA: ${msa_file_abs_path}"
echo "Tree (original): ${original_tree_file_abs_path}"
echo "Tree (unmarked for MEME): ${unmarked_tree_file_abs_path_meme}"
echo "Output JSON (MEME): ${output_json_abs_path_meme}"

# Make sure necessary parent directories for output exist
mkdir -p "$(dirname "${output_json_abs_path_meme}")"

singularity exec ${SINGULARITY_BIND_OPTS} "${HYPHY_IMAGE}" hyphy CPU=32 meme \\
    --alignment "${msa_file_abs_path}" \\
    --tree "${unmarked_tree_file_abs_path_meme}" \\
    --branches All \\
    --output "${output_json_abs_path_meme}" \\
    --code Universal \\
    --pvalue 0.1 < /dev/null # 使用默认的 p-value，也可以按需修改

exit_code=$?
echo "HyPhy MEME for ${gene_id} finished with exit code \${exit_code}."
# --- End sub-script for ${gene_id} (MEME) ---
EOF

        # 4. 使子脚本可执行
        chmod +x "${sub_script_path_meme}"

        # 5. 提交子脚本
        echo "Submitting sbatch script for ${gene_id} (MEME): ${sub_script_path_meme}"
        sbatch "${sub_script_path_meme}"
        echo "-----------------------------------------------------"

    else
        echo "Skipping non-file item: ${msa_file_abs_path}"
    fi
done

echo "All HyPhy MEME sub-jobs have been generated and submitted." 