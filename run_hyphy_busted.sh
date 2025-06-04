#!/bin/bash

#SBATCH --job-name=hyphy_busted
#SBATCH --partition=your_partition # 请替换为您的集群分区名称
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1    # 运行一个主脚本任务
#SBATCH --cpus-per-task=64     # 为主脚本请求所有CPU，脚本内部管理并发
#SBATCH --mem=128G             # 根据您的需求调整内存
#SBATCH --time=72:00:00        # 根据您的需求调整运行时间
#SBATCH --output=hyphy_busted_%j.out
#SBATCH --error=hyphy_busted_%j.err

# 加载 singularity 模块 (如果需要)
# module load singularity # 取消注释并根据您的集群环境修改

# 定义输入输出目录
MSA_DIR="/home/yuhangjia/data/AlternativeSplicing/evo_rate_test_RNA_splicing_term/GO_0008380_msa_codon_clipkit_for_paml"
TREE_DIR="/home/yuhangjia/data/AlternativeSplicing/evo_rate_test_RNA_splicing_term/gene_trees_from_M0_remarked_v2"
OUTPUT_DIR="/home/yuhangjia/data/AlternativeSplicing/evo_rate_test_RNA_splicing_term/BUSTED"
HYPHY_IMAGE="/usr/local/biotools/h/hyphy:2.5.65--he91c24d_0" # Singularity镜像路径

# Singularity 绑定路径 (根据您的集群和数据位置进行修改)
# 如果您的MSA、TREE或OUTPUT目录不在Singularity默认绑定的路径下(如 /home/$USER, /tmp, $PWD),
# 您需要在这里明确指定绑定。
# 示例: SINGULARITY_BIND_OPTS="-B /path/to/data_on_host:/data_in_container -B /another/path:/another/path"
# 但为了明确和以防万一，您可以绑定一个共同的父目录，例如:
# SINGULARITY_BIND_OPTS="-B /home/yuhangjia/data:/home/yuhangjia/data"
SINGULARITY_BIND_OPTS="-B /lustre10:/lustre10" # Singularity 绑定路径, 根据需要修改, e.g. "-B /path1:/path1 -B /path2:/path2"

# 创建输出目录 (如果不存在)
mkdir -p "${OUTPUT_DIR}"

# 并发控制参数
# 总共请求了 64 CPUs，每个 hyphy 进程使用 4 CPUs
# 所以 MAX_JOBS = 64 / 4 = 16
MAX_JOBS=16
job_count=0

echo "Starting HyPhy BUSTED analysis..."
echo "MSA Directory: ${MSA_DIR}"
echo "Tree Directory: ${TREE_DIR}"
echo "Output Directory: ${OUTPUT_DIR}"
echo "Max concurrent jobs: ${MAX_JOBS}"

# DEBUG: List all files matching the pattern in MSA_DIR
echo "DEBUG: Files found in ${MSA_DIR} matching *_codon.clipkit.fasta:"
ls -1 "${MSA_DIR}"/*_codon.clipkit.fasta 2>/dev/null | tee /dev/stderr | wc -l
echo "DEBUG: ---- End of file list ----"

# 遍历 MSA 文件 (假设文件以 .fasta 结尾)
# 用户提供的文件名格式是 OG0002009_codon.clipkit.fasta
for msa_file in "${MSA_DIR}"/*_codon.clipkit.fasta; do
    echo "DEBUG: Current msa_file variable is: [$msa_file]"
    if [ -f "$msa_file" ]; then
        base_name=$(basename "$msa_file")
        # 从 "OG0002009_codon.clipkit.fasta" 提取 "OG0002009"
        gene_id=$(echo "$base_name" | sed 's/_codon\.clipkit\.fasta$//')

        # 构建对应的 tree 文件名
        # 用户提供的文件名格式是 OG0001125_from_M0_marked.treefile
        tree_file="${TREE_DIR}/${gene_id}_from_M0_marked.treefile"
        output_json="${OUTPUT_DIR}/${gene_id}.BUSTED.json"

        # 检查 tree 文件是否存在
        if [ ! -f "$tree_file" ]; then
            echo "WARNING: Tree file not found for ${gene_id}: ${tree_file}. Skipping."
            continue
        fi

        echo "Processing ${gene_id}: MSA='${base_name}', Tree='$(basename "${tree_file}")'"

        # 运行 HyPhy BUSTED 命令
        # 使用 --CPU 4 来指定每个 hyphy 进程使用的核心数
        # 命令在后台运行 (&)
        singularity exec ${SINGULARITY_BIND_OPTS} "${HYPHY_IMAGE}" hyphy busted \\
            --alignment "${msa_file}" \\
            --tree "${tree_file}" \\
            --output "${output_json}" \\
            --CPU 4 &

        job_count=$((job_count + 1))
        echo "Launched job for ${gene_id}. Current job count: ${job_count}"

        # 如果达到最大并发作业数，则等待任一作业完成
        if [ "${job_count}" -ge "${MAX_JOBS}" ]; then
            echo "Reached max concurrent jobs (${MAX_JOBS}). Waiting for a slot..."
            wait -n # 等待任何一个后台任务结束
            job_count=$((job_count - 1)) # 减少计数器，因为一个作业已完成
            echo "Slot freed. Current job count: ${job_count}"
        fi
    else
        echo "Skipping non-file item: ${msa_file}"
    fi
done

# 等待所有剩余的后台作业完成
echo "All jobs launched. Waiting for remaining jobs to complete..."
wait
echo "All HyPhy BUSTED analyses completed."
echo "Output files are located in ${OUTPUT_DIR}" 