#!/bin/bash
echo "Submitting 147 individual jobs to the SHORT partition..."
while IFS= read -r OG_ID; do
    if [ -z "${OG_ID}" ]; then continue; fi
    SUB_SCRIPT="/lustre10/home/yuhangjia/data/AlternativeSplicing/evo_rate_test_RNA_splicing_term/hyphy_relax/submission_scripts/${OG_ID}.sh"

    # Define file paths
    MSA_FILE="/lustre10/home/yuhangjia/data/AlternativeSplicing/evo_rate_test_RNA_splicing_term/GO_0008380_msa_codon_clipkit_for_paml/${OG_ID}_codon.clipkit.fasta"
    TREE_FILE_RAW_PATTERN="/lustre10/home/yuhangjia/data/AlternativeSplicing/all_OGs_for_relax_test/subtrees_with_foreground/{OG}_marked.treefile"
    TREE_FILE=${TREE_FILE_RAW_PATTERN/\{OG\}/${OG_ID}}
    JSON_OUTPUT_FILE="/lustre10/home/yuhangjia/data/AlternativeSplicing/evo_rate_test_RNA_splicing_term/hyphy_relax/relax_json_results/${OG_ID}.RELAX.json"

    # Check for files before creating script
    if [ ! -f "${MSA_FILE}" ]; then echo "WARNING: MSA file not found, skipping ${OG_ID}: ${MSA_FILE}"; continue; fi
    if [ ! -f "${TREE_FILE}" ]; then echo "WARNING: Tree file not found, skipping ${OG_ID}: ${TREE_FILE}"; continue; fi

    # Create a specific sbatch script for this OG
    cat << SBATCH_EOF > "${SUB_SCRIPT}"
#!/bin/bash
#SBATCH --job-name=relax_s_${OG_ID}
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=128G
#SBATCH --output=/lustre10/home/yuhangjia/data/AlternativeSplicing/evo_rate_test_RNA_splicing_term/hyphy_relax/slurm_logs/short_${OG_ID}_%j.out
#SBATCH --error=/lustre10/home/yuhangjia/data/AlternativeSplicing/evo_rate_test_RNA_splicing_term/hyphy_relax/slurm_logs/short_${OG_ID}_%j.err

echo "Running HyPhy RELAX for ${OG_ID}"
singularity exec -B /lustre10:/lustre10 "/usr/local/biotools/h/hyphy:2.5.65--he91c24d_0" hyphy CPU=${SLURM_CPUS_PER_TASK} relax \
    --alignment "${MSA_FILE}" \
    --tree "${TREE_FILE}" \
    --test test \
    --output "${JSON_OUTPUT_FILE}" \
    --code Universal < /dev/null
echo "Job for ${OG_ID} finished with exit code $?."
SBATCH_EOF

    # Submit the job
    sbatch "${SUB_SCRIPT}"
    sleep 0.1 # Small delay to avoid overwhelming the SLURM controller
done < "/lustre10/home/yuhangjia/data/AlternativeSplicing/evo_rate_test_RNA_splicing_term/prefix.txt"
echo "All jobs for the SHORT partition have been submitted."