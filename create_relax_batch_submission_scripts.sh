#!/bin/bash
# shellcheck disable=SC2086,SC2046

# ==============================================================================
# SCRIPT: create_relax_batch_submission_scripts.sh
# DESCRIPTION: This script automates the process of running HyPhy RELAX analysis
#              on a large number of orthogroups. It first calculates the length
#              of each sequence alignment, then sorts and partitions the
#              orthogroups into four batches based on length. Finally, it
#              generates tailored SLURM submission scripts for each batch,
#              targeting different partitions (medium, short, epyc, rome) with
#              specific resource allocations and job submission strategies
#              (SLURM array vs. individual jobs).
#
# USAGE: ./create_relax_batch_submission_scripts.sh
#
# REQUIREMENTS:
#   - HyPhy singularity image.
#   - A directory with aligned codon FASTA files.
#   - A directory with corresponding newick tree files (foreground branches
#     marked with '{test}').
#   - SLURM workload manager on the cluster.
#
# AUTHOR: Gemini AI
# DATE: 2024-07-29
# ==============================================================================

set -e # Exit immediately if a command exits with a non-zero status.
set -u # Treat unset variables as an error when substituting.

# --- Configuration Section ---
# PLEASE VERIFY AND ADJUST THESE PATHS BEFORE RUNNING THE SCRIPT.

# Directory containing the aligned codon sequences (e.g., *_codon.clipkit.fasta).
MSA_DIR="aligned_codon_replace_stop_codon_70_coverage"

# Directory containing the corresponding tree files.
# IMPORTANT: Each tree file must have its foreground branches pre-marked with '{test}'.
# The script expects tree files to be named following the pattern: ${orthogroup_id}_some_suffix.tree
# Please adjust the 'TREE_FILE_PATTERN' variable below if your naming scheme is different.
TREE_DIR="subtrees_with_foreground" # Assuming tree files are in the current directory.

# Pattern to find the tree file for a given orthogroup ID.
# The "{OG}" placeholder will be replaced with the actual orthogroup ID (e.g., OG0000001).
# Example: if your trees are named 'OG0000001_marked.tree', use "{OG}_marked.tree".
TREE_FILE_PATTERN="{OG}_marked.treefile"

# Path to the HyPhy singularity image file.
HYPHY_IMAGE="/usr/local/biotools/h/hyphy:2.5.65--he91c24d_0"

# Singularity bind options. Adjust if you need to mount other paths.
SINGULARITY_BIND_OPTS="-B /lustre10:/lustre10"

# Main output directory where all generated scripts and results will be stored.
# A timestamp is used to ensure each run has a unique directory.
MAIN_OUTPUT_DIR="hyphy_relax_submission_$(date +%Y%m%d_%H%M%S)"


# --- Main Script Logic ---

echo "Starting HyPhy RELAX batch preparation script."
echo "-----------------------------------------------------"

# 1. Validate input directories and files
# -----------------------------------------------------
if [ ! -d "${MSA_DIR}" ]; then
    echo "ERROR: MSA directory not found at '${MSA_DIR}'."
    echo "Please set the 'MSA_DIR' variable to the correct path."
    exit 1
fi
if [ ! -d "${TREE_DIR}" ]; then
    echo "ERROR: Tree directory not found at '${TREE_DIR}'."
    echo "Please set the 'TREE_DIR' variable to the correct path."
    exit 1
fi
if [ ! -f "${HYPHY_IMAGE}" ]; then
    echo "ERROR: HyPhy singularity image not found at '${HYPHY_IMAGE}'."
    echo "Please ensure the path is correct."
    exit 1
fi

MSA_FILES_COUNT=$(ls -1q "${MSA_DIR}"/*_codon.clipkit.fasta 2>/dev/null | wc -l)
if [ "${MSA_FILES_COUNT}" -eq 0 ]; then
    echo "ERROR: No MSA files matching '*_codon.clipkit.fasta' were found in '${MSA_DIR}'."
    exit 1
fi
echo "Found ${MSA_FILES_COUNT} MSA files to process."
echo "-----------------------------------------------------"


# 2. Create output directories
# -----------------------------------------------------
echo "Creating output directory structure under '${MAIN_OUTPUT_DIR}'..."
BATCH_LISTS_DIR="${MAIN_OUTPUT_DIR}/batch_lists"
SUBMISSION_SCRIPTS_DIR="${MAIN_OUTPUT_DIR}/submission_scripts"
JSON_OUTPUT_DIR="${MAIN_OUTPUT_DIR}/relax_json_results"
SLURM_LOGS_DIR="${MAIN_OUTPUT_DIR}/slurm_logs"

mkdir -p "${BATCH_LISTS_DIR}"
mkdir -p "${SUBMISSION_SCRIPTS_DIR}"
mkdir -p "${JSON_OUTPUT_DIR}"
mkdir -p "${SLURM_LOGS_DIR}"
echo "Directories created successfully."
echo "-----------------------------------------------------"


# 3. Calculate sequence lengths and create a sorted list
# -----------------------------------------------------
LENGTH_FILE="${BATCH_LISTS_DIR}/orthogroup_lengths.txt"
echo "Calculating sequence lengths for all ${MSA_FILES_COUNT} orthogroups..."
echo "This may take a moment..."

# Ensure the length file is empty before starting
> "${LENGTH_FILE}"

for msa_file in "${MSA_DIR}"/*_codon.clipkit.fasta; do
    if [ -f "${msa_file}" ]; then
        # Extract the length of the first sequence (including gaps)
        length=$(awk '/^>/ {if(p){exit}; p=1;next} {seq=seq $0} END{print length(seq)}' "${msa_file}")
        # Extract the base name (orthogroup ID)
        base_name=$(basename "${msa_file}" _codon.clipkit.fasta)
        echo "${length} ${base_name}" >> "${LENGTH_FILE}"
    fi
done

SORTED_LIST_FILE="${BATCH_LISTS_DIR}/orthogroups_sorted_by_length.txt"
# Sort by length (first column), numerically and in reverse order
sort -k1,1nr "${LENGTH_FILE}" > "${SORTED_LIST_FILE}"

# Clean up intermediate file
rm "${LENGTH_FILE}"

TOTAL_OGS=$(wc -l < "${SORTED_LIST_FILE}")
echo "Finished calculating lengths and sorted ${TOTAL_OGS} orthogroups."
echo "-----------------------------------------------------"


# 4. Split the sorted list into batches
# -----------------------------------------------------
echo "Partitioning ${TOTAL_OGS} orthogroups into 4 batches based on new logic..."
echo "  - The SHORTEST 50% will be assigned to the 'short' partition."
echo "  - The remaining (longest 50%) will be split among 'medium', 'epyc', and 'rome'."

# Calculate batch sizes
# Note: bc rounds down, so there might be a few OGs left for the last batch (Rome).
N_SHORT=$(printf "%.0f" $(echo "scale=0; ${TOTAL_OGS} * 0.50" | bc))
N_MEDIUM=$(printf "%.0f" $(echo "scale=0; ${TOTAL_OGS} * 0.05" | bc)) # Longest 5% of total
N_EPYC=$(printf "%.0f" $(echo "scale=0; ${TOTAL_OGS} * 0.25" | bc)) # Next 25% of total

# Define batch list files
LIST_MEDIUM="${BATCH_LISTS_DIR}/batch_medium.txt"
LIST_SHORT="${BATCH_LISTS_DIR}/batch_short.txt"
LIST_EPYC="${BATCH_LISTS_DIR}/batch_epyc.txt"
LIST_ROME="${BATCH_LISTS_DIR}/batch_rome.txt"

# The sorted list is from LONGEST to SHORTEST.
# awk '{print $2}' is used to extract just the orthogroup ID.

# SHORT batch: The last N_SHORT lines (shortest OGs)
tail -n ${N_SHORT} "${SORTED_LIST_FILE}" | awk '{print $2}' > "${LIST_SHORT}"

# MEDIUM batch: The first N_MEDIUM lines (longest OGs)
head -n ${N_MEDIUM} "${SORTED_LIST_FILE}" | awk '{print $2}' > "${LIST_MEDIUM}"

# EPYC batch: The next N_EPYC lines after the MEDIUM batch
tail -n +$((N_MEDIUM + 1)) "${SORTED_LIST_FILE}" | head -n ${N_EPYC} | awk '{print $2}' > "${LIST_EPYC}"

# ROME batch: The remaining OGs between EPYC and SHORT batches.
# These are the lines from (N_MEDIUM + N_EPYC + 1) down to the start of the SHORT batch.
# We take all OGs and clip off the top (medium, epyc) and the bottom (short).
OFFSET_START=$((N_MEDIUM + N_EPYC))
NUM_LONGER_HALF=$((TOTAL_OGS - N_SHORT))
head -n ${NUM_LONGER_HALF} "${SORTED_LIST_FILE}" | tail -n +$((OFFSET_START + 1)) | awk '{print $2}' > "${LIST_ROME}"


# Recalculate exact counts from the generated files for SLURM arrays and reporting.
N_MEDIUM=$(wc -l < "${LIST_MEDIUM}")
N_SHORT=$(wc -l < "${LIST_SHORT}")
N_EPYC=$(wc -l < "${LIST_EPYC}")
N_ROME=$(wc -l < "${LIST_ROME}")

echo "Batch sizes:"
echo "  - Medium (Longest ~5%):   ${N_MEDIUM} orthogroups"
echo "  - Epyc   (Longer ~25%):   ${N_EPYC} orthogroups"
echo "  - Rome   (Intermediate):  ${N_ROME} orthogroups"
echo "  - Short  (Shortest ~50%): ${N_SHORT} orthogroups"
echo "Batch lists created in '${BATCH_LISTS_DIR}'."
echo "-----------------------------------------------------"


# --- Generate Submission Scripts ---

# Absolute paths for use inside submission scripts
MSA_DIR_ABS=$(realpath ${MSA_DIR})
TREE_DIR_ABS=$(realpath ${TREE_DIR})
JSON_OUTPUT_DIR_ABS=$(realpath ${JSON_OUTPUT_DIR})
SLURM_LOGS_DIR_ABS=$(realpath ${SLURM_LOGS_DIR})
BATCH_LISTS_DIR_ABS=$(realpath ${BATCH_LISTS_DIR})
SUBMISSION_SCRIPTS_DIR_ABS=$(realpath ${SUBMISSION_SCRIPTS_DIR})


# 5. Generate submission script for 'medium' partition (SLURM Array)
# --------------------------------------------------------------------
echo "Generating submission script for MEDIUM partition..."
SCRIPT_MEDIUM="${SUBMISSION_SCRIPTS_DIR}/run_relax_medium_array.sh"

cat << EOF > "${SCRIPT_MEDIUM}"
#!/bin/bash
#SBATCH --job-name=relax_medium
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --array=1-${N_MEDIUM}%50 # Run up to 50 jobs concurrently
#SBATCH --output=${SLURM_LOGS_DIR_ABS}/medium_array_%A_%a.out
#SBATCH --error=${SLURM_LOGS_DIR_ABS}/medium_array_%A_%a.err

# Get the orthogroup ID for this task
OG_ID=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "${BATCH_LISTS_DIR_ABS}/batch_medium.txt")
echo "SLURM_ARRAY_TASK_ID: \${SLURM_ARRAY_TASK_ID}, Orthogroup ID: \${OG_ID}"

# Define file paths
MSA_FILE="${MSA_DIR_ABS}/\${OG_ID}_codon.clipkit.fasta"
TREE_FILE_RAW_PATTERN="${TREE_DIR_ABS}/${TREE_FILE_PATTERN}"
TREE_FILE=\${TREE_FILE_RAW_PATTERN/\{OG\}/\${OG_ID}}
JSON_OUTPUT_FILE="${JSON_OUTPUT_DIR_ABS}/\${OG_ID}.RELAX.json"

# Check for files
if [ ! -f "\${MSA_FILE}" ]; then echo "ERROR: MSA file not found: \${MSA_FILE}"; exit 1; fi
if [ ! -f "\${TREE_FILE}" ]; then echo "ERROR: Tree file not found: \${TREE_FILE}"; exit 1; fi

# Run HyPhy RELAX
echo "Running HyPhy RELAX for \${OG_ID}"
singularity exec ${SINGULARITY_BIND_OPTS} "${HYPHY_IMAGE}" hyphy CPU=\${SLURM_CPUS_PER_TASK} relax \\
    --alignment "\${MSA_FILE}" \\
    --tree "\${TREE_FILE}" \\
    --test test \\
    --output "\${JSON_OUTPUT_FILE}" \\
    --code Universal < /dev/null

echo "Job for \${OG_ID} finished with exit code \$?."
EOF
chmod +x "${SCRIPT_MEDIUM}"
echo " -> Created ${SCRIPT_MEDIUM}"


# 6. Generate submission script for 'short' partition (Individual Jobs)
# ----------------------------------------------------------------------
echo "Generating submission script for SHORT partition..."
SCRIPT_SHORT="${SUBMISSION_SCRIPTS_DIR}/run_relax_short_individual.sh"
SUB_SCRIPTS_DIR_SHORT_ABS="${SUBMISSION_SCRIPTS_DIR_ABS}/short_sub_scripts"
mkdir -p "${SUB_SCRIPTS_DIR_SHORT_ABS}"

cat << EOF > "${SCRIPT_SHORT}"
#!/bin/bash
echo "Submitting ${N_SHORT} individual jobs to the SHORT partition..."
while IFS= read -r OG_ID; do
    if [ -z "\${OG_ID}" ]; then continue; fi
    SUB_SCRIPT="${SUB_SCRIPTS_DIR_SHORT_ABS}/\${OG_ID}.sh"
    
    # Define file paths
    MSA_FILE="${MSA_DIR_ABS}/\${OG_ID}_codon.clipkit.fasta"
    TREE_FILE_RAW_PATTERN="${TREE_DIR_ABS}/${TREE_FILE_PATTERN}"
    TREE_FILE=\${TREE_FILE_RAW_PATTERN/\{OG\}/\${OG_ID}}
    JSON_OUTPUT_FILE="${JSON_OUTPUT_DIR_ABS}/\${OG_ID}.RELAX.json"
    
    # Check for files before creating script
    if [ ! -f "\${MSA_FILE}" ]; then echo "WARNING: MSA file not found, skipping \${OG_ID}: \${MSA_FILE}"; continue; fi
    if [ ! -f "\${TREE_FILE}" ]; then echo "WARNING: Tree file not found, skipping \${OG_ID}: \${TREE_FILE}"; continue; fi

    # Create a specific sbatch script for this OG
    cat << SBATCH_EOF > "\${SUB_SCRIPT}"
#!/bin/bash
#SBATCH --job-name=relax_s_\${OG_ID}
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=64G
#SBATCH --output=${SLURM_LOGS_DIR_ABS}/short_\${OG_ID}_%j.out
#SBATCH --error=${SLURM_LOGS_DIR_ABS}/short_\${OG_ID}_%j.err

echo "Running HyPhy RELAX for \${OG_ID}"
singularity exec ${SINGULARITY_BIND_OPTS} "${HYPHY_IMAGE}" hyphy CPU=\${SLURM_CPUS_PER_TASK} relax \\
    --alignment "\${MSA_FILE}" \\
    --tree "\${TREE_FILE}" \\
    --test test \\
    --output "\${JSON_OUTPUT_FILE}" \\
    --code Universal < /dev/null
echo "Job for \${OG_ID} finished with exit code \$?."
SBATCH_EOF

    # Submit the job
    sbatch "\${SUB_SCRIPT}"
    sleep 0.1 # Small delay to avoid overwhelming the SLURM controller
done < "${BATCH_LISTS_DIR_ABS}/batch_short.txt"
echo "All jobs for the SHORT partition have been submitted."
EOF
chmod +x "${SCRIPT_SHORT}"
echo " -> Created ${SCRIPT_SHORT}"


# 7. Generate submission script for 'epyc' partition (SLURM Array)
# ------------------------------------------------------------------
echo "Generating submission script for EPYC partition..."
SCRIPT_EPYC="${SUBMISSION_SCRIPTS_DIR}/run_relax_epyc_array.sh"

cat << EOF > "${SCRIPT_EPYC}"
#!/bin/bash
#SBATCH --job-name=relax_epyc
#SBATCH --partition=epyc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --array=1-${N_EPYC}%100 # Run up to 100 jobs concurrently
#SBATCH --output=${SLURM_LOGS_DIR_ABS}/epyc_array_%A_%a.out
#SBATCH --error=${SLURM_LOGS_DIR_ABS}/epyc_array_%A_%a.err

# Get the orthogroup ID for this task
OG_ID=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "${BATCH_LISTS_DIR_ABS}/batch_epyc.txt")
echo "SLURM_ARRAY_TASK_ID: \${SLURM_ARRAY_TASK_ID}, Orthogroup ID: \${OG_ID}"

# Define file paths
MSA_FILE="${MSA_DIR_ABS}/\${OG_ID}_codon.clipkit.fasta"
TREE_FILE_RAW_PATTERN="${TREE_DIR_ABS}/${TREE_FILE_PATTERN}"
TREE_FILE=\${TREE_FILE_RAW_PATTERN/\{OG\}/\${OG_ID}}
JSON_OUTPUT_FILE="${JSON_OUTPUT_DIR_ABS}/\${OG_ID}.RELAX.json"

# Check for files
if [ ! -f "\${MSA_FILE}" ]; then echo "ERROR: MSA file not found: \${MSA_FILE}"; exit 1; fi
if [ ! -f "\${TREE_FILE}" ]; then echo "ERROR: Tree file not found: \${TREE_FILE}"; exit 1; fi

# Run HyPhy RELAX
echo "Running HyPhy RELAX for \${OG_ID}"
singularity exec ${SINGULARITY_BIND_OPTS} "${HYPHY_IMAGE}" hyphy CPU=\${SLURM_CPUS_PER_TASK} relax \\
    --alignment "\${MSA_FILE}" \\
    --tree "\${TREE_FILE}" \\
    --test test \\
    --output "\${JSON_OUTPUT_FILE}" \\
    --code Universal < /dev/null

echo "Job for \${OG_ID} finished with exit code \$?."
EOF
chmod +x "${SCRIPT_EPYC}"
echo " -> Created ${SCRIPT_EPYC}"


# 8. Generate submission script for 'rome' partition (SLURM Array)
# ----------------------------------------------------------------
echo "Generating submission script for ROME partition..."
SCRIPT_ROME="${SUBMISSION_SCRIPTS_DIR}/run_relax_rome_array.sh"

cat << EOF > "${SCRIPT_ROME}"
#!/bin/bash
#SBATCH --job-name=relax_rome
#SBATCH --partition=rome
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --array=1-${N_ROME}%100 # Run up to 100 jobs concurrently
#SBATCH --output=${SLURM_LOGS_DIR_ABS}/rome_array_%A_%a.out
#SBATCH --error=${SLURM_LOGS_DIR_ABS}/rome_array_%A_%a.err

# Get the orthogroup ID for this task
OG_ID=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "${BATCH_LISTS_DIR_ABS}/batch_rome.txt")
echo "SLURM_ARRAY_TASK_ID: \${SLURM_ARRAY_TASK_ID}, Orthogroup ID: \${OG_ID}"

# Define file paths
MSA_FILE="${MSA_DIR_ABS}/\${OG_ID}_codon.clipkit.fasta"
TREE_FILE_RAW_PATTERN="${TREE_DIR_ABS}/${TREE_FILE_PATTERN}"
TREE_FILE=\${TREE_FILE_RAW_PATTERN/\{OG\}/\${OG_ID}}
JSON_OUTPUT_FILE="${JSON_OUTPUT_DIR_ABS}/\${OG_ID}.RELAX.json"

# Check for files
if [ ! -f "\${MSA_FILE}" ]; then echo "ERROR: MSA file not found: \${MSA_FILE}"; exit 1; fi
if [ ! -f "\${TREE_FILE}" ]; then echo "ERROR: Tree file not found: \${TREE_FILE}"; exit 1; fi

# Run HyPhy RELAX
echo "Running HyPhy RELAX for \${OG_ID}"
singularity exec ${SINGULARITY_BIND_OPTS} "${HYPHY_IMAGE}" hyphy CPU=\${SLURM_CPUS_PER_TASK} relax \\
    --alignment "\${MSA_FILE}" \\
    --tree "\${TREE_FILE}" \\
    --test test \\
    --output "\${JSON_OUTPUT_FILE}" \\
    --code Universal < /dev/null

echo "Job for \${OG_ID} finished with exit code \$?."
EOF
chmod +x "${SCRIPT_ROME}"
echo " -> Created ${SCRIPT_ROME}"

# --- Final Instructions ---
echo "-----------------------------------------------------"
echo "✅ All submission scripts have been generated successfully!"
echo ""
echo " IMPORTANT: Before proceeding, please check the configuration at the top of this script"
echo "            and inside the generated scripts in '${SUBMISSION_SCRIPTS_DIR}' to ensure all paths and settings are correct."
echo ""
echo "Next Steps:"
echo "1. Navigate to the submission scripts directory:"
echo "   cd ${SUBMISSION_SCRIPTS_DIR}"
echo ""
echo "2. Execute the scripts to submit your jobs to SLURM:"
echo "   ./run_relax_medium_array.sh"
echo "   ./run_relax_short_individual.sh"
echo "   ./run_relax_epyc_array.sh"
echo "   ./run_relax_rome_array.sh"
echo ""
echo "3. Monitor your jobs using 'squeue', 'sinfo', etc."
echo "   JSON results will be saved in: ${JSON_OUTPUT_DIR_ABS}"
echo "   SLURM logs will be saved in:   ${SLURM_LOGS_DIR_ABS}"
echo "-----------------------------------------------------" 