#!/bin/bash
#SBATCH --job-name=orthofinder_66
#SBATCH --output=orthofinder_66_%j.out
#SBATCH --error=orthofinder_66_%j.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G

ulimit -s unlimited
echo "Running on $(hostname)"
echo "Starting at $(date)"

# Activate OrthoFinder conda environment if applicable
#conda activate orthofinder

# Define paths
input_dir="/home/yuhangjia/data/AlternativeSplicing/orthofinder_66species_input/longest_proteins"
output_dir="/home/yuhangjia/data/AlternativeSplicing/orthofinder_66species_output"

# Run OrthoFinder version 2.5.5
singularity exec /home/yuhangjia/data/AlternativeSplicing/orthofinder_A_M_D/orthofinder_latest.sif orthofinder -f "${input_dir}" -t 32 -a 8 -M msa -o "${output_dir}"

echo "Ending at $(date)"