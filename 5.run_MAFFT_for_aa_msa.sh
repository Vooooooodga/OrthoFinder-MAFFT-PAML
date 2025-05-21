#!/bin/bash

#$ -S /bin/bash

#$ -cwd

#$ -pe def_slot 64

#$ -l s_vmem=20G

#$ -l mem_req=20G

  

ulimit -s unlimited

echo "Running on $(hostname)"

echo "Starting at $(date)"

  

# 输入和输出目录定义

input="./translated_proteins_deg"

output="./aligned_translated_proteins_deg"

  

# 确保输出目录存在

mkdir -p $output

  

# 检查输入文件

echo "Checking input files in $input:"

ls $input

  

# 为目录中的每个文件运行MAFFT

for file in $input/*.fa; do

echo "Processing $file..."

outfile="$output/$(basename "${file%.fa}")_aligned.fa"

echo "Output will be saved to $outfile"

singularity exec /usr/local/biotools/m/mafft:7.525--h031d066_0 mafft --auto --thread 64 "$file" > "$outfile"

done

  

echo "Ending at $(date)"