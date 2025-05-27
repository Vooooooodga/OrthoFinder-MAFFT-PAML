#!/bin/bash
#SBATCH --job-name=iqtree_66
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=iqtree_infer_spec_tree_%j.out
#SBATCH --error=iqtree_infer_spec_tree_%j.err

singularity exec /usr/local/biotools/i/iqtree:2.3.6--hdbdd923_0 iqtree \
-s ./concatenated_aa.fasta \
-p ./partitions_aa.txt \
-g ./astral.tree \
-B 1000 -T AUTO --prefix 66species_astral_topo2 -o Drosophila_melanogaster

echo "Ending at $(date)"