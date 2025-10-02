#!/bin/bash -l

# --------------------------------------------------------------------------------
# SGE Job Directives
# --------------------------------------------------------------------------------
# Project name
#$ -P cancergrp

# Hard time limit: hh:mm:ss.
#$ -l h_rt=12:00:00

# Job name
#$ -N SequenceSimilarityCalculation

# Merge stdout and stderr into a single file (.o#jobID).
#$ -j y

# Email notification (job begins (b), ends (e), is aborted (a), suspended (s), or never (n))
# 'be' for both start and end.
#$ -m e

#$ -M phitro@bu.edu

# Request parallel environment for shared memory applications (like multiprocessing.Pool).
#$ -pe omp 32


# Request memory per core. Total memory = N * mem_per_core.
#$ -l mem_per_core=8G

# --------------------------------------------------------------------------------
# Set up Conda environment
# --------------------------------------------------------------------------------

echo "Loading Conda environment..."
module load miniconda

# Activate specific Conda environment
conda activate parasail_env

# Verify environment is active (optional, good for debugging)
echo "Conda environment active: $CONDA_DEFAULT_ENV (prefix: $CONDA_PREFIX)"
echo "Current Python environment: $(which python)"
echo "Python version: $(python --version)"

# --------------------------------------------------------------------------------
# Run Python script
# --------------------------------------------------------------------------------

echo "Job started on host: $(hostname)"
echo "Requested number of CPU cores (NSLOTS): $NSLOTS"

# Run Python script with arguments.
python TF-sequence-similarity-calculation.py \
    --num_processes $NSLOTS \
    --input /projectnb/cancergrp/Philipp/data/pairs_htf_htf_with_sequences.csv \
    --output /projectnb/cancergrp/Philipp/data/pairs_htf_htf_with_sequences_with_sequence-similarity.csv \

# --------------------------------------------------------------------------------
# Post-job cleanup
# --------------------------------------------------------------------------------

echo "Python script finished."

# Deactivate the Conda environment. 
if [ -n "$CONDA_PREFIX" ]; then
    conda deactivate
    echo "Conda environment deactivated."
fi
