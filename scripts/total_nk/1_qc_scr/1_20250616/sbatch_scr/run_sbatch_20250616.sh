#!/bin/bash

#SBATCH --job-name=merged_dge_complete                    # Job name
#SBATCH --nodes=1                           # Number of nodes
#SBATCH --output=merged_dge_complete.out              # Output log file (%j = job ID)
#SBATCH --error=merged_dge_complete.err               # Error log file

# set working directory
cd /mnt/lustre/RDS-live/moon/ephemeral/NK_project/outputs

# load docker in the HPC
module load docker

# define paths
HOST_DATA_DIR="/mnt/lustre/RDS-live/moon/ephemeral/NK_project"
CONTAINER_DATA_DIR="/home"
R_SCRIPT="/scripts/total_nk/1_qc_scr/1_20250616/3_merged_dge_complete_20250626.r"

# run the docker container with the R script with chat gpt.
docker run --rm \
    -v ${HOST_DATA_DIR}:${CONTAINER_DATA_DIR} \
    satijalab/seurat:5.0.0 \
    Rscript ${CONTAINER_DATA_DIR}/${R_SCRIPT}