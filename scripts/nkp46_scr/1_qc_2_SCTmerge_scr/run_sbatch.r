#!/bin/bash

#SBATCH --job-name=nkp46                # Job name
#SBATCH --nodes=1                           # Number of nodes
#SBATCH --output=nkp46_qc_SCT.out              # Output log file (%j = job ID)
#SBATCH --error=nkp46_qc_SCT.err               # Error log file

# set working directory
cd /mnt/lustre/RDS-live/moon/ephemeral/NK_project/outputs

# load docker in the HPC
module load docker

# define paths
HOST_DATA_DIR="/mnt/lustre/RDS-live/moon/ephemeral/NK_project"
CONTAINER_DATA_DIR="/home"
R_SCRIPT="/scripts/nkp46_scr/1_qc_2_SCTmerge_scr/1_qc_2_SCTmerge.r"

# run the docker container with the R script with chat gpt.
docker run --rm \
    -v ${HOST_DATA_DIR}:${CONTAINER_DATA_DIR} \
    satijalab/seurat:5.0.0 \
    Rscript ${CONTAINER_DATA_DIR}/${R_SCRIPT}