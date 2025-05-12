#!/bin/bash

#SBATCH --job-name=nkp46_analysis    # Job name
#SBATCH --nodes=1                           # Number of nodes
#SBATCH --output=seurat_%j.out              # Output log file (%j = job ID)
#SBATCH --error=seurat_%j.err               # Error log file

# Set working directory (optional, where logs will be written)
cd /mnt/lustre/RDS-live/moon/ephemeral/NK_project/scripts/nkp46_scr

# Load Docker module in the HPC environment
module load docker

# Define paths
HOST_DATA_DIR="/mnt/lustre/RDS-live/moon/ephemeral/NK_project"
CONTAINER_DATA_DIR="/home"
R_SCRIPT="/scripts/nkp46_scr/nkp46_dimplot_20250331_pm.r"

# Run the Docker container with the R script
docker run --rm \
  -v ${HOST_DATA_DIR}:${CONTAINER_DATA_DIR} \
  satijalab/seurat:5.0.0 \
  Rscript ${CONTAINER_DATA_DIR}/${R_SCRIPT}
