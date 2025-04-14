#!/bin/bash

#SBATCH --job-name=totalNK_wo_51_52         # Job name
#SBATCH --nodes=1                           # Number of nodes
#SBATCH --output=seurat_%j.out              # Output log file (%j = job ID)
#SBATCH --error=seurat_%j.err               # Error log file

# Record the start time in seconds since epoch
START_TIME=$(date +%s)
echo "Job started at: $(date)"

# Set working directory (optional, where logs will be written)
cd /mnt/lustre/RDS-live/moon/ephemeral/NK_project/scripts/total_nk/2_gex_scr/without_Animal_51_52

# Load Docker module in the HPC environment
module load docker

# Define paths
HOST_DATA_DIR="/mnt/lustre/RDS-live/moon/ephemeral/NK_project"
CONTAINER_DATA_DIR="/home"
R_SCRIPT="scripts/total_nk/2_gex_scr/without_Animal_51_52/gene_count.r"

# Run the Docker container with the R script and remove the container after execution
docker run --rm \
  -v ${HOST_DATA_DIR}:${CONTAINER_DATA_DIR} \
  satijalab/seurat:5.0.0 \
  Rscript ${CONTAINER_DATA_DIR}/${R_SCRIPT} || {
    echo "❌ Error: Docker run failed"
    exit 1
  }

# Record the end time
END_TIME=$(date +%s)
echo "Job ended at: $(date)"

# Calculate the elapsed time in seconds
ELAPSED_TIME=$((END_TIME - START_TIME))

# Convert elapsed time to a human-readable format (days, hours, minutes, seconds)
DAYS=$((ELAPSED_TIME / 86400))
HOURS=$(( (ELAPSED_TIME % 86400) / 3600 ))
MINUTES=$(( (ELAPSED_TIME % 3600) / 60 ))
SECONDS=$((ELAPSED_TIME % 60))

echo "Total running time: ${DAYS} days, ${HOURS} hours, ${MINUTES} minutes, ${SECONDS} seconds"

echo "✅ Job completed successfully"