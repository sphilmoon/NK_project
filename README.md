# NK_project
Gene expression analysis of NK cell subsets based on the presence or absence of the NKp46 marker.

## Using docker container for Seurat R package

1. **Running a published image and remove the container once the job is finished.**
    ```bash
    # To run seurat for GEX analysis
    docker run -it --rm -v /mnt/lustre/RDS-live/moon/ephemeral/NK_project:/home satijalab/seurat:5.0.0 bash