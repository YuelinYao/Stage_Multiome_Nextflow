# Planned option: running with containers

The workflow currently runs with Conda (`-profile conda` or `-profile slurm,conda`). A container option would make it portable to other clusters. It is not implemented yet; this page records the plan and the cluster constraints found so far.

## Cluster constraint

Running Singularity as an ordinary user requires either a setuid installation or Linux user namespaces. The `singularity/3.8.6` module is a Conda installation without setuid, so it needs user namespaces:

| Node | `cat /proc/sys/user/max_user_namespaces` | `singularity exec docker://alpine:3.19 …` |
| --- | --- | --- |
| sledge01 | 15000 | Works |
| sledge02, sledge03, sledge04, sledge14, lambda ... | 0 | Fails: `Failed to create user namespace: user namespace disabled` |


Until the HPC administrators enable user namespaces or a setuid Singularity/Apptainer on more nodes, container jobs would have to run on `sledge01` only, so Conda remains the practical choice on this cluster.

## Implementation plan

1.  **Build one image** with the same software as the `stage_multiome` environment: R 4.3, Seurat 5.3.0, Signac 1.14.0, DoubletFinder 2.0.4, cellXY, SCP, glmGamPoi, Rsamtools, EnsDb.Mmusculus.v79, BSgenome.Mmusculus.UCSC.mm10, JASPAR2020, TFBSTools, chromVAR, motifmatchr and MACS3. Use a `containers/Dockerfile` (or Singularity `.def`) based on `micromamba`: install the Conda packages from `environment.yml`, then the R-installed packages listed in its first line, all with pinned versions. Build on a machine with Docker (building needs root or user namespaces), push to a registry, and pull once on the cluster with `singularity pull stage_multiome.sif docker://…` into a shared folder. Expect an image of roughly 5–8 GB, mostly from the mm10 genome package.
    
2.  **Add a `singularity` profile to `nextflow.config`**, alongside `conda`:
    
    ```groovy
    params.container = null            // path to the .sif, or a docker:// address
    singularity {                      // use: -profile slurm,singularity
      singularity.enabled    = true
      singularity.autoMounts = true
      singularity.runOptions = '-B /suffolk -B /nfs4'
      process.container      = params.container
    }
    ```
    
    The bind mounts are required: the R scripts read ARC, mgatk and fragment files at their original absolute paths and run scripts from the pipeline folder. Do not combine the `singularity` and `conda` profiles.
    
3.  **Relax host checks in `main.nf`.** The workflow requires `macs3_path` to be an executable on the host and checks `conda_env`. In container mode, MACS3 is inside the image, so `macs3_path` should default to `macs3` on the container’s `PATH`, and both host checks should be skipped.
    
4.  **Cluster settings in your copy of `site.config`:** send jobs to nodes that allow containers (currently `process.clusterOptions = '--nodelist=sledge01'`), make Singularity available in each job (for example by adding `module add singularity/3.8.6` to the existing `beforeScript`), keep `SINGULARITY_TMPDIR` on local disk (`/tmp`) and the image cache (`NXF_SINGULARITY_CACHEDIR`) on shared storage.
    
5.  **Validate** by running `tests/*.R` inside the image (`singularity exec stage_multiome.sif Rscript tests/…`), then running the README vignette’s three test samples with both `-profile slurm,conda` and `-profile slurm,singularity` and comparing the outputs (for example `wnn_mgatk.rds` metadata and the heteroplasmy and depth tables). Differences point to package versions that do not match.
    
6.  **Document** the profile, image location and pull command in the main [README](../README.md).
    

## Considerations

-   Steps 2–4 and 6 are small (about 20 lines of configuration and code plus documentation). Most of the effort and risk is in step 1: reproducing the exact package versions, especially the packages installed from R rather than Conda.
-   While only `sledge01` allows containers, all container jobs queue for one node and large runs will be slow. The benefit is a portable workflow that produces identical results on other clusters.
