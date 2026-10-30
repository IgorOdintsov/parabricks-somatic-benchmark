# Cluster paths for the SLURM scripts in pipeline/. Copy to config.sh (git-ignored), edit, and submit jobs
# from the repository root: the scripts source "${SLURM_SUBMIT_DIR}/config.sh".
CONTAINERS=/path/to/containers   # Apptainer images: parabricks_4.4.0-1.sif, gatk_4.3.0.0.sif, bwamem2_2.2.1.sif,
                                 # bwa_0.7.18.sif, bcftools_1.21.sif
REF=/path/to/ref/hg19.fa         # with .fai, .dict, BWA (fq2bam) and BWA-MEM2 indexes (pipeline/build_ref_index.sh)
SCRATCH=/path/to/scratch         # per-case working directories (BAMs, intermediate VCFs)
RESULTS=/path/to/results         # small per-case outputs kept after the run (timing, VCFs, call tables, metrics)
BIND="--bind /path/to/data,${SCRATCH},${RESULTS},$(dirname ${REF})"   # every path the containers must see, FASTQs included
