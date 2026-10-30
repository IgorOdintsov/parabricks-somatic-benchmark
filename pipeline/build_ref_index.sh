#!/bin/bash
#SBATCH --job-name=ref_index
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=160G
#SBATCH --time=08:00:00
#SBATCH --output=logs/ref_index_%j.log

# hg19 indexes for both arms. Usage: sbatch build_ref_index.sh <DIR_WITH_hg19.fa_fai_dict>
#   BWA 0.7.18 index  -> used by Parabricks fq2bam
#   BWA-MEM2 2.2.1    -> CPU arm
set -euo pipefail
source "${SLURM_SUBMIT_DIR}/config.sh"   # CONTAINERS, REF, SCRATCH, RESULTS, BIND (see config.example.sh)
REFDIR=$(dirname ${REF})
SRC=${1:?usage: sbatch build_ref_index.sh <DIR_WITH_hg19.fa_fai_dict>}
mkdir -p ${REFDIR}
cp -L ${SRC}/hg19.fa ${REFDIR}/hg19.fa
cp ${SRC}/hg19.fa.fai ${SRC}/hg19.dict ${REFDIR}/
cd ${REFDIR}
( T0=$(date +%s); apptainer exec ${BIND} ${CONTAINERS}/bwa_0.7.18.sif bwa index -a bwtsw hg19.fa > bwa_index.log 2>&1
  echo "bwa index seconds: $(( $(date +%s) - T0 ))" ) &
( T0=$(date +%s); apptainer exec ${BIND} ${CONTAINERS}/bwamem2_2.2.1.sif bwa-mem2 index hg19.fa > bwamem2_index.log 2>&1
  echo "bwa-mem2 index seconds: $(( $(date +%s) - T0 ))" ) &
wait
ls -la ${REFDIR}
echo "Done."
