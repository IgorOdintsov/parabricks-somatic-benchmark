#!/bin/bash
#SBATCH --job-name=bench_gpu
#SBATCH --partition=gpu
#SBATCH --gres=gpu:nvidia_l40s:1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=01:30:00
#SBATCH --output=logs/bench_gpu_%A_%a.log

# AMP 2026 benchmark, GPU arm. Same hardware as the CPU arm: L40S GPU nodes (AMD EPYC 9354P),
# 32 CPUs, 128 GB, plus 1 L40S. Config matched to bench_cpu.sh: lanes concatenated, RG with PU,
# Parabricks 4.4.0-1 fq2bam --bwa-options="-K 10000000" (no -Y), mutectcaller genome-wide with downsampling off,
# FilterMutectCalls 4.3.0.0, bcftools norm. Every step timed.
# Usage: sbatch --array=1-N bench_gpu.sh <CASE_LIST> <RUN_DIR>   (case = line N of list, FASTQs in RUN_DIR/<case>)
set -euo pipefail
source "${SLURM_SUBMIT_DIR}/config.sh"   # CONTAINERS, REF, SCRATCH, RESULTS, BIND (see config.example.sh)
S=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$1")
FQDIR=$2/${S}
WORK=${SCRATCH}/${S}/gpu
KEEP=${RESULTS}/${S}
G="apptainer exec ${BIND} ${CONTAINERS}/gatk_4.3.0.0.sif"
B="apptainer exec ${BIND} ${CONTAINERS}/bcftools_1.21.sif"
PB="apptainer exec --nv ${BIND} ${CONTAINERS}/parabricks_4.4.0-1.sif pbrun"
RG="@RG\tID:${S}\tSM:${S}\tPL:ILLUMINA\tLB:${S}\tPU:${S}"  # fq2bam requires PU
T=gpu_nods
TIMING=${KEEP}/${S}_gpu_timing.tsv
rm -rf ${WORK}; mkdir -p ${WORK} ${KEEP}
[ $(ls ${FQDIR}/*_L00[12]_R[12]_001.fastq.gz 2>/dev/null | wc -l) -eq 4 ] || { echo "Expected 4 FASTQs in ${FQDIR}"; exit 1; }

# node context: who else was on the node (the node is shared, not exclusive)
node() { echo -e "${S}\tgpu\t$1\t$(hostname)\t$(scontrol show node $(hostname) | grep -oE 'CPUAlloc=[0-9]+')\t$(squeue -w $(hostname) -h | wc -l) jobs"; }
step() {  # name cmd...: wall-clock seconds per step
  local N=$1; shift; local T0=$(date +%s)
  "$@"
  echo -e "${S}\tgpu\t${N}\t$(( $(date +%s) - T0 ))" >> ${TIMING}
}
echo -e "sample\tarm\tstep\tseconds" > ${TIMING}
lscpu | grep "Model name"; nvidia-smi --query-gpu=name,driver_version --format=csv,noheader
node start | tee ${KEEP}/${S}_gpu_node.tsv

concat() {
  cat ${FQDIR}/*_L001_R1_001.fastq.gz ${FQDIR}/*_L002_R1_001.fastq.gz > ${WORK}/${S}_R1.fastq.gz
  cat ${FQDIR}/*_L001_R2_001.fastq.gz ${FQDIR}/*_L002_R2_001.fastq.gz > ${WORK}/${S}_R2.fastq.gz
}
fq2bam() {  # align + sort + mark duplicates + index
  ${PB} fq2bam --ref ${REF} --in-fq ${WORK}/${S}_R1.fastq.gz ${WORK}/${S}_R2.fastq.gz "${RG}" \
    --bwa-options="-K 10000000" --out-bam ${WORK}/${S}_gpu.bam \
    --out-duplicate-metrics ${KEEP}/${S}_gpu_markdup_metrics.txt \
    --num-gpus 1 --tmp-dir ${WORK} > ${WORK}/fq2bam.log 2>&1
}
mutect2() {
  ${PB} mutectcaller --ref ${REF} --in-tumor-bam ${WORK}/${S}_gpu.bam --tumor-name ${S} \
    --mutectcaller-options="-max-reads-per-alignment-start 0" \
    --out-vcf ${WORK}/${S}_${T}.vcf.gz --num-gpus 1 > ${WORK}/${T}_mutect.log 2>&1
}
filter() {
  ${B} bcftools index -f -t ${WORK}/${S}_${T}.vcf.gz
  ${G} gatk --java-options "-Xmx8g" FilterMutectCalls -R ${REF} -V ${WORK}/${S}_${T}.vcf.gz \
    --stats ${WORK}/${S}_${T}.vcf.gz.stats -O ${WORK}/${S}_${T}.filtered.vcf.gz > ${WORK}/filter_${T}.log 2>&1
  ${B} bcftools norm -m -any -f ${REF} -Oz -o ${KEEP}/${S}_${T}.filtered.norm.vcf.gz ${WORK}/${S}_${T}.filtered.vcf.gz
}

step concat  concat
step fq2bam  fq2bam
step mutect2 mutect2
step filter  filter
node end | tee -a ${KEEP}/${S}_gpu_node.tsv

${B} bcftools query -f "${T}\t%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/TLOD\t%INFO/ECNT\t[%DP]\t[%AD{1}]\t[%AF]\n" \
  ${KEEP}/${S}_${T}.filtered.norm.vcf.gz | gzip > ${KEEP}/${S}_${T}.calls.tsv.gz
${G} samtools flagstat -@ 8 ${WORK}/${S}_gpu.bam > ${KEEP}/${S}_gpu_flagstat.txt
rm -f ${WORK}/${S}_R1.fastq.gz ${WORK}/${S}_R2.fastq.gz
cp ${WORK}/*.log ${KEEP}/ 2>/dev/null || true
cat ${TIMING}
echo "Done."
