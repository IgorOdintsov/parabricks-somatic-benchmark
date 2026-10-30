#!/bin/bash
#SBATCH --job-name=bench_cpu
#SBATCH --partition=gpu   # same node type as the GPU arm; no GPU requested
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH --output=logs/bench_cpu_%A_%a.log

# AMP 2026 benchmark, CPU arm. Same hardware as the GPU arm: L40S GPU nodes (AMD EPYC 9354P),
# 32 CPUs, 128 GB; no GPU allocated. Config matched to bench_gpu.sh: lanes concatenated, RG with PU,
# bwa-mem2 2.2.1 -K 10000000 (no -Y) -> samtools sort -> GATK MarkDuplicates, GATK 4.3.0.0 Mutect2
# genome-wide with downsampling off, FilterMutectCalls 4.3.0.0, bcftools norm. Every step timed.
# Usage: sbatch --array=1-N bench_cpu.sh <CASE_LIST> <RUN_DIR>   (case = line N of list, FASTQs in RUN_DIR/<case>)
set -euo pipefail
source "${SLURM_SUBMIT_DIR}/config.sh"   # CONTAINERS, REF, SCRATCH, RESULTS, BIND (see config.example.sh)
S=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$1")
FQDIR=$2/${S}
WORK=${SCRATCH}/${S}/cpu
KEEP=${RESULTS}/${S}
G="apptainer exec ${BIND} ${CONTAINERS}/gatk_4.3.0.0.sif"
B="apptainer exec ${BIND} ${CONTAINERS}/bcftools_1.21.sif"
M2="apptainer exec ${BIND} ${CONTAINERS}/bwamem2_2.2.1.sif"
RG="@RG\tID:${S}\tSM:${S}\tPL:ILLUMINA\tLB:${S}\tPU:${S}"  # fq2bam requires PU
T=cpu_nods
TIMING=${KEEP}/${S}_cpu_timing.tsv
rm -rf ${WORK}; mkdir -p ${WORK} ${KEEP}
[ $(ls ${FQDIR}/*_L00[12]_R[12]_001.fastq.gz 2>/dev/null | wc -l) -eq 4 ] || { echo "Expected 4 FASTQs in ${FQDIR}"; exit 1; }

# node context: who else was on the node (the node is shared, not exclusive)
node() { echo -e "${S}\tcpu\t$1\t$(hostname)\t$(scontrol show node $(hostname) | grep -oE 'CPUAlloc=[0-9]+')\t$(squeue -w $(hostname) -h | wc -l) jobs"; }
step() {  # name cmd...: wall-clock seconds per step
  local N=$1; shift; local T0=$(date +%s)
  "$@"
  echo -e "${S}\tcpu\t${N}\t$(( $(date +%s) - T0 ))" >> ${TIMING}
}
echo -e "sample\tarm\tstep\tseconds" > ${TIMING}
lscpu | grep "Model name"; node start | tee ${KEEP}/${S}_cpu_node.tsv

concat() {
  cat ${FQDIR}/*_L001_R1_001.fastq.gz ${FQDIR}/*_L002_R1_001.fastq.gz > ${WORK}/${S}_R1.fastq.gz
  cat ${FQDIR}/*_L001_R2_001.fastq.gz ${FQDIR}/*_L002_R2_001.fastq.gz > ${WORK}/${S}_R2.fastq.gz
}
align_sort() {
  ${M2} bwa-mem2 mem -t 32 -K 10000000 -R "${RG}" ${REF} ${WORK}/${S}_R1.fastq.gz ${WORK}/${S}_R2.fastq.gz \
    2> ${WORK}/bwamem2.log | ${G} samtools sort -@ 8 -m 2G -T ${WORK}/sorttmp -o ${WORK}/${S}_cpu_sorted.bam -
}
markdup() {
  ${G} gatk --java-options "-Xmx32g" MarkDuplicates -I ${WORK}/${S}_cpu_sorted.bam -O ${WORK}/${S}_cpu.bam \
    -M ${KEEP}/${S}_cpu_markdup_metrics.txt --TMP_DIR ${WORK} > ${WORK}/markdup.log 2>&1
}
mutect2() {
  ${G} gatk --java-options "-Xmx16g" Mutect2 -R ${REF} -I ${WORK}/${S}_cpu.bam -tumor ${S} \
    --max-reads-per-alignment-start 0 -O ${WORK}/${S}_${T}.vcf.gz > ${WORK}/${T}_mutect2.log 2>&1
}
filter() {
  ${B} bcftools index -f -t ${WORK}/${S}_${T}.vcf.gz
  ${G} gatk --java-options "-Xmx8g" FilterMutectCalls -R ${REF} -V ${WORK}/${S}_${T}.vcf.gz \
    --stats ${WORK}/${S}_${T}.vcf.gz.stats -O ${WORK}/${S}_${T}.filtered.vcf.gz > ${WORK}/filter_${T}.log 2>&1
  ${B} bcftools norm -m -any -f ${REF} -Oz -o ${KEEP}/${S}_${T}.filtered.norm.vcf.gz ${WORK}/${S}_${T}.filtered.vcf.gz
}

step concat     concat
step align_sort align_sort
step markdup    markdup
step index      ${G} samtools index -@ 8 ${WORK}/${S}_cpu.bam
step mutect2    mutect2
step filter     filter
node end | tee -a ${KEEP}/${S}_cpu_node.tsv

${B} bcftools query -f "${T}\t%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/TLOD\t%INFO/ECNT\t[%DP]\t[%AD{1}]\t[%AF]\n" \
  ${KEEP}/${S}_${T}.filtered.norm.vcf.gz | gzip > ${KEEP}/${S}_${T}.calls.tsv.gz
${G} samtools flagstat -@ 8 ${WORK}/${S}_cpu.bam > ${KEEP}/${S}_cpu_flagstat.txt
rm -f ${WORK}/${S}_cpu_sorted.bam ${WORK}/${S}_R1.fastq.gz ${WORK}/${S}_R2.fastq.gz
cp ${WORK}/*.log ${KEEP}/ 2>/dev/null || true
cat ${TIMING}
echo "Done."
