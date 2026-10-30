#!/bin/bash
#SBATCH --job-name=bam_compare
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/bam_compare_%j.log

# Read-level comparison of the benchmark CPU and GPU BAMs for one case (same logic as fq2bam_test.sh).
# key = QNAME + R1/R2 + supplementary bit (+ chrom:pos for supplementary, since a read can have several);
# value = flag, chrom, pos, MAPQ, CIGAR. Secondary records excluded.
# Usage: sbatch bam_compare.sh <SAMPLE>
set -euo pipefail
source "${SLURM_SUBMIT_DIR}/config.sh"   # CONTAINERS, REF, SCRATCH, RESULTS, BIND (see config.example.sh)
S=$1
W=${SCRATCH}/${S}
KEEP=${RESULTS}/${S}
G="apptainer exec ${BIND} ${CONTAINERS}/gatk_4.3.0.0.sif"
TMP=${W}/bamcmp; mkdir -p ${TMP}

tab() {  # bam -> sorted tsv
  ${G} samtools view -@ 4 -F 0x100 $1 | \
    awk -F'\t' 'BEGIN{OFS="\t"} {r=and($2,64)?1:2; s=and($2,2048)?"S/"$3":"$4:"P"; print $1"/"r"/"s, $2, $3, $4, $5, $6}' | \
    LC_ALL=C sort -k1,1 -S 12G --parallel=6 -T ${TMP} > $2
}
tab ${W}/cpu/${S}_cpu.bam ${TMP}/cpu.tsv &
tab ${W}/gpu/${S}_gpu.bam ${TMP}/gpu.tsv &
wait
LC_ALL=C join -t $'\t' -a1 -a2 -e NA -o 0,1.2,1.3,1.4,1.5,1.6,2.2,2.3,2.4,2.5,2.6 ${TMP}/cpu.tsv ${TMP}/gpu.tsv | \
  awk -F'\t' -v S=${S} '
    { n++
      if ($2=="NA") { gonly++; next } if ($7=="NA") { conly++; next }
      both++
      if ($2==$7 && $3==$8 && $4==$9 && $5==$10 && $6==$11) same++
      else {
        if ($3!=$8 || $4!=$9) dpos++
        if ($6!=$11) dcig++
        if ($5!=$10) dmq++
        if (and($2,1024)!=and($7,1024)) ddup++
        if (xor(and($2,0xFBFF),and($7,0xFBFF))) dflag++
      }
    }
    END { printf "sample=%s records_union=%d both=%d cpu_only=%d gpu_only=%d identical=%d (%.6f) diff_position=%d diff_cigar=%d diff_mapq=%d diff_dupflag=%d diff_other_flags=%d\n",
          S, n, both, conly, gonly, same, same/both, dpos, dcig, dmq, ddup, dflag }' | tee ${KEEP}/${S}_bam_comparison.txt
rm -rf ${TMP}
