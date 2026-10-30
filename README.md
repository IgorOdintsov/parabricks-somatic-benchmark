# GPU vs CPU somatic variant calling: a matched Parabricks benchmark

A head-to-head benchmark of NVIDIA Parabricks against a standard CPU pipeline (BWA-MEM2 + GATK Mutect2) on clinical
tumor-only cancer panel sequencing, from FASTQ to filtered VCF.

> **Status: results embargoed.** This work will be presented at the **AMP 2026 Annual Meeting** (Seattle,
> November 10–14, 2026) as an oral presentation and a poster: *"GPU-Accelerated Alignment and Somatic Variant
> Calling in Cancer Panel Sequencing: A Head-to-Head Benchmark Against a Standard CPU Pipeline."* Results, figures
> and data are not included here until then; please refer to the poster and talk. This repository holds the
> analysis code only.

## Question

Parabricks runs GPU implementations of BWA-MEM2 and GATK Mutect2. Clinical adoption needs more than "similar"
output. **With every setting matched, does the GPU pipeline produce the same alignments and the same variant calls
as the CPU pipeline, and how much faster is it on identical hardware?**

## Design

- **Samples:** clinical tumor-only targeted panel cases from a single sequencing run (2 lanes, paired-end).
- **Hardware:** both arms run on the same node type (L40S GPU nodes, AMD EPYC 9354P) with 32 vCPUs (16 cores) and
  128 GB per job. The GPU arm adds one NVIDIA L40S; the CPU arm requests no GPU.
- **Identical input:** lanes concatenated once per case, one read group (with PU) used by both arms.
- **Timing:** every step is wall-clock timed. FASTQ concatenation and post-filtering are identical on both arms and
  excluded from the totals.

### Matched pipelines

| Step | CPU arm (`pipeline/bench_cpu.sh`) | GPU arm (`pipeline/bench_gpu.sh`) |
|---|---|---|
| Alignment | BWA-MEM2 2.2.1 `mem -K 10000000` (no `-Y`) | Parabricks 4.4.0-1 `fq2bam --bwa-options="-K 10000000"` |
| Sort | samtools sort | (inside fq2bam) |
| Duplicate marking | GATK 4.3.0.0 MarkDuplicates | (inside fq2bam) |
| Calling | GATK 4.3.0.0 Mutect2, genome-wide, `--max-reads-per-alignment-start 0` | Parabricks `mutectcaller`, genome-wide, `-max-reads-per-alignment-start 0` |
| Filtering | GATK 4.3.0.0 FilterMutectCalls | GATK 4.3.0.0 FilterMutectCalls |
| Normalization | `bcftools norm -m -any` | `bcftools norm -m -any` |

GATK 4.3.0.0 is the version bundled in Parabricks 4.4.0, so the caller and filter versions are the same on both arms.
Downsampling is turned off on both arms so that read selection cannot differ between the two implementations.
The batch size `-K` is passed explicitly, as the Parabricks documentation advises for reproducing BWA-MEM2 output.

Each of these settings matters: aligner flags (`-Y`), batch size, duplicate marking tool, downsampling and the
GATK version of the filter all change the call set. A CPU vs GPU comparison is only valid when the whole pipeline is
matched, not just the caller.

## What is compared

1. **Runtime:** per step and end to end, per case; speedup vs read depth.
2. **Alignments:** read-by-read comparison of the two BAMs (`pipeline/bam_compare.sh`). Each primary or
   supplementary record is keyed by read name, mate and position, and compared on flag, position, MAPQ and CIGAR.
3. **Variant calls:** Jaccard index of the call sets (CHROM, POS, REF, ALT after normalization), genome-wide and on
   the panel target region (IDT V5 exons ±10 bp), for all calls and for PASS calls.
4. **Allele fractions:** VAF agreement for calls shared by both arms on the panel.
5. **Cost:** cloud list price per hour (AWS on-demand, m6i.8xlarge vs g6e.8xlarge) × measured wall-clock time.

## Repository layout

```
config.example.sh          site paths (copy to config.sh; git-ignored)
pipeline/
  build_ref_index.sh       hg19 BWA (for fq2bam) and BWA-MEM2 indexes
  bench_cpu.sh             CPU arm, SLURM array job, one case per task
  bench_gpu.sh             GPU arm, SLURM array job, one case per task
  bam_compare.sh           read-level CPU vs GPU BAM comparison for one case
analysis/
  bench_summary.py         per-case timing, BAM and VCF concordance -> analysis/tables/bench_*.tsv
figures/
  00_config.R              palette, theme, paths
  10_bench_runtime.R       runtime figures
  11_bench_concordance.R   BAM / call / VAF concordance figures
  12_bench_headline.R      headline numbers table
  14_bench_cost_bam.R      cost per case
  run_bench_figures.R      runs all of the above
```

## Running it

Requirements: SLURM, Apptainer, and the container images listed in `config.example.sh`. A Parabricks license is not
required for version 4.x. Locally: Python 3 (standard library only) and R with data.table, ggplot2, patchwork, scales.

```bash
cp config.example.sh config.sh                 # edit paths
mkdir -p logs
sbatch pipeline/build_ref_index.sh /path/to/dir_with_hg19.fa   # once

# cases.txt: one case ID per line; FASTQs in <RUN_DIR>/<case>/*_L00{1,2}_R{1,2}_001.fastq.gz
N=$(wc -l < cases.txt)
sbatch --array=1-$N pipeline/bench_gpu.sh cases.txt /path/to/RUN_DIR
sbatch --array=1-$N pipeline/bench_cpu.sh cases.txt /path/to/RUN_DIR
for S in $(cat cases.txt); do sbatch pipeline/bam_compare.sh $S; done   # after both arms finish
```

Then copy each `${RESULTS}/<case>/` directory (small files only: `*.tsv`, `*.tsv.gz`, `*.txt`) to `results/<case>/`,
put the panel target BED at `analysis/tables/panel_targets.bed`, and run:

```bash
python -I analysis/bench_summary.py .
Rscript figures/run_bench_figures.R            # figures -> analysis/figures/
```

## Contact

Igor Odintsov, Clinical Cancer Genomics Laboratory, University of California San Francisco.

Code and documentation were developed with the assistance of Claude Code (Anthropic); analysis design and results were
reviewed by the authors.
