# Headline numbers for poster/talk text, recomputed from the benchmark tables every run.
# Output: analysis/tables/bench_headline.tsv (metric, value), also printed.

source("figures/00_config.R")

sm   <- fread(file.path(TAB_DIR, "bench_summary.tsv"))
tc   <- fread(file.path(TAB_DIR, "bench_target_calls.tsv"), colClasses = list(character = c("cpu_filter", "gpu_filter")))
disc <- fread(file.path(TAB_DIR, "bench_discordant.tsv"))
pooled <- function(t) sm[, (sum(get(paste0(t, "_n"))) - sum(get(paste0(t, "_disc")))) / sum(get(paste0(t, "_n")))]
rng <- function(x, d = 1) sprintf(paste0("%.", d, "f (%.", d, "f–%.", d, "f)"), median(x), min(x), max(x))
sh <- tc[!is.na(cpu_AF) & !is.na(gpu_AF)]
bam <- sm[!is.na(bam_records)]

h <- data.table(metric = c(
  "cases", "read pairs, millions: median (range)",
  "end-to-end speedup: median (range)", "Mutect2 speedup: median (range)", "alignment+dedup speedup: median (range)",
  "CPU minutes per case: median (range)", "GPU minutes per case: median (range)",
  "CPU hours, all cases", "GPU hours, all cases",
  "flagstat identical (cases)", "BAM records compared (cases)", "BAM records compared", "BAM records not identical",
  "BAM records in one arm only",
  "genome-wide raw Jaccard, pooled", "genome-wide PASS Jaccard, pooled", "panel raw Jaccard, pooled", "panel PASS Jaccard, pooled",
  "panel calls (union), raw", "panel discordant calls", "genome-wide discordant calls", "discordant calls with depth <= 3",
  "discordant PASS calls with depth >= 10", "panel VAF r (shared calls)", "panel VAF identical", "panel max |dVAF|"),
  value = c(
  nrow(sm), rng(sm$read_pairs / 1e6), rng(sm$total_speedup), rng(sm$m2_speedup), rng(sm$bam_speedup),
  rng(sm$cpu_total_s / 60, 0), rng(sm$gpu_total_s / 60), sprintf("%.1f", sum(sm$cpu_total_s) / 3600),
  sprintf("%.1f", sum(sm$gpu_total_s) / 3600),
  sprintf("%d / %d", sum(sm$flagstat_identical), nrow(sm)), nrow(bam), comma(sum(bam$bam_records)),
  comma(sum(bam$bam_nonidentical)), comma(sum(bam$bam_only_one_arm)),
  sprintf("%.5f", pooled("genome_raw")), sprintf("%.5f", pooled("genome_pass")),
  sprintf("%.5f", pooled("target_raw")), sprintf("%.5f", pooled("target_pass")),
  comma(sum(sm$target_raw_n)), sum(sm$target_raw_disc), nrow(disc),
  sprintf("%d (%.1f%%)", disc[DP <= 3, .N], 100 * disc[, mean(DP <= 3)]), disc[filter == "PASS" & DP >= 10, .N],
  sprintf("%.6f", cor(sh$cpu_AF, sh$gpu_AF)), sprintf("%.2f%%", 100 * mean(sh$cpu_AF == sh$gpu_AF)),
  sprintf("%.4f", max(abs(sh$cpu_AF - sh$gpu_AF)))))
fwrite(h, file.path(TAB_DIR, "bench_headline.tsv"), sep = "\t")
print(h, row.names = FALSE)
