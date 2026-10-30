# Cost per case from the benchmark runtimes, and the BAM-identity summary panel.
# Prices: AWS us-east-1 Linux on-demand (same as 04_cost.R; verify at aws.amazon.com/ec2/pricing before presenting).
# m6i.8xlarge (32 vCPU, 128 GiB) for CPU; g6e.8xlarge (32 vCPU, 1x L40S) for GPU. Runtime = FASTQ -> unfiltered VCF
# (bench_summary cpu_total_s / gpu_total_s); post-filtering is identical on both arms and excluded.

source("figures/00_config.R")

sm <- fread(file.path(TAB_DIR, "bench_summary.tsv"))
N  <- nrow(sm)
PRICE <- c(CPU = 1.536, GPU = 4.5286)
INST  <- c(CPU = "CPU\nm6i.8xlarge", GPU = "GPU (L40S)\ng6e.8xlarge")

# Cost -------------------------------------------------------------------------------------------
cost <- rbind(sm[, .(sample, arm = "CPU", usd = cpu_total_s / 3600 * PRICE["CPU"])],
              sm[, .(sample, arm = "GPU", usd = gpu_total_s / 3600 * PRICE["GPU"])])
cs <- cost[, .(median_usd = median(usd), min_usd = min(usd), max_usd = max(usd)), by = arm]
cs[, `:=`(usd_hr = PRICE[arm], lab = factor(INST[arm], levels = INST))]
saving <- cs[arm == "CPU", median_usd] / cs[arm == "GPU", median_usd]
fwrite(cs[, !"lab"], file.path(TAB_DIR, "bench_cost.tsv"), sep = "\t")
print(cs)
cost[, lab := factor(INST[arm], levels = INST)]

p_cost <- ggplot(cs, aes(lab, median_usd, fill = arm)) +
  geom_col(width = 0.55) +
  geom_point(data = cost, aes(lab, usd), inherit.aes = FALSE, position = position_jitter(width = 0.08, seed = 3),
             size = 1.6, shape = 21, fill = SURFACE, colour = INK_2, stroke = 0.4) +
  geom_text(aes(label = sprintf("$%.2f\n($%.2f/hr)", median_usd, usd_hr)), position = position_nudge(x = 0.33),
            hjust = 0, size = 3.2, colour = INK, lineheight = 1) +
  scale_fill_manual(values = arm_cols, guide = "none") +
  scale_x_discrete(expand = expansion(add = c(0.5, 0.9))) +
  scale_y_continuous(labels = label_dollar(), expand = expansion(mult = c(0, 0.08)), name = "Compute cost per case (median)") +
  labs(x = NULL, title = sprintf("%.1f× lower compute cost per case", saving),
       subtitle = sprintf("GPU instance costs %.1f× more per hour but finishes %.0f× faster",
                          PRICE["GPU"] / PRICE["CPU"], median(sm$total_speedup)),
       caption = sprintf("AWS us-east-1 on-demand list prices × measured wall-clock time (FASTQ → VCF); n = %d", N)) +
  theme_amp() + theme(panel.grid.major.x = element_blank())
save_fig(p_cost, "fig_bench_cost", 5.5, 4.5)

# BAM identity tiles -----------------------------------------------------------------------------
bam <- sm[!is.na(bam_records)]
tot <- sum(bam$bam_records); nonid <- sum(bam$bam_nonidentical); one <- sum(bam$bam_only_one_arm)
fmt_big <- function(x) if (x >= 1e9) sprintf("%.2f billion", x / 1e9) else sprintf("%.0f million", x / 1e6)
tiles <- data.table(
  x = 1:3,
  big = c(fmt_big(tot), comma(one), comma(nonid)),
  small = c(sprintf("aligned reads compared\nread by read (%d cases)", nrow(bam)),
            "reads present in only\none of the two BAMs",
            sprintf("records differ (CIGAR only):\n1 in %s; positions, MAPQ and\nduplicate flags all identical",
                    if (nonid > 0) comma(round(tot / nonid, -6)) else "none")))
p_bam <- ggplot(tiles) +
  geom_text(aes(x, 1, label = big), size = 9, fontface = "bold", colour = INK, vjust = 0) +
  geom_text(aes(x, 0.82, label = small), size = 3.4, colour = INK_2, vjust = 1, lineheight = 0.95) +
  scale_x_continuous(limits = c(0.4, 3.6)) + scale_y_continuous(limits = c(0.3, 1.35)) +
  labs(title = "CPU and GPU produce the same BAM",
       subtitle = sprintf("BWA-MEM2 + samtools sort + GATK MarkDuplicates vs Parabricks fq2bam; flagstat identical in %d / %d cases",
                          sum(sm$flagstat_identical), N)) +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold", colour = INK, size = rel(1.15), margin = margin(b = 4)),
        plot.subtitle = element_text(colour = INK_2, margin = margin(b = 4)),
        plot.background = element_rect(fill = SURFACE, colour = NA), plot.margin = margin(10, 12, 6, 10))
save_fig(p_bam, "fig_bench_bam", 10, 2.8)
