# Benchmark runtime figures (v2). Both arms on the same L40S nodes (AMD EPYC 9354P,
# 32 vCPUs = 16 cores, 128 GB); GPU arm adds 1x L40S. Matched configs; Mutect2 genome-wide, downsampling off.
# Input: analysis/tables/bench_{summary,timing_long}.tsv from analysis/bench_summary.py
# (rebuild with `Rscript figures/run_bench_figures.R`).
# FASTQ concatenation and post-filtering are identical on both arms and excluded from all totals.

source("figures/00_config.R")

sm <- fread(file.path(TAB_DIR, "bench_summary.tsv"))
tl <- fread(file.path(TAB_DIR, "bench_timing_long.tsv"))
N  <- nrow(sm)
sm[, `:=`(mpairs = read_pairs / 1e6, cpu_min = cpu_total_s / 60, gpu_min = gpu_total_s / 60)]
HW <- "Same L40S nodes, 32 vCPUs (16 cores), 128 GB; GPU arm + 1× L40S"

# A. Cases finished over time: all cases run back to back on one node per arm ---------------------
cum <- rbind(sm[, .(arm = "CPU", h = cumsum(cpu_total_s) / 3600)], sm[, .(arm = "GPU", h = cumsum(gpu_total_s) / 3600)])
cum <- cum[, .(h = c(0, h), n = 0:.N), by = arm]
cpu_h <- sum(sm$cpu_total_s) / 3600; gpu_h <- sum(sm$gpu_total_s) / 3600
cpu_at_gpu <- cum[arm == "CPU" & n > 0 & h <= gpu_h, .N]
p_a <- ggplot(cum, aes(h, n, colour = arm)) +
  geom_step(linewidth = 1) +
  annotate("segment", x = gpu_h, xend = gpu_h, y = 0, yend = N, colour = INK_2, linetype = "22", linewidth = 0.35) +
  annotate("text", x = gpu_h + 3, y = N, hjust = 0, vjust = 1, size = 3.1, colour = INK, lineheight = 1,
           label = sprintf("At %.1f h: GPU done with all %d,
CPU done with %d", gpu_h, N, cpu_at_gpu)) +
  annotate("text", x = cpu_h, y = N + 3, hjust = 1, size = 3.1, colour = COL_CPU, fontface = "bold",
           label = sprintf("CPU: %.1f days", cpu_h / 24)) +
  scale_colour_manual(values = arm_cols, labels = c(CPU = "CPU (GATK)", GPU = "GPU (Parabricks)"), name = NULL) +
  scale_x_continuous(breaks = seq(0, 120, 24), labels = function(x) paste0(x, " h")) +
  scale_y_continuous(breaks = seq(0, 60, 20), limits = c(0, N + 5)) +
  labs(x = "Wall-clock time, cases back to back on one node", y = "Cases finished",
       title = sprintf("All %d cases: %.1f h vs %.1f days", N, gpu_h, cpu_h / 24)) +
  theme_amp() + theme(legend.position = "top", legend.justification = "left")

# B. Where the time goes: median minutes per step -----------------------------------------------
STEP <- c(align_sort = "Align + sort", markdup = "Mark duplicates + index", index = "Mark duplicates + index",
          fq2bam = "fq2bam (align, sort, mark dup)", mutect2 = "Mutect2")
st <- tl[step %in% names(STEP)][, step_lab := STEP[step]]
st <- st[, .(seconds = sum(seconds)), by = .(sample, arm, step_lab)][, .(minutes = median(seconds) / 60), by = .(arm, step_lab)]
st[, arm := factor(toupper(arm), levels = c("GPU", "CPU"))]
st[, step_lab := factor(step_lab, levels = rev(c("Align + sort", "Mark duplicates + index",
                                                 "fq2bam (align, sort, mark dup)", "Mutect2")))]
st[, fill := paste(arm, ifelse(step_lab == "Mutect2", "call", "bam"))]
tot <- sm[, .(CPU = median(cpu_min), GPU = median(gpu_min))]
tot_b <- data.table(arm = factor(c("CPU", "GPU"), levels = c("GPU", "CPU")), minutes = c(tot$CPU, tot$GPU))
tot_b[, txt := ifelse(arm == "CPU", sprintf("%.0f min", minutes),
                      sprintf("%.1f min  (%.0f× faster)", minutes, median(sm$total_speedup)))]
step_cols <- c("CPU bam" = COL_CPU_L, "CPU call" = COL_CPU, "GPU bam" = COL_GPU_L, "GPU call" = COL_GPU)
p_b <- ggplot(st, aes(minutes, arm, fill = fill)) +
  geom_col(width = 0.6, colour = SURFACE, linewidth = 0.6, position = position_stack(reverse = TRUE)) +
  geom_text(data = tot_b, aes(minutes, arm, label = txt), inherit.aes = FALSE, hjust = -0.08, size = 3.3, colour = INK) +
  scale_fill_manual(values = step_cols, breaks = c("CPU bam", "CPU call", "GPU bam", "GPU call"),
                    labels = c("CPU alignment + dedup", "CPU Mutect2", "GPU fq2bam", "GPU Mutect2")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.35))) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(x = "Median minutes per case", y = NULL, title = "Where the time goes") +
  theme_amp() + theme(panel.grid.major.y = element_blank())

# C. Speedup vs depth ---------------------------------------------------------------------------
sp <- melt(sm, id.vars = c("sample", "mpairs"), measure.vars = c("total_speedup", "m2_speedup", "bam_speedup"),
           variable.name = "part", value.name = "speedup")
sp[, part := factor(part, levels = c("m2_speedup", "total_speedup", "bam_speedup"),
                    labels = c("Mutect2", "End to end", "Alignment + dedup"))]
lab_c <- sp[, .(mpairs = max(mpairs), speedup = median(speedup[mpairs >= quantile(mpairs, 0.5)])), by = part]
SP_COL <- c("Mutect2" = COL_GPU, "End to end" = NAVY, "Alignment + dedup" = COL_GPU_L)   # step hues as in panel B
p_c <- ggplot(sp, aes(mpairs, speedup, shape = part, fill = part)) +
  geom_point(size = 2.3, stroke = 0.4, colour = SURFACE) +
  geom_text(data = lab_c, aes(label = part), hjust = 0, nudge_x = 1.5, size = 3.2, colour = INK) +
  scale_shape_manual(values = c("Mutect2" = 21, "End to end" = 23, "Alignment + dedup" = 24), guide = "none") +
  scale_fill_manual(values = SP_COL, guide = "none") +
  scale_y_continuous(limits = c(0, NA), labels = function(x) paste0(x, "×"), expand = expansion(mult = c(0, 0.08))) +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.4))) +
  labs(x = "Read pairs (millions)", y = "GPU speedup (CPU time ÷ GPU time)", title = "Speedup grows with depth") +
  theme_amp()

cap <- sprintf("%d cases. %s. Matched pipelines: BWA-MEM2 2.2.1 vs Parabricks 4.4.0 fq2bam (-K 10000000, no -Y); GATK 4.3.0.0 Mutect2 vs Parabricks mutectcaller, genome-wide, downsampling off.", N, HW)
p <- (p_a | p_b | p_c) + plot_layout(widths = c(1, 1.1, 1)) +
  plot_annotation(title = sprintf("Parabricks runs the same somatic pipeline %.0f× faster on the same node",
                                  median(sm$total_speedup)),
                  caption = cap, theme = theme_amp())
save_fig(p, "fig_bench_runtime", 15, 5.2)
save_fig(((p_a | p_c) / p_b) + plot_layout(heights = c(1.5, 1)), "fig_bench_runtime_poster", 11.5, 6.8)   # stacked for the poster column
save_fig(p_a + labs(caption = paste(strwrap(cap, 120), collapse = "
")) + theme(plot.caption = element_text(size = 7)),
         "fig_bench_runtime_A", 7, 5)
save_fig(p_b, "fig_bench_runtime_B", 7, 4)
save_fig(p_c, "fig_bench_runtime_C", 7, 5)

fwrite(st[order(arm, step_lab)], file.path(TAB_DIR, "bench_runtime_by_step.tsv"), sep = "\t")
cat(sprintf("runtime: %d cases, median total speedup %.1fx (range %.1f-%.1f)\n", N,
            median(sm$total_speedup), min(sm$total_speedup), max(sm$total_speedup)))
