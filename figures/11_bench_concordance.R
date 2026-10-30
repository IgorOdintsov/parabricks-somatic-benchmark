# Benchmark concordance figures (v2): CPU vs GPU, matched pipelines, downsampling off.
# Input: analysis/tables/bench_{summary,target_calls,discordant,dp_all}.tsv from analysis/bench_summary.py.
# Target = IDT V5 exons +/-10 bp (529 genes). Jaccard = shared / union of call sets (CHROM, POS, REF, ALT after
# bcftools norm -m -any). PASS = GATK 4.3.0.0 FilterMutectCalls on both arms.

source("figures/00_config.R")

sm   <- fread(file.path(TAB_DIR, "bench_summary.tsv"))
tc   <- fread(file.path(TAB_DIR, "bench_target_calls.tsv"), colClasses = list(character = c("cpu_filter", "gpu_filter")))
disc <- fread(file.path(TAB_DIR, "bench_discordant.tsv"))
dpa  <- fread(file.path(TAB_DIR, "bench_dp_all.tsv"))
N <- nrow(sm)

# A. Concordance per case and tier ---------------------------------------------------------------
TIERS <- c(genome_raw = "Genome-wide\nall calls", genome_pass = "Genome-wide\nPASS",
           target_raw = "Panel\nall calls", target_pass = "Panel\nPASS")
cc <- rbindlist(lapply(names(TIERS), function(t)
  sm[, .(sample, tier = t, J = get(paste0(t, "_J")), n = get(paste0(t, "_n")), disc = get(paste0(t, "_disc")))]))
pool <- cc[, .(J = (sum(n) - sum(disc)) / sum(n), n = sum(n), disc = sum(disc)), by = tier]
cc[, tier := factor(TIERS[tier], levels = TIERS)]; pool[, tier := factor(TIERS[tier], levels = TIERS)]
big <- function(x) ifelse(x >= 1e6, sprintf("%.2fM", x / 1e6), comma(x))   # millions abbreviated so labels don't collide
pool[, txt := ifelse(n >= 1e6, sprintf("%.3f%%\n%s differ of %s", 100 * J, comma(disc), big(n)),
                     sprintf("%.3f%%\n%s / %s", 100 * J, comma(n - disc), comma(n)))]
y_min <- min(99.94, floor(100 * min(cc$J) * 100) / 100)   # axis follows the lowest case
y_lo <- y_min - 0.025; y_lab <- y_min - 0.015
p_a <- ggplot(cc, aes(tier, 100 * J)) +
  geom_hline(yintercept = 100, colour = GRID, linewidth = 0.5) +
  geom_jitter(width = 0.12, height = 0, size = 2.2, alpha = 0.75, colour = COL_GPU, stroke = 0) +
  geom_point(data = pool, shape = 23, size = 3.6, fill = INK, colour = SURFACE, stroke = 0.6) +
  geom_text(data = pool, aes(y = y_lab, label = txt), size = 2.9, colour = INK_2, lineheight = 0.9) +
  scale_y_continuous(limits = c(y_lo, 100.002), breaks = seq(100 - 0.02 * floor((100 - y_lo) / 0.02), 100, 0.02),
                     labels = function(x) paste0(x, "%")) +
  labs(x = NULL, y = "Concordance (Jaccard)", title = "CPU vs GPU call concordance",
       subtitle = "Dots: cases. Diamonds: pooled, with call counts below.") +
  theme_amp() + theme(panel.grid.major.x = element_blank())

# B. VAF agreement on the panel ----------------------------------------------------------------
sh <- tc[!is.na(cpu_AF) & !is.na(gpu_AF)]
r  <- cor(sh$cpu_AF, sh$gpu_AF)
p_b <- ggplot(sh, aes(cpu_AF, gpu_AF)) +
  geom_abline(slope = 1, intercept = 0, colour = GRID, linewidth = 0.6) +
  geom_point(size = 1.1, alpha = 0.35, colour = COL_GPU, stroke = 0) +
  annotate("text", x = 0.02, y = 0.98, hjust = 0, vjust = 1, size = 3.2, colour = INK,
           label = sprintf("%s shared panel calls\nr %s\nidentical VAF: %.2f%%\nmax |ΔVAF| = %.3f",
                           comma(nrow(sh)), if (r > 0.99999) "> 0.99999" else sprintf("= %.5f", r), 100 * mean(sh$cpu_AF == sh$gpu_AF), max(abs(sh$cpu_AF - sh$gpu_AF)))) +
  scale_x_continuous(labels = percent, limits = c(0, 1)) + scale_y_continuous(labels = percent, limits = c(0, 1)) +
  coord_equal() +
  labs(x = "VAF, CPU (GATK)", y = "VAF, GPU (Parabricks)", title = "Allele fractions agree") +
  theme_amp()

# C. What the residual discordant calls are: depth ---------------------------------------------
BIN_LAB <- c("1", "2", "3–4", "5–9", "10–19", "20–49", "50–99", "100–199", "200–499", "500–999", "≥1000")
dbin <- function(d) findInterval(d, c(1, 2, 3, 5, 10, 20, 50, 100, 200, 500, 1000))
dd <- rbind(data.table(set = "All calls", bin = dbin(dpa$DP_bin), n = dpa$n),
            disc[, .(set = "Discordant calls", bin = dbin(DP), n = 1)])
dd <- dd[, .(n = sum(n)), by = .(set, bin)][, frac := n / sum(n), by = set]
dd[, bin := factor(BIN_LAB[bin], levels = BIN_LAB)]
low <- disc[, mean(DP <= 3)]
pass_hi <- disc[filter == "PASS" & DP >= 10, .N]
panel_disc <- sum(sm$target_raw_disc); panel_pass_disc <- sum(sm$target_pass_disc)
panel_txt <- if (panel_disc == 0) "none on the panel" else sprintf("%d on the panel, all filtered", panel_disc)
p_c <- ggplot(dd, aes(bin, frac, fill = set)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  scale_fill_manual(values = c("All calls" = "#bdbbb4", "Discordant calls" = INK)) +
  scale_y_continuous(labels = percent, expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Read depth at the call", y = "Share of calls", title = "Residual differences are mostly single-read calls",
       subtitle = sprintf("%d genome-wide discordant calls: %.0f%% have depth ≤ 3;\n%d PASS with depth ≥ 10; %s.",
                          nrow(disc), 100 * low, pass_hi, panel_txt)) +
  theme_amp() + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major.x = element_blank())

cap <- sprintf(paste0("%d cases. Matched pipelines, GATK 4.3.0.0 / Parabricks 4.4.0, genome-wide calling, ",
                      "downsampling off; panel = IDT V5 exons ±10 bp (529 genes)."), N)
p <- (p_a | p_b | p_c) + plot_layout(widths = c(1.15, 0.9, 1.1)) +
  plot_annotation(title = if (panel_disc == 0) "Same calls, same allele fractions: 100% concordance on the panel"
                          else if (panel_pass_disc == 0) "Same calls, same allele fractions: 100% PASS concordance on the panel"
                          else "Same calls, same allele fractions",
                  caption = cap, theme = theme_amp())
save_fig(p, "fig_bench_concordance", 16, 5.6)
save_fig(p_a, "fig_bench_concordance_A", 6.5, 5)
save_fig(p_b, "fig_bench_concordance_B", 5.5, 5.5)
save_fig(p_c, "fig_bench_concordance_C", 7, 5)

# Poster flow: BAM -> calls -> VAF ------------------------------------------------------------------
# BAM concordance = share of aligned records identical in position, CIGAR, MAPQ and flags (bam_compare.sh).
bam <- sm[!is.na(bam_records)][, idn := 1 - bam_nonidentical / bam_records]
tot <- sum(bam$bam_records); nonid <- sum(bam$bam_nonidentical); one <- sum(bam$bam_only_one_arm)
bam_lo <- 100 - ceiling((100 - 100 * min(bam$idn)) * 1.25 * 1e6) / 1e6   # axis follows the least concordant case
p_bam <- ggplot(bam, aes(1, 100 * idn)) +
  geom_hline(yintercept = 100, colour = GRID, linewidth = 0.5) +
  geom_jitter(width = 0.12, height = 0, size = 2.2, alpha = 0.75, colour = COL_GPU, stroke = 0, seed = 1) +
  annotate("point", x = 1, y = 100 * (1 - nonid / tot), shape = 23, size = 3.6, fill = INK, colour = SURFACE, stroke = 0.6) +
  annotate("text", x = 1.45, y = 100 - (100 - bam_lo) * 0.08, hjust = 0, vjust = 1, size = 3.6, colour = INK, lineheight = 1.1,
           label = sprintf("%s reads compared\n%s in only one BAM\n%s differ (CIGAR only)\n%d of %d cases identical\nrecord for record",
                           if (tot >= 1e9) sprintf("%.2f billion", tot / 1e9) else comma(tot), comma(one), comma(nonid),
                           sum(bam$bam_nonidentical == 0 & bam$bam_only_one_arm == 0), nrow(bam))) +
  scale_x_continuous(limits = c(0.7, 2.6), breaks = 1, labels = "All aligned\nreads") +
  scale_y_continuous(limits = c(bam_lo, 100.000002), breaks = function(l) pretty(l, 4), labels = function(x) sprintf("%.6f%%", x)) +
  labs(x = NULL, y = "Identical records", title = "1. Same alignments",
       subtitle = "Position, CIGAR, MAPQ, flags per read") +
  theme_amp() + theme(panel.grid.major.x = element_blank())
p3 <- (p_bam | (p_a + labs(title = "2. Same variant calls")) | (p_b + labs(title = "3. Same allele fractions"))) +
  plot_layout(widths = c(0.85, 1.25, 0.95)) +
  plot_annotation(title = "CPU and GPU pipelines agree at every level: reads, calls and allele fractions",
                  caption = sprintf(paste0("%d cases. Matched pipelines, GATK 4.3.0.0 / Parabricks 4.4.0, genome-wide ",
                                           "calling, downsampling off; panel = IDT V5 exons ±10 bp (529 genes). ",
                                           "Diamonds: pooled over all cases. %s"), N,
                                    if (panel_disc > 0) sprintf("The %d panel calls that differ are all filtered (VAF < 2%%).", panel_disc) else ""),
                  theme = theme_amp())
save_fig(p3, "fig_bench_concordance3", 16, 5.6)

fwrite(pool[, .(tier = gsub("\n", " ", tier), J, n, disc)], file.path(TAB_DIR, "bench_concordance_pooled.tsv"), sep = "\t")
print(pool[, .(tier = gsub("\n", " ", tier), J = round(J, 5), n, disc)])

# Poster (48 x 32 in, half-width column): alignments + calls on top, VAF scatter separately beside the bullets
p12 <- (p_bam | (p_a + labs(title = "2. Same variant calls"))) + plot_layout(widths = c(0.8, 1.3)) +
  plot_annotation(caption = sprintf("%d cases. Diamonds: pooled over all cases.%s", N,
                                    if (panel_disc > 0) sprintf(" The %d panel calls that differ are all filtered (VAF < 2%%).", panel_disc) else ""),
                  theme = theme_amp())
save_fig(p12, "fig_bench_concordance12_poster", 11, 4.6)
save_fig(p_b + labs(title = "3. Same allele fractions"), "fig_bench_vaf_poster", 4.6, 4.6)
