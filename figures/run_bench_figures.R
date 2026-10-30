# Rebuild every benchmark figure and the headline table. Run from the project root after
# copying per-case outputs into results/ and running `python -I analysis/bench_summary.py .` (writes analysis/tables/bench_*.tsv):
#   Rscript figures/run_bench_figures.R
for (f in c("10_bench_runtime.R", "11_bench_concordance.R", "12_bench_headline.R",
            "14_bench_cost_bam.R")) {
  message("== ", f)
  source(file.path("figures", f), local = new.env())
}
