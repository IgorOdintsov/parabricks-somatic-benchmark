# Shared configuration for the AMP 2026 Parabricks figures.
# Source from the project root:  source("figures/00_config.R")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

PROJ     <- normalizePath(".")
FIG_DIR  <- file.path(PROJ, "analysis", "figures")
TAB_DIR  <- file.path(PROJ, "analysis", "tables")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)


# Colors: UCSF palette option 1 (UCSF Blue GPU, UCSF Orange CPU, UCSF Navy bands)
COL_CPU   <- "#f48024"
COL_GPU   <- "#178ccb"
COL_CPU_L <- "#f9c08f"   # lighter step of CPU hue for the second stacked segment
COL_GPU_L <- "#9fd2ee"   # light UCSF blue
INK       <- "#0b0b0b"
INK_2     <- "#52514e"
GRID      <- "#e4e3df"
SURFACE   <- "#ffffff"
NAVY      <- "#052049"   # UCSF Navy

arm_cols <- c(CPU = COL_CPU, GPU = COL_GPU)

theme_amp <- function(base_size = 11) {
  theme_minimal(base_size = base_size, base_family = "sans") %+replace%
    theme(
      text             = element_text(colour = INK),
      axis.text        = element_text(colour = INK_2),
      axis.title       = element_text(colour = INK_2),
      panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.line.x      = element_line(colour = INK_2, linewidth = 0.3),
      plot.title       = element_text(face = "bold", hjust = 0, size = rel(1.15),
                                      margin = margin(b = 6)),
      plot.subtitle    = element_text(colour = INK_2, hjust = 0, margin = margin(b = 8)),
      plot.caption     = element_text(colour = INK_2, hjust = 0, size = rel(0.8)),
      legend.position  = "top",
      legend.justification = "left",
      legend.title     = element_blank(),
      strip.text       = element_text(colour = INK, face = "bold"),
      plot.background  = element_rect(fill = SURFACE, colour = NA),
      plot.margin      = margin(10, 12, 10, 10)
    )
}

save_fig <- function(p, name, w, h) {
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 600, bg = SURFACE)  # 600 dpi: poster print
  ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), p, width = w, height = h, bg = SURFACE, device = cairo_pdf)
  invisible(p)
}

