# 06.3_figures.R
# Phase 1: Coordinate-independent results figure generation

# PURPOSE
# turn the validated EDA/statistics RDS objects into interpretable Phase 1
#  figures
# this script does not redo filtering or statistical calculations

# REQUIRED UPSTREAM FILES
#  rds_objects/coordinate_independent_eda.rds
#  rds_objects/coordinate_independent_stats.rds

# FIGURE PHILOSOPHY
#  every figure is compatible with assembly-specific coordinates
# figures either summarize one dataset independently or compare matched chromosome
#  labels


# install package
library(optparse)

# parse command-line options and construct output paths
option_list <- list(
  make_option("--output_dir", type = "character"),
  make_option("--repeat_filtered", type = "character", default = "FALSE")
)
opt <- parse_args(OptionParser(option_list = option_list))

repeat_filtered <- toupper(opt$repeat_filtered) == "TRUE"
RUN_LABEL <- if (repeat_filtered) "repeat_filtered" else "all_retained_sites"

OUTPUT_DIR <- file.path(opt$output_dir, RUN_LABEL)
RDS_DIR <- file.path(OUTPUT_DIR, "rds_objects")
FIG_DIR <- file.path(OUTPUT_DIR, "figures", "results")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# confirm prerequisite analysis objects exist, then load them
eda_file <- file.path(RDS_DIR, "coordinate_independent_eda.rds")
stats_file <- file.path(RDS_DIR, "coordinate_independent_stats.rds")

if (!file.exists(eda_file)) {
  stop(paste("ERROR: Missing", eda_file))
}
if (!file.exists(stats_file)) {
  stop(paste("ERROR: Missing", stats_file))
}

eda <- readRDS(eda_file)
stats <- readRDS(stats_file)

cat("COORDINATE-INDEPENDENT RESULT FIGURES\n")
cat("Run type:", RUN_LABEL, "\n\n")


# == FIGURE: WHOLE-DATASET MEAN METHYLATION ==
# each bar is the mean % 5mC across retained CpGs in one dataset/subgenome
# descriptive, does not require any cross-assembly coordinate match
s <- eda$sample_summary

pdf(file.path(FIG_DIR, "sample_mean_methylation.pdf"), width = 9, height = 6)
bp <- barplot(
  s$mean_pct_5mC,
  names.arg = s$sample,
  las = 2,
  ylim = c(0, 100),
  ylab = "Mean 5mC (%)",
  main = paste0("Mean 5mC by dataset (", RUN_LABEL, ")")
)

# add numerical mean above each bar for easier reading/export
text(
  bp,
  s$mean_pct_5mC + 2,
  labels = sprintf("%.1f", s$mean_pct_5mC),
  cex = 0.8
)
dev.off()


# == FUNCTION ==
# create three figures for one parent/subgenome comparison
plot_comparison <- function(comp, fig_dir, run_label) {
  # comp$matched contains one row per chromosome label found in both the parent
  # and the corresponding hybrid-derived subgenome.
  m <- comp$matched
  prefix <- comp$name


  # parent vs hybrid chromosome-mean scatter plot
  # each point is ONE matched chromosome, not one CpG
  # the y=x diagonal indicates identical chromosome-wide mean methylation
  pdf(
    file.path(fig_dir, paste0(prefix, "_chromosome_scatter.pdf")),
    width = 7,
    height = 7
  )
  plot(
    m$parent_mean,
    m$hybrid_mean,
    xlim = c(0, 100),
    ylim = c(0, 100),
    pch = 19,
    xlab = "Parental chromosome mean 5mC (%)",
    ylab = "Hybrid subgenome chromosome mean 5mC (%)",
    main = paste0(
      prefix,
      "\nchromosome-level comparison (",
      run_label,
      ")"
    )
  )
  abline(0, 1, lty = 2)

  # label points with harmonized chromosome name (Chr_1, Chr_2, ...)
  text(
    m$parent_mean,
    m$hybrid_mean,
    labels = m$chrom,
    pos = 3,
    cex = 0.65
  )
  dev.off()

  
  # chromosome-level effect-size plot
  # effect = hybrid mean 5mC - parental mean 5mC (%)
  # positive bars = higher in hybrid; negative bars = lower in hybrid
  ord <- order(m$mean_diff_hybrid_minus_parent)
  x <- m[ord, ]

  # guarantee a useful y-axis even when all observed effects are small
  lim <- max(abs(x$mean_diff_hybrid_minus_parent), 10, na.rm = TRUE)

  pdf(
    file.path(fig_dir, paste0(prefix, "_chromosome_effects.pdf")),
    width = 10,
    height = 6
  )
  bp <- barplot(
    x$mean_diff_hybrid_minus_parent,
    names.arg = x$chrom,
    las = 2,
    cex.names = 0.7,
    ylim = c(-lim * 1.15, lim * 1.15),
    ylab = "Hybrid - parent mean 5mC (percentage points)",
    main = paste0(
      prefix,
      "\nchromosome-level effect sizes (",
      run_label,
      ")"
    )
  )

  # zero = identical mean
  # +/-10 lines correspond to the default descriptive
  # parental-pattern threshold used by 06.2_stats.R
  abline(h = 0)
  abline(h = c(-10, 10), lty = 2)
  dev.off()


  # parental-pattern conservation classification
  # categories describe WHOLE-CHROMOSOME mean methylation only
  #  parent_like = within the configured threshold of parent
  #  higher_than_parent = hybrid chromosome mean above threshold
  #  lower_than_parent = hybrid chromosome mean below threshold
  #  they are NOT locus-level inheritance categories
  p <- comp$patterns
  desired <- c("parent_like", "higher_than_parent", "lower_than_parent")

  # start with zeros so categories absent from a comparison still appear
  vals <- setNames(rep(0, length(desired)), desired)
  vals[p$pattern] <- p$n_chromosomes

  pdf(
    file.path(fig_dir, paste0(prefix, "_parental_pattern.pdf")),
    width = 7,
    height = 5
  )
  bp <- barplot(
    vals,
    ylab = "Matched chromosomes",
    main = paste0(
      prefix,
      "\nparental-pattern conservation (",
      run_label,
      ")"
    ),
    las = 2
  )
  text(bp, vals, labels = vals, pos = 3, cex = 0.8)
  dev.off()
}

# generate figures only for comparisons that were valid and actually run
# two comparison figure sets for all-retained-sites analysis
# deppii comparison only for repeat-filtered analysis
for (comp in stats$comparisons) {
  plot_comparison(comp, FIG_DIR, RUN_LABEL)
}


# == FIGURE: CHROMOSOME-LEVEL SPEARMAN CORRELATION SUMMARY ==
# rho: whether relative chromosome methylation ranking is conserved
#  between each parent and its corresponding hybrid subgenome
# NA can occur when too few matched chromosomes are available or when one side
#  has no variation in chromosome means
# NA values are plotted at zero only so a placeholder bar location exists
cs <- stats$comparison_summary
rho <- cs$spearman_rho_across_chromosomes
rho_plot <- ifelse(is.na(rho), 0, rho)

pdf(file.path(FIG_DIR, "chromosome_correlation_summary.pdf"), width = 8, height = 5)
# use the comparison names present in this run rather than assuming two tests
comparison_labels <- gsub("_", " ", cs$comparison)

bp <- barplot(
  rho_plot,
  names.arg = comparison_labels,
  ylim = c(-1, 1),
  las = 2,
  ylab = "Spearman rho across matched chromosomes",
  main = paste0(
    "Chromosome-level methylation conservation (",
    RUN_LABEL,
    ")"
  )
)
abline(h = 0)
text(
  bp,
  rho_plot,
  labels = ifelse(is.na(rho), "NA", sprintf("%.3f", rho)),
  pos = ifelse(rho_plot >= 0, 3, 1),
  cex = 0.85
)
dev.off()

cat("Figure generation complete.\n")
if (repeat_filtered) {
  cat("NOTE: No angusticeps-vs-cozumelus_ang comparison figures were generated because no matched repeat-filtered angusticeps parent is available.\n")
}
cat("Figures written to:", FIG_DIR, "\n")
cat("No CpG-level volcano, MA, DMC, DMR, or inheritance figures were generated.\n")
