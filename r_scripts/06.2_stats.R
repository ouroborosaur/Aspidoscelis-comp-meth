# 06.2_stats.R
# coordinate-independent chromosome-level stats

# compare each sexual parent with the corresponding ancestry-annotated
#  cozumelus subgenome (without assuming CpG-coordinate homology)
#  A. angusticeps ANMO5228 vs cozumelus angusticeps-derived subgenome
#  A. deppii ANMO5280 vs cozumelus deppii-derived subgenome
#  A. deppii LNG1403 vs cozumelus deppii-derived subgenome

# STATISTICS / SUMMARIES
#  whole-dataset descriptive effect sizes
#  chromosome mean hybrid-minus-parent differences
#  Spearman correlation across matched chromosome means
#  paired Wilcoxon test across chromosome means (EXPLORATORY)
#  bootstrap 95% CI for mean chromosome-level difference
#  descriptive parental-pattern conservation categories


# install package
library(optparse)

# parse command-line options
option_list <- list(
  make_option("--output_dir", type = "character"),
  make_option("--repeat_filtered", type = "character", default = "FALSE"),

  # number of chromosome-resampling bootstrap iterations
  make_option("--bootstrap_reps", type = "integer", default = 5000),

  # difference in percentage points used to call a chromosome "parent_like"
  # Example, threshold = 10:
  #  hybrid-parent = +7  -> parent_like
  #  hybrid-parent = +15 -> higher_than_parent
  #  hybrid-parent = -15 -> lower_than_parent
  make_option("--pattern_threshold", type = "double", default = 10)
)
opt <- parse_args(OptionParser(option_list = option_list))

repeat_filtered <- toupper(opt$repeat_filtered) == "TRUE"
RUN_LABEL <- if (repeat_filtered) "repeat_filtered" else "all_retained_sites"

OUTPUT_DIR <- file.path(opt$output_dir, RUN_LABEL)
RDS_DIR <- file.path(OUTPUT_DIR, "rds_objects")
STATS_DIR <- file.path(OUTPUT_DIR, "statistics")
dir.create(STATS_DIR, recursive = TRUE, showWarnings = FALSE)

# load validated EDA object
eda_file <- file.path(RDS_DIR, "coordinate_independent_eda.rds")
if (!file.exists(eda_file)) {
  stop(paste("ERROR: EDA RDS not found:", eda_file))
}

eda <- readRDS(eda_file)
chrom <- eda$chrom_summary
sample_summary <- eda$sample_summary

cat("COORDINATE-INDEPENDENT STATISTICAL SUMMARY\n")
cat("Run type:", RUN_LABEL, "\n")
cat("Matched unit: chromosome label, NOT CpG coordinate.\n\n")

# == BUILD ONE PARENTAL CHROMOSOME TABLE PER COMPARISON ==
# angusticeps: one parental sample
# its chromosome means are used directly
# deppii: two parental samples
# analyze each individual separately
make_parent_chrom <- function(samples) {
  x <- chrom[chrom$sample %in% samples, ]

  if (length(samples) == 1) {
    out <- x[, c("chrom", "mean_pct_5mC", "median_pct_5mC", "n_sites")]
    names(out) <- c(
      "chrom", "parent_mean", "parent_median", "parent_n_sites"
    )
    return(out)
  }

  # equal weight to each biological parental sample for mean/median methylation
  mean_tab <- aggregate(mean_pct_5mC ~ chrom, data = x, FUN = mean)
  med_tab <- aggregate(median_pct_5mC ~ chrom, data = x, FUN = mean)

  # n_sites is summed only as descriptive information about the number of
  # retained CpGs represented across the parental samples
  n_tab <- aggregate(n_sites ~ chrom, data = x, FUN = sum)

  out <- Reduce(
    function(a, b) merge(a, b, by = "chrom", all = TRUE),
      list(mean_tab, med_tab, n_tab)
  )
  names(out) <- c("chrom", "parent_mean", "parent_median", "parent_n_sites")
  return(out)
}

# build corresponding hybrid-subgenome chromosome table
make_hybrid_chrom <- function(sample) {
  x <- chrom[
    chrom$sample == sample,
    c("chrom", "mean_pct_5mC", "median_pct_5mC", "n_sites")
  ]
  names(x) <- c("chrom", "hybrid_mean", "hybrid_median", "hybrid_n_sites")
  return(x)
}

# == BOOTSTRAP HELPER ==
# mean chromosome-level hybrid-minus-parent difference
# bootstrap resamples matched CHROMOSOMES with replacement
# "How stable is the mean chromosome-level difference across the observed set
#   of chromosomes?"
# does not estimate uncertainty among biological individuals
# (number of biological samples is limited and chromosomes are not replicates)
bootstrap_mean_diff <- function(diffs, B = 5000, seed = 20261021) {
  diffs <- diffs[is.finite(diffs)]

  # confidence interval cannot be meaningfully generated from <2 chromosomes
  if (length(diffs) < 2) {
    return(c(lower = NA_real_, upper = NA_real_))
  }

  # fixed seed for reproducibility
  set.seed(seed)

  boot <- replicate(
    B,
    mean(sample(diffs, length(diffs), replace = TRUE))
  )

  return(quantile(
    boot,
    c(0.025, 0.975),
    na.rm = TRUE,
    names = FALSE
  ))
}

# run one ancestry-matched chromosome comparison
run_chrom_comparison <- function(parent_samples, hybrid_sample, comparison_name) {

  # create one parental and one hybrid row per harmonized chromosome label
  p <- make_parent_chrom(parent_samples)
  h <- make_hybrid_chrom(hybrid_sample)

  # inner merge (all=FALSE) keeps only chromosome labels present in BOTH sides
  # chromosome-name matching only; no start/end coordinates used
  m <- merge(p, h, by = "chrom", all = FALSE)
  m <- m[complete.cases(m$parent_mean, m$hybrid_mean), ]

  # effect size in methylation %
  # positive = higher mean 5mC in hybrid-derived chromosome
  # negative = lower mean 5mC in the hybrid-derived chromosome
  m$mean_diff_hybrid_minus_parent <- m$hybrid_mean - m$parent_mean
  m$median_diff_hybrid_minus_parent <- m$hybrid_median - m$parent_median

  # descriptive parental-pattern classification
  # "parent_like" means similar at the whole-chromosome mean level only
  threshold <- opt$pattern_threshold
  m$parental_pattern <- ifelse(
    abs(m$mean_diff_hybrid_minus_parent) <= threshold,
    "parent_like",
    ifelse(
      m$mean_diff_hybrid_minus_parent > threshold,
      "higher_than_parent",
      "lower_than_parent"
    )
  )

  n_chr <- nrow(m)

  # spearman correlation across matched chromosomes
  # Do chromosomes that are relatively highly methylated in the parent also
  #  tend to be relatively highly methylated in the corresponding subgenome?
  # spearman uses ranks, does not require a linear relationship
  rho <- NA_real_
  rho_p <- NA_real_
  if (n_chr >= 3 && sd(m$parent_mean) > 0 && sd(m$hybrid_mean) > 0) {
    sp <- suppressWarnings(cor.test(
      m$parent_mean,
      m$hybrid_mean,
      method = "spearman",
      exact = FALSE
    ))
    rho <- unname(sp$estimate)
    rho_p <- sp$p.value
  }

  # exploratory paired Wilcoxon signed-rank test across chromosomes
  # pairing is by chromosome label
  # tests whether chromosome-level hybrid means tend to be shifted relative to
  #  their corresponding parental means
  
  # IMPORTANT: chromosomes are not biological replicates
  # p-value is saved with the suffix "_exploratory"
  # interpret as supportive/descriptive rather than population-level inference
  wilcox_v <- NA_real_
  wilcox_p <- NA_real_
  if (n_chr >= 2 && any(m$mean_diff_hybrid_minus_parent != 0)) {
    wt <- suppressWarnings(wilcox.test(
      m$hybrid_mean,
      m$parent_mean,
      paired = TRUE,
      exact = FALSE
    ))
      wilcox_v <- unname(wt$statistic)
      wilcox_p <- wt$p.value
  }

  # bootstrap interval for mean chromosome-level effect
  ci <- bootstrap_mean_diff(
    m$mean_diff_hybrid_minus_parent,
    opt$bootstrap_reps
  )

  # one-row summary for this parent/subgenome comparison
  summary <- data.frame(
    comparison = comparison_name,
    parent_samples = paste(parent_samples, collapse = ";"),
    hybrid_sample = hybrid_sample,
    matched_chromosomes = n_chr,
    mean_parent_pct_5mC = mean(m$parent_mean, na.rm = TRUE),
    mean_hybrid_pct_5mC = mean(m$hybrid_mean, na.rm = TRUE),
    mean_difference_hybrid_minus_parent = mean(
      m$mean_diff_hybrid_minus_parent,
      na.rm = TRUE
    ),
    median_difference_hybrid_minus_parent = median(
      m$mean_diff_hybrid_minus_parent,
      na.rm = TRUE
    ),
    bootstrap_95ci_lower = ci[1],
    bootstrap_95ci_upper = ci[2],
    spearman_rho_across_chromosomes = rho,
    spearman_p_exploratory = rho_p,
    paired_wilcoxon_V = wilcox_v,
    paired_wilcoxon_p_exploratory = wilcox_p,
    stringsAsFactors = FALSE
  )

  # count how many chromosomes fall into each descriptive pattern category
  pattern_counts <- as.data.frame(
    table(m$parental_pattern),
    stringsAsFactors = FALSE
  )
  names(pattern_counts) <- c("pattern", "n_chromosomes")
  pattern_counts$comparison <- comparison_name
  pattern_counts$pct_chromosomes <-
    100 * pattern_counts$n_chromosomes / n_chr

  # return detailed & summarized results for figure generation & auditing
  return(list(
    name = comparison_name,
    matched = m,
    summary = summary,
    patterns = pattern_counts
  ))
}

# decide which ancestry-matched comparisons are valid for this run
# ALL RETAINED SITES:
#   A. angusticeps parent vs cozumelus angusticeps-derived subgenome
#   B. two deppii parents vs cozumelus deppii-derived subgenome
#
# REPEAT FILTERED:
#   A. angusticeps parent is unavailable
#      cozumelus_ang was retained in EDA for descriptive summaries only
#   B. repeat-filtered files exist for both deppii parents & corresponding 
#      corresponding cozumelus deppii-derived subgenome
comparisons <- list()

if (!repeat_filtered) {
  comparisons$angusticeps_vs_cozumelus_ang_subgenome <- run_chrom_comparison(
    parent_samples = "angusticeps_ANMO5228",
    hybrid_sample = "cozumelus_ang",
    comparison_name = "angusticeps_vs_cozumelus_ang_subgenome"
  )
} else {
  cat(
    "Skipping angusticeps vs cozumelus_ang statistical comparison:\n",
    "  no repeat-filtered angusticeps parental dataset is available.\n",
    "  cozumelus_ang remains descriptive-only in this run.\n\n",
    sep = ""
  )
}


# both samples in one comparison
# comparisons$deppii_vs_cozumelus_dep_subgenome <- run_chrom_comparison(
#  parent_samples = c("deppii_ANMO5280", "deppii_LNG1403"),
#  hybrid_sample = "cozumelus_dep",
#  comparison_name = "deppii_vs_cozumelus_dep_subgenome"
# )

comparisons$deppii_ANMO5280_vs_cozumelus_dep_subgenome <-
  run_chrom_comparison(
    parent_samples = "deppii_ANMO5280",
    hybrid_sample = "cozumelus_dep",
    comparison_name =
      "deppii_ANMO5280_vs_cozumelus_dep_subgenome"
  )

comparisons$deppii_LNG1403_vs_cozumelus_dep_subgenome <-
  run_chrom_comparison(
    parent_samples = "deppii_LNG1403",
    hybrid_sample = "cozumelus_dep",
    comparison_name =
      "deppii_LNG1403_vs_cozumelus_dep_subgenome"
  )

# whole-dataset descriptive effect sizes for VALID comparisons only
# These effect sizes summarize complete datasets/subgenomes. No CpG-level
# p-value is calculated because CpGs are not treated as biological replicates.
get_sample_mean <- function(sample) {
  value <- sample_summary$mean_pct_5mC[sample_summary$sample == sample]
  if (length(value) != 1) {
    stop(paste("ERROR: Expected exactly one sample summary for", sample))
  }
  value
}

overall_rows <- list()

if (!repeat_filtered) {
  overall_rows[["angusticeps"]] <- data.frame(
    comparison = "angusticeps_vs_cozumelus_ang_subgenome",
    parent_mean_pct_5mC = get_sample_mean("angusticeps_ANMO5228"),
    hybrid_mean_pct_5mC = get_sample_mean("cozumelus_ang"),
    stringsAsFactors = FALSE
  )
}

# both deppii samples in one analysis
# overall_rows[["deppii"]] <- data.frame(
#     comparison = "deppii_vs_cozumelus_dep_subgenome",
#     parent_mean_pct_5mC = mean(c(
#         get_sample_mean("deppii_LNG1403"),
#         get_sample_mean("deppii_ANMO5280")
#     )),
#     hybrid_mean_pct_5mC = get_sample_mean("cozumelus_dep"),
#     stringsAsFactors = FALSE
# )

overall_rows[["deppii_ANMO5280"]] <- data.frame(
  comparison = "deppii_ANMO5280_vs_cozumelus_dep_subgenome",
  parent_mean_pct_5mC = get_sample_mean("deppii_ANMO5280"),
  hybrid_mean_pct_5mC = get_sample_mean("cozumelus_dep"),
  stringsAsFactors = FALSE
)

overall_rows[["deppii_LNG1403"]] <- data.frame(
  comparison = "deppii_LNG1403_vs_cozumelus_dep_subgenome",
  parent_mean_pct_5mC = get_sample_mean("deppii_LNG1403"),
  hybrid_mean_pct_5mC = get_sample_mean("cozumelus_dep"),
  stringsAsFactors = FALSE
)

overall_effects <- do.call(rbind, overall_rows)
rownames(overall_effects) <- NULL
overall_effects$difference_hybrid_minus_parent <-
  overall_effects$hybrid_mean_pct_5mC - overall_effects$parent_mean_pct_5mC

# combine results dynamically
# Using lapply() over the comparison list means the output automatically contains
# two comparisons for all-retained-sites runs and one for repeat-filtered runs.
comparison_summary <- do.call(
  rbind,
  lapply(comparisons, function(x) x$summary)
)
rownames(comparison_summary) <- NULL

matched_chromosomes <- do.call(
  rbind,
  lapply(comparisons, function(x) {
    transform(x$matched, comparison = x$name)
  })
)
rownames(matched_chromosomes) <- NULL

pattern_summary <- do.call(
  rbind,
  lapply(comparisons, function(x) x$patterns)
)
rownames(pattern_summary) <- NULL

# small design table travels with the results and makes intentional omissions
# explicit when CSVs are revisited

# analysis lumping deppii samples
# analysis_design <- data.frame(
#   run_label = RUN_LABEL,
#   comparison = c(
#     "angusticeps_vs_cozumelus_ang_subgenome",
#     "deppii_vs_cozumelus_dep_subgenome"
#   ),
#   status = c(
#     if (repeat_filtered) "not_run_no_matched_repeat_filtered_parent" else "run",
#     "run"
#   ),
#   reason = c(
#     if (repeat_filtered) {
#       "Repeat-filtered angusticeps parent unavailable; cozumelus_ang is descriptive-only."
#     } else {
#       "Matched all-retained-sites parent and hybrid-subgenome datasets available."
#     },
#       "Matched parent and hybrid-subgenome preprocessing available for this run type."
#   ),
#   stringsAsFactors = FALSE
# )

analysis_design <- data.frame(
  run_label = RUN_LABEL,

  comparison = c(
    "angusticeps_vs_cozumelus_ang_subgenome",
    "deppii_ANMO5280_vs_cozumelus_dep_subgenome",
    "deppii_LNG1403_vs_cozumelus_dep_subgenome"
  ),

  status = c(
    if (repeat_filtered)
      "not_run_no_matched_repeat_filtered_parent"
    else
      "run",
      "run",
      "run"
  ),

  reason = c(
    if (repeat_filtered) {
      paste(
        "Repeat-filtered angusticeps parent unavailable;",
        "cozumelus_ang is descriptive-only."
      )
    } else {
      paste(
        "Matched all-retained-sites parent and",
        "hybrid-subgenome datasets available."
      )
    },

    "Matched ANMO5280 and cozumelus_dep preprocessing available.",
    "Matched LNG1403 and cozumelus_dep preprocessing available."
  ),

  stringsAsFactors = FALSE
)

# write human-readable results
write.csv(
  overall_effects,
  file.path(STATS_DIR, "overall_effect_sizes.csv"),
  row.names = FALSE
)
write.csv(
  comparison_summary,
  file.path(STATS_DIR, "chromosome_comparison_summary.csv"),
  row.names = FALSE
)
write.csv(
  matched_chromosomes,
  file.path(STATS_DIR, "matched_chromosome_effects.csv"),
  row.names = FALSE
)
write.csv(
  pattern_summary,
  file.path(STATS_DIR, "parental_pattern_by_chromosome.csv"),
  row.names = FALSE
)
write.csv(
  analysis_design,
  file.path(STATS_DIR, "analysis_design.csv"),
  row.names = FALSE
)

# save complete results for figure generation and future review
# The interpretation_note is intentionally stored inside the RDS object so the
# statistical limitations travel with the results even if the script is revisited
# much later.
stats_results <- list(
  run_label = RUN_LABEL,
  repeat_filtered = repeat_filtered,
  analysis_design = analysis_design,
  overall_effects = overall_effects,
  comparison_summary = comparison_summary,
  matched_chromosomes = matched_chromosomes,
  pattern_summary = pattern_summary,
  comparisons = comparisons,
  interpretation_note = paste(
    "Analyses are coordinate-independent. Chromosome labels are matched by annotated ancestry/name.",
    "No CpG sites or fixed genomic windows are assumed homologous across assemblies.",
    "Spearman and paired Wilcoxon p-values use chromosomes as matched observational units and are exploratory,",
    "not substitutes for biological replication. Bootstrap CIs resample chromosomes and describe across-chromosome consistency.",
    if (repeat_filtered) {
      paste(
        "The angusticeps comparison is intentionally absent because",
        "no repeat-filtered angusticeps parent is available.",
        "The two deppii individuals are compared independently with",
        "the cozumelus deppii-derived subgenome."
      )
    } else {
      paste(
        "Three parent-subgenome comparisons are available:",
        "one angusticeps comparison and two independent deppii comparisons."
      )
    }
  )
)
saveRDS(
  stats_results,
  file.path(RDS_DIR, "coordinate_independent_stats.rds")
)

# print a compact summary into the SLURM/log output for quick inspection
cat("Statistical summary complete.\n")
cat("Matched chromosomes:\n")
print(comparison_summary[, c(
  "comparison",
  "matched_chromosomes",
  "mean_difference_hybrid_minus_parent",
  "spearman_rho_across_chromosomes",
  "paired_wilcoxon_p_exploratory"
)])
cat("\nIMPORTANT: No CpG-level DMC/DMR tests were performed.\n")
