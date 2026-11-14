# main_tables.R
# OUTPUT
# Table 1. Per-sample descriptive methylation summary
# Table 2. Whole-dataset parent vs hybrid effect sizes
# Table 3. Chromosome-level parent vs hybrid statistics

# import package
library(flextable)

# settings
RESULTS_ROOT <- "/scratch/maelstr1/thesis/03_analysis/results"
RUN <- "repeat_filtered"      #  "all_retained_sites" or "repeat_filtered"

RESULTS_DIR <- file.path(
  RESULTS_ROOT,
  RUN
)

STATS_DIR <- file.path(
  RESULTS_DIR,
  "statistics"
)

TABLE_DIR <- file.path(
  RESULTS_ROOT,
  "tables"
)

dir.create(
  TABLE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# read input files
sample_file <- file.path(
  RESULTS_DIR,
  "per_sample_methylation_summary.csv"
)

effect_file <- file.path(
  STATS_DIR,
  "overall_effect_sizes.csv"
)

chrom_stats_file <- file.path(
  STATS_DIR,
  "chromosome_comparison_summary.csv"
)

# confirm required files exist
required_files <- c(
  sample_file,
  effect_file,
  chrom_stats_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  
  stop(
    paste(
      "ERROR: Missing required file(s):\n",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  )
}

# load data
ps <- read.csv(
  sample_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

eff <- read.csv(
  effect_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

chrom_stats <- read.csv(
  chrom_stats_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat(
  "\nCreating tables for run:",
  RUN,
  "\n\n"
)


# DISPLAY LABELS
sample_labels <- c(
  
  "angusticeps_ANMO5228" = "A. angusticeps (ANMO5228)",
  "cozumelus_ang"        = "A. cozumelus, angusticeps-derived subgenome",
  "deppii_ANMO5280"      = "A. deppii (ANMO5280)",
  "deppii_LNG1403"       = "A. deppii (LNG1403)",
  "cozumelus_dep"        = "A. cozumelus, deppii-derived subgenome"
)

comparison_labels <- c(
  
  "angusticeps_vs_cozumelus_ang_subgenome" =
    "A. angusticeps (ANMO5228) vs. A. cozumelus, angusticeps-derived subgenome",
  
  "deppii_ANMO5280_vs_cozumelus_dep_subgenome" =
    "A. deppii (ANMO5280) vs. A. cozumelus, deppii-derived subgenome",
  
  "deppii_LNG1403_vs_cozumelus_dep_subgenome" =
    "A deppii (LNG1403) vs. A. cozumelus, deppii-derived subgenome"
)

# expected atasets for each run
if (RUN == "all_retained_sites") {
  
  expected_samples <- c(
    "angusticeps_ANMO5228",
    "cozumelus_ang",
    "deppii_ANMO5280",
    "deppii_LNG1403",
    "cozumelus_dep"
  )
  
  expected_comparisons <- c(
    "angusticeps_vs_cozumelus_ang_subgenome",
    "deppii_ANMO5280_vs_cozumelus_dep_subgenome",
    "deppii_LNG1403_vs_cozumelus_dep_subgenome"
  )
  
} else if (RUN == "repeat_filtered") {
  
  # no repeat-filtered A. angusticeps parental dataset
  expected_samples <- c(
    "cozumelus_ang",
    "deppii_ANMO5280",
    "deppii_LNG1403",
    "cozumelus_dep"
  )
  
  expected_comparisons <- c(
    "deppii_ANMO5280_vs_cozumelus_dep_subgenome",
    "deppii_LNG1403_vs_cozumelus_dep_subgenome"
  )
  
} else {
  
  stop(
    paste(
      "ERROR: Unknown RUN:",
      RUN
    )
  )
}


# validate input
missing_samples <- setdiff(
  expected_samples,
  ps$sample
)

if (length(missing_samples) > 0) {
  
  stop(
    paste(
      "ERROR: Missing expected sample(s):",
      paste(
        missing_samples,
        collapse = ", "
      )
    )
  )
}


missing_effect_comparisons <- setdiff(
  expected_comparisons,
  eff$comparison
)

if (length(missing_effect_comparisons) > 0) {
  
  stop(
    paste(
      "ERROR: overall_effect_sizes.csv is missing comparison(s):",
      paste(
        missing_effect_comparisons,
        collapse = ", "
      ),
      "\nDid you run 06.2_stats.R?"
    )
  )
}

missing_chrom_comparisons <- setdiff(
  expected_comparisons,
  chrom_stats$comparison
)

if (length(missing_chrom_comparisons) > 0) {
  
  stop(
    paste(
      "ERROR: chromosome_comparison_summary.csv is missing comparison(s):",
      paste(
        missing_chrom_comparisons,
        collapse = ", "
      )
    )
  )
}


# == TABLE 1. Per-sample descriptive methylation summary ==
ps_table <- ps[
  match(
    expected_samples,
    ps$sample
  ),
  ,
  drop = FALSE
]

table1 <- data.frame(
  
  Sample = 
    unname(
      sample_labels[
        ps_table$sample
        ]
      ),
  
  `Sites (millions)` = 
    round(
      ps_table$total_sites / 1e6, 
      1
    ),
  
  `Mean 5mC (%)` = 
    round(
      ps_table$mean_pct_5mC, 
      1
    ),
  
  `Median 5mC (%)` = 
    round(
      ps_table$median_pct_5mC, 
      1
    ),
  
  `SD (%)` = 
    round(
      ps_table$sd_pct_5mC, 
      1
    ),
  
  `Sites ≥80% (%)` = 
    round(
      ps_table$pct_sites_high_80, 
      1
    ),
  
  `Sites ≤20% (%)`= 
    round(
      ps_table$pct_sites_low_20, 
      1
    ),
  
  `Median coverage (×)` = 
    round(
      ps$median_valid_cov,
      1
    ),
  
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# == TABLE 2: Whole-dataset parent vs hybrid effect sizes ==
# positive difference: hybrid has higher mean 5mC
# negative difference: hybrid has lower mean 5mC
eff_table <- eff[
  match(
    expected_comparisons,
    eff$comparison
  ),
  ,
  drop = FALSE
]

table2 <- data.frame(
  Comparison =
    unname(
      comparison_labels[
        eff_table$comparison
      ]
    ),
  
  `Parent mean 5mC (%)` =
    round(
      eff_table$parent_mean_pct_5mC,
      1
    ),
  
  `Hybrid mean 5mC (%)` =
    round(
      eff_table$hybrid_mean_pct_5mC,
      1
    ),
  
  `Difference (hybrid - parent, pp)` =
    round(
      eff_table$difference_hybrid_minus_parent,
      1
    ),
  
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# == TABLE 3. Chromosome-level comparison statistics
cs <- chrom_stats[
  match(
    expected_comparisons,
    chrom_stats$comparison
  ),
  ,
  drop = FALSE
]

# format bootstrap CI (one column)
bootstrap_ci <- paste0(
  "[",
  sprintf(
    "%.2f",
    cs$bootstrap_95ci_lower
  ),
  ", ",
  sprintf(
    "%.2f",
    cs$bootstrap_95ci_upper
  ),
  "]"
)

# helper for readable p-values
format_p <- function(x) {
  
  out <- rep(
    NA_character_,
    length(x)
  )
  
  valid <- is.finite(x)
  
  out[
    valid & x < 0.001
  ] <- "<0.001"
  
  out[
    valid & x >= 0.001
  ] <- sprintf(
    "%.3f",
    x[
      valid & x >= 0.001
    ]
  )
  
  out
}

table3 <- data.frame(
  
  Comparison =
    unname(
      comparison_labels[
        cs$comparison
      ]
    ),
  
  `Matched chromosomes` =
    cs$matched_chromosomes,
  
  `Mean difference (pp)` =
    round(
      cs$mean_difference_hybrid_minus_parent,
      2
    ),
  
  `Bootstrap 95% CI (pp)` =
    bootstrap_ci,
  
  `Spearman rho` =
    round(
      cs$spearman_rho_across_chromosomes,
      3
    ),
  
  `Spearman p` =
    format_p(
      cs$spearman_p_exploratory
    ),
  
  `Paired Wilcoxon p` =
    format_p(
      cs$paired_wilcoxon_p_exploratory
    ),
  
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# helpers for flextable formatting
format_table <- function(
    df,
    first_col_width,
    total_width = 6.5,
    font_size = 9
) {
  
  ft <- flextable(
    df
  )
  
  # font size
  ft <- fontsize(
    ft,
    size = font_size,
    part = "all"
  )
  
  # header formatting
  ft <- bold(
    ft,
    part = "header"
  )
  
  ft <- align(
    ft,
    part = "header",
    align = "center"
  )
  
  # left-align first column
  ft <- align(
    ft,
    j = 1,
    align = "left",
    part = "all"
  )
  
  # center numerical/result columns
  if (ncol(df) > 1) {
    
    ft <- align(
      ft,
      j = 2:ncol(df),
      align = "center",
      part = "all"
    )
  }
  
  # fixed layout prevents extremely wide comparison columns
  ft <- set_table_properties(
    ft,
    layout = "fixed"
  )
  
  # first column receives more room for wrapping
  ft <- width(
    ft,
    j = 1,
    width = first_col_width
  )
  
  # divide remaining width among remaining columns
  if (ncol(df) > 1) {
    
    other_width <-
      (
        total_width -
          first_col_width
      ) /
      (
        ncol(df) - 1
      )
    
    ft <- width(
      ft,
      j = 2:ncol(df),
      width = other_width
    )
  }
  
  # padding for readability
  ft <- padding(
    ft,
    padding.top = 3,
    padding.bottom = 3,
    padding.left = 4,
    padding.right = 4,
    part = "all"
  )
  
  # vertically center cells
  ft <- valign(
    ft,
    valign = "center",
    part = "all"
  )
  
  return(ft)
}

# CREATE FORMATTED TABLES
ft1 <- format_table(
  table1,
  first_col_width = 2.0,
  total_width = 7.2,
  font_size = 8.5
)


ft2 <- format_table(
  table2,
  first_col_width = 3.1,
  total_width = 7.2,
  font_size = 9
)


ft3 <- format_table(
  table3,
  first_col_width = 2.7,
  total_width = 7.2,
  font_size = 8.5
)

# EXPORT (.docx)
table1_file <- file.path(
  TABLE_DIR,
  paste0(
    "Table1_per_sample_",
    RUN,
    ".docx"
  )
)

table2_file <- file.path(
  TABLE_DIR,
  paste0(
    "Table2_effect_sizes_",
    RUN,
    ".docx"
  )
)

table3_file <- file.path(
  TABLE_DIR,
  paste0(
    "Table3_chromosome_statistics_",
    RUN,
    ".docx"
  )
)


save_as_docx(
  ft1,
  path = table1_file
)

save_as_docx(
  ft2,
  path = table2_file
)

save_as_docx(
  ft3,
  path = table3_file
)

# EXPORT (.csv)
write.csv(
  table1,
  file.path(
    TABLE_DIR,
    paste0(
      "Table1_per_sample_",
      RUN,
      ".csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  table2,
  file.path(
    TABLE_DIR,
    paste0(
      "Table2_effect_sizes_",
      RUN,
      ".csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  table3,
  file.path(
    TABLE_DIR,
    paste0(
      "Table3_chromosome_statistics_",
      RUN,
      ".csv"
    )
  ),
  row.names = FALSE
)

# print to log
cat(
  "\n============================================================\n"
)

cat(
  "TABLE GENERATION COMPLETE\n"
)

cat(
  "============================================================\n\n"
)


cat(
  "TABLE 1: Per-sample summary\n\n"
)

print(
  table1
)


cat(
  "\nTABLE 2: Whole-dataset effect sizes\n\n"
)

print(
  table2
)


cat(
  "\nTABLE 3: Chromosome-level statistics\n\n"
)

print(
  table3
)


cat(
  "\nTables saved to:\n",
  TABLE_DIR,
  "\n\n"
)


cat(
  "Word files:\n",
  "  ",
  table1_file,
  "\n",
  "  ",
  table2_file,
  "\n",
  "  ",
  table3_file,
  "\n",
  sep = ""
)