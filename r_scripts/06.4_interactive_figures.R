# 06.4_interactive_figures.R

# install packages
library(ggplot2)
library(patchwork)

# ggrepel for chromosome labels if requested
HAS_GGREPEL <- requireNamespace(
  "ggrepel",
  quietly = TRUE
)

# == SETTINGS ==
RESULTS_ROOT <- "/scratch/maelstr1/thesis/03_analysis/results"
# for CpG-level histograms
PROCESSED_DIR <- "/scratch/maelstr1/thesis/02_processed"

# interactive output directory
OUTPUT_DIR <- file.path(RESULTS_ROOT, "interactive_figures")
dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# scatter-plot chromosome labels
# FALSE for less visual clutter
LABEL_SCATTER_CHROMOSOMES <- FALSE

# common effect-size axis
# use the same scale to make all-sites and repeat-filtered chromosome
# effect figures directly comparable
# increase limits if any observed effect falls outside defined range

EFFECT_X_LIMITS <- c(-5, 5)
EFFECT_X_BREAKS <- seq(
  -5, 5,
  by = 1
)

# CpG histogram settings
METHYLATION_BINWIDTH <- 5
# FALSE: do not load any large BED files
# TRUE: generate CpG methylation/coverage histograms
# BED files are loaded one at a time, as needed, to reduce memory use
MAKE_CPG_HISTOGRAMS <- FALSE


# == DISPLAY LABELS ==
SAMPLE_LABELS <- c(
  "angusticeps_ANMO5228" =
    "A. angusticeps",
  
  "deppii_ANMO5280" =
    "A. deppii sample 1",
  
  "deppii_LNG1403" =
    "A. deppii sample 2",
  
  "cozumelus_ang" =
    "A. cozumelus, A. angusticeps subgenome",
  
  "cozumelus_dep" =
    "A. cozumelus, A. deppii subgenome"
)


# == COLORS ==
GROUP_COLORS <- c(
  "A. angusticeps" = "#4477AA",
  "A. deppii sample 1" = "#AA3377",
  "A. deppii sample 2" = "#CC6677",
  "A. cozumelus, A. angusticeps subgenome" = "#228833",
  "A. cozumelus, A. deppii subgenome" = "#66AA55"
)

# colors used only for comparisons between preprocessing treatments
RUN_COLORS <- c(
  "All retained sites" = "#4477AA",
  "Repeat filtered" = "#CC6677"
)


# == ITALICIZED AXIS LABELS ==
SAMPLE_AXIS_LABELS <- c(
  
  "A. angusticeps" =
    expression(
      italic("A. angusticeps")
    ),
  
  "A. deppii sample 1" =
    expression(
      italic("A. deppii") ~ "ANMO5280"
    ),
  
  "A. deppii sample 2" =
    expression(
      italic("A. deppii") ~ "LNG1403"
    ),
  
  "A. cozumelus, A. angusticeps subgenome" =
    expression(
      atop(
        italic("A. cozumelus"),
        italic("A. angusticeps") ~ "subgenome"
      )
    ),
  
  "A. cozumelus, A. deppii subgenome" =
    expression(
      atop(
        italic("A. cozumelus"),
        italic("A. deppii") ~ "subgenome"
      )
    )
)


# == GENERAL HELPER FUNCTIONS ==
save_pdf <- function(
    plot_object,
    filename,
    width = 7,
    height = 5
) {
  
  full_path <- file.path(OUTPUT_DIR, filename)
  
  ggsave(
    filename = full_path,
    plot = plot_object,
    width = width,
    height = height,
    device = "pdf"
  )
  
  cat("Saved:", full_path, "\n")
}


# convert raw sample ID to display label
sample_display_name <- function(x) {
  # convert to character explicitly in case x is stored as factor
  x <- as.character(x)

  # match raw sample IDs against the names of SAMPLE_LABELS
  idx <- match(
    x,
    names(SAMPLE_LABELS)
  )

  # stop if any raw sample ID has no display label
  if (any(is.na(idx))) {

    missing_ids <- unique(
      x[is.na(idx)]
    )

    stop(
      paste0(
        "ERROR: The following sample IDs do not have entries ",
        "in SAMPLE_LABELS:\n",
        paste(
          missing_ids,
          collapse = "\n"
        )
      )
    )
  }

  # return display labels in the same order as the input sample IDs
  unname(
    SAMPLE_LABELS[idx]
  )
}


# == LOAD ONE COMPLETED ANALYSIS ==
load_analysis <- function(
  results_root,
  run_name
) {
  # construct directories
  results_dir <- file.path(results_root, run_name)
  rds_dir <- file.path(results_dir, "rds_objects")
  
  # construct file paths
  stats_file <- file.path(rds_dir, "coordinate_independent_stats.rds")
  
  sample_summary_file <- file.path(
    results_dir, "per_sample_methylation_summary.csv"
  )
  
  chrom_summary_file <- file.path(
    results_dir, "per_chromosome_methylation_summary.csv"
  )
  
  # confirm required files exist
  required_files <- c(
    stats_file,
    sample_summary_file,
    chrom_summary_file
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
  
  
  # load stats RDS
  stats <- readRDS(
    stats_file
  )
  
  # load EDA summary tables
  sample_summary <- read.csv(sample_summary_file, stringsAsFactors = FALSE)
  chrom_summary <- read.csv(chrom_summary_file, stringsAsFactors = FALSE)
  
  # check intended files were loaded
  if (!identical(stats$run_label, run_name)) {
    
    stop(
      paste0(
        "ERROR: Requested run '",
        run_name,
        "' but stats RDS reports run_label = '",
        stats$run_label,
        "'."
      )
    )
  }
  
  # return objects
  return(
    list(
      run_name = run_name,
      results_dir = results_dir,
      stats = stats,
      sample_summary = sample_summary,
      chrom_summary = chrom_summary
    )
  )
}


# == LOAD COMPLETED ANALYSES ==
all_analysis <- load_analysis(RESULTS_ROOT, "all_retained_sites")
repeat_analysis <- load_analysis(RESULTS_ROOT, "repeat_filtered")

# == CREATE CONVENIENT SHORT NAMES ==
all_stats <- all_analysis$stats
repeat_stats <- repeat_analysis$stats

all_sample_summary <- all_analysis$sample_summary
repeat_sample_summary <- repeat_analysis$sample_summary

all_chrom_summary <- all_analysis$chrom_summary
repeat_chrom_summary <- repeat_analysis$chrom_summary


# == VERIFY LOADED ANALYSES ==
cat(
  "\nAll-retained-sites run label:",
  all_stats$run_label,
  "\n"
)

cat(
  "Repeat-filtered run label:",
  repeat_stats$run_label,
  "\n\n"
)

cat(
  "All-retained-sites comparisons:\n"
)

cat(
  paste(
    names(all_stats$comparisons),
    collapse = "\n"
  ),
  "\n\n"
)

cat(
  "Repeat-filtered comparisons:\n"
)

cat(
  paste(
    names(repeat_stats$comparisons),
    collapse = "\n"
  ),
  "\n\n"
)


# == VERIFY EXPECTED COMPARISONS ==
# deppii samples lumped together
# required_all_comparisons <- c(
#   "angusticeps_vs_cozumelus_ang_subgenome",
#   "deppii_vs_cozumelus_dep_subgenome"
# )

required_all_comparisons <- c(
    "angusticeps_vs_cozumelus_ang_subgenome",
    "deppii_ANMO5280_vs_cozumelus_dep_subgenome",
    "deppii_LNG1403_vs_cozumelus_dep_subgenome"
)

missing_all_comparisons <- setdiff(
  required_all_comparisons,
  names(all_stats$comparisons)
)


if (length(missing_all_comparisons) > 0) {
  stop(
    paste(
      "ERROR: all-retained-sites stats are missing comparison(s):",
      paste(
        missing_all_comparisons,
        collapse = ", "
      )
    )
  )
}

required_repeat_comparisons <- c(
  "deppii_ANMO5280_vs_cozumelus_dep_subgenome",
  "deppii_LNG1403_vs_cozumelus_dep_subgenome"
)

missing_repeat_comparisons <- setdiff(
  required_repeat_comparisons,
  names(repeat_stats$comparisons)
)

if (length(missing_repeat_comparisons) > 0) {
  stop(
    paste(
      "ERROR: repeat-filtered stats are missing comparison(s):",
      paste(
        missing_repeat_comparisons,
        collapse = ", "
      )
    )
  )
}


# == EXTRACT COMPLETE COMPARISON OBJECTS ==
# lumped deppii stats commented out
ang_all <- all_stats$comparisons[["angusticeps_vs_cozumelus_ang_subgenome"]]
# dep_all <- all_stats$comparisons[["deppii_vs_cozumelus_dep_subgenome"]]
dep1_all <- all_stats$comparisons[["deppii_ANMO5280_vs_cozumelus_dep_subgenome"]]
dep2_all <- all_stats$comparisons[["deppii_LNG1403_vs_cozumelus_dep_subgenome"]]

# dep_repeat <- repeat_stats$comparisons[["deppii_vs_cozumelus_dep_subgenome"]]
dep1_repeat <- repeat_stats$comparisons[["deppii_ANMO5280_vs_cozumelus_dep_subgenome"]]
dep2_repeat <- repeat_stats$comparisons[["deppii_LNG1403_vs_cozumelus_dep_subgenome"]]

# confirm the objects contain expected data
cat(
  "ang_all components:",
  paste(
    names(ang_all),
    collapse = ", "
  ),
  "\n"
)

cat(
  "dep1_all components:",
  paste(names(dep1_all), collapse = ", "),
  "\n"
)

cat(
  "dep2_all components:",
  paste(names(dep2_all), collapse = ", "),
  "\n"
)

cat(
  "dep1_repeat components:",
  paste(names(dep1_repeat), collapse = ", "),
  "\n"
)

cat(
  "dep2_repeat components:",
  paste(names(dep2_repeat), collapse = ", "),
  "\n\n"
)


# == PREPARE SAMPLE SUMMARY FOR PLOTTING ==
prepare_sample_summary <- function(
    summary_table
) {

  x <- summary_table

  # define biological/sample order using RAW SAMPLE IDs
  desired_sample_order <- c(
    "angusticeps_ANMO5228",
    "deppii_ANMO5280",
    "deppii_LNG1403",
    "cozumelus_ang",
    "cozumelus_dep"
  )

  # convert raw sample IDs to display labels
  x$display_name <- sample_display_name(
    x$sample
  )

  # determine which samples are actually present in this analysis
  # different samples available for all-retained vs repeat-filtered analyses
  present_samples <- desired_sample_order[
    desired_sample_order %in% x$sample
  ]

  # retrieve corresponding display labels
  present_display_levels <- unname(
    SAMPLE_LABELS[
      present_samples
    ]
  )

  # convert display names to an ordered factor
  x$display_name <- factor(
    x$display_name,
    levels = present_display_levels
  )

  # defensive check
  if (any(is.na(x$display_name))) {

    stop(
      paste(
        "ERROR: One or more display labels became NA during factor ordering."
      )
    )
  }

  return(x)
}


# == PLOT: MEAN METHYLATION BY DATASET ==
plot_sample_means <- function(sample_summary, subtitle_text) {
  x <- prepare_sample_summary(sample_summary)

  ggplot(x,
    aes(
      x = display_name,
      y = mean_pct_5mC,
      color = display_name
    )
  ) +
  
  # vertical line from zero to observed mean
  geom_segment(
    aes(
      xend = display_name,
      y = 0,
      yend = mean_pct_5mC
    ),
    linewidth = 0.6,
    show.legend = FALSE
  ) +
    
  # mean methylation point
  geom_point(
    size = 3.5,
    show.legend = FALSE
  ) +
    
  # exact mean above each point
  geom_text(
    aes(
      label = sprintf(
        "%.1f",
        mean_pct_5mC
      )
    ),
    vjust = -0.8,
    color = "black",
    size = 4,
    show.legend = FALSE
  ) +
    
  scale_color_manual(values = GROUP_COLORS) +
  scale_x_discrete(labels = SAMPLE_AXIS_LABELS) +
  scale_y_continuous(
    limits = c(0, 90),
    breaks = seq(
      0,
      80,
      by = 20
    ),
    expand = expansion(
      mult = c(0, 0.03)
    )
  ) +
    
  labs(
    x = NULL,
    y = "Mean 5mC (%)",
    title = "Mean methylation by dataset",
    subtitle = subtitle_text
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    axis.text.x =
      element_text(
        angle = 0,
        hjust = 0.5,
        vjust = 1,
        size = 10
      ),
    plot.margin = 
      margin(
        t = 10,
        r = 15,
        b = 20,
        l = 15
      ), 
    plot.title =
      element_text(
        size = 16
        ),
    plot.subtitle =
      element_text(
        size = 12
        )
      )
}


# == PLOT: PARENT VS HYBRID CHROMOSOME SCATTER ==
plot_chromosome_scatter <- function(
  comp,
  title_expression,
  subtitle_text
) {
  # validate comparison object
  if (is.null(comp)) {
    stop(
      "ERROR: comparison object is NULL."
    )
  }
  if (
    is.null(comp$matched) ||
    is.null(comp$summary)
  ) {
    stop(
      paste(
        "ERROR: comparison object must contain $matched and $summary."
      )
    )
  }
  m <- comp$matched
  # retrieve Spearman rho
  rho <- comp$summary$spearman_rho_across_chromosomes
  # check to prevent error in if(is.finite(rho)): argument is of length zero
  # occurs when expected field is absent
  if (length(rho) != 1) {
    rho <- NA_real_
  }
  
  # scatter plot
  p <- ggplot(m,
    aes(
      x = parent_mean,
      y = hybrid_mean
    )
  ) +
    
  # y = x means identical chromosome-wide methylation.
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.5
  ) +
    
  geom_point(
    size = 2.7,
    alpha = 0.85
  ) +
    
  # preserve equal scaling on x and y
  coord_equal() +
  
  labs(
    x = "Parental chromosome mean 5mC (%)",
    y = "Hybrid subgenome chromosome mean 5mC (%)",
    title = title_expression,
    subtitle = subtitle_text
  ) +
  
  theme_classic(base_size = 12)

  # add Spearman rho if available
  if (
    length(rho) == 1 &&
    is.finite(rho)
  ) {
    p <- p +
      annotate(
        "text",
        x = Inf,
        y = -Inf,
        label = sprintf(
          "Spearman rho = %.3f",
          rho
        ),
        hjust = 1.08,
        vjust = -0.8,
        size = 4
      )
    }
    
  # optional chromosome labels
  if (LABEL_SCATTER_CHROMOSOMES) {
    if (HAS_GGREPEL) {
      p <- p +
        ggrepel::geom_text_repel(
          aes(
            label = chrom
          ),
          size = 3,
          max.overlaps = Inf
        )
    } else {
      warning(
        paste(
          "LABEL_SCATTER_CHROMOSOMES=TRUE",
          "but ggrepel is not installed.",
          "Chromosome labels were omitted."
        )
      )
    }
  }
  return(p)
}


# == PLOT: CHROMOSOME EFFECT SIZES ==
# effect: hybrid chromosome mean - parental chromosome mean
# positive: hybrid has higher mean 5mC
#negative: hybrid has lower mean 5mC
plot_chromosome_effects <- function(
  comp,
  title_expression,
  subtitle_text
) {
  
    # confirm comparison object contains matched chromosome data
  if (
    is.null(comp) ||
    is.null(comp$matched)
  ) {
    stop(
      "ERROR: comparison object is missing $matched data."
    )
  }

  # extract matched chromosome table
  x <- comp$matched

  # extract chromosome number from labels
  x$chrom_number <- as.numeric(
    sub(
      "^Chr_",
      "",
      x$chrom
    )
  )

  # check every chromosome label was converted successfully
  if (any(is.na(x$chrom_number))) {
    warning(
      paste(
        "Some chromosome labels could not be converted",
        "to numerical chromosome values."
      )
    )
  }

  # sort chromosomes numerically
  x <- x[
    order(x$chrom_number),
    ,
    drop = FALSE
  ]

  # set chromosome plotting order
  # ggplot places the first factor level at the bottom of a discrete y-axis
  # reverse the levels so Chr_1 displays at the top
  x$chrom <- factor(
    x$chrom,
    levels = rev(x$chrom)
  )

  # build chromosome-effect plot
  # effect = hybrid chromosome mean - parental chromosome mean
  # positive value: hybrid subgenome has higher mean 5mC
  # negative value: hybrid subgenome has lower mean 5mC
  ggplot(
    x,
    aes(
      x = mean_diff_hybrid_minus_parent,
      y = chrom
    )
  ) +

  # zero indicates no parent-hybrid difference
  geom_vline(
    xintercept = 0,
    linewidth = 0.5
  ) +

  # draw line from zero to each chromosome's observed effect
  geom_segment(
    aes(
      x = 0,
      xend = mean_diff_hybrid_minus_parent,
      yend = chrom
    ),
    linewidth = 0.65
  ) +

  # mark chromosome effect
  geom_point(
    size = 2.6
  ) +

  # use the same effect scale for all comparison figures
  scale_x_continuous(
    limits = EFFECT_X_LIMITS,
    breaks = EFFECT_X_BREAKS
  ) +

  labs(
    x = paste(
      "Hybrid - Parent mean 5mC (%)"
    ),
    y = NULL,
    title = title_expression,
    subtitle = subtitle_text
    ) +

    theme_classic(
      base_size = 12
    )
}


# == ALL-RETAINED-SITES MEAN METHYLATION ==
p_mean_all <- plot_sample_means(
  all_sample_summary,
  subtitle_text = "All retained sites"
)

save_pdf(
  p_mean_all,
  "all_sites_mean_methylation_by_dataset.pdf",
  width = 9,
  height = 6
)

# == REPEAT-FILTERED MEAN METHYLATION ==
p_mean_repeat <- plot_sample_means(
  repeat_sample_summary,
  subtitle_text = "Repeat-filtered sites"
)

save_pdf(
  p_mean_repeat,
  "repeat_filtered_mean_methylation_by_dataset.pdf",
  width = 8,
  height = 6
)


# == ALL-SITES ANGUSTICEPS COMPARISON ==
p_ang_scatter_all <- plot_chromosome_scatter(
  ang_all,
  title_expression = expression(
    italic("A. angusticeps") ~
    "vs." ~
    italic("A. cozumelus") * "," ~
    italic("A. angusticeps") ~
    "subgenome"
  ),

  subtitle_text =
  "All retained sites"
)

save_pdf(
    p_ang_scatter_all,
    "all_sites_angusticeps_vs_cozumelus_ang_subgenome_scatter.pdf",
    width = 6,
    height = 6
)

p_ang_effect_all <- plot_chromosome_effects(
  ang_all,
  title_expression = expression(
    italic("A. cozumelus") * "," ~
    italic("A. angusticeps") ~
    "subgenome"
  ),

  subtitle_text ="All retained sites"
)

save_pdf(
  p_ang_effect_all,
  "all_sites_angusticeps_vs_cozumelus_ang_subgenome_effects.pdf",
  width = 7,
  height = 7
)


# == ALL-SITES DEPPII COMPARISON ==
# lumped deppii samples
# p_dep_scatter_all <- plot_chromosome_scatter(
#   dep_all,
#   title_expression = expression(
#     italic("A. deppii") ~
#     "vs." ~
#     italic("A. cozumelus") * "," ~
#     italic("A. deppii") ~
#     "subgenome"
#   ),
#   subtitle_text = "All retained sites"
# )

p_dep1_scatter_all <- plot_chromosome_scatter(
  dep1_all,
  title_expression = expression(
    italic("A. deppii") ~ "ANMO5280 vs." ~
    italic("A. cozumelus") * "," ~
    italic("A. deppii") ~ "subgenome"
  ),
  subtitle_text = "All retained sites"
)

save_pdf(
  p_dep1_scatter_all,
  "all_sites_deppii_ANMO5280_vs_cozumelus_dep_subgenome_scatter.pdf",
  width = 6,
  height = 6
)

p_dep2_scatter_all <- plot_chromosome_scatter(
  dep2_all,
  title_expression = expression(
    italic("A. deppii") ~ "LNG1403 vs." ~
    italic("A. cozumelus") * "," ~
    italic("A. deppii") ~ "subgenome"
  ),
  subtitle_text = "All retained sites"
)

save_pdf(
  p_dep2_scatter_all,
  "all_sites_deppii_LNG1403_vs_cozumelus_dep_subgenome_scatter.pdf",
  width = 6,
  height = 6
)

# p_dep_effect_all <- plot_chromosome_effects(
#   dep_all,
#   title_expression = expression(
#     italic("A. cozumelus") * "," ~
#     italic("A. deppii") ~
#     "subgenome"
#   ),
#   subtitle_text = "All retained sites"
# )

p_dep1_effect_all <- plot_chromosome_effects(
  dep1_all,
  title_expression = expression(
    italic("A. deppii") ~ "subgenome -" ~
    italic("A. deppii") ~ "ANMO5280"
  ),
  subtitle_text = "All retained sites"
)

save_pdf(
  p_dep1_effect_all,
  "all_sites_deppii_ANMO5280_vs_cozumelus_dep_subgenome_effects.pdf",
  width = 7,
  height = 7.5
)

p_dep2_effect_all <- plot_chromosome_effects(
  dep2_all,
  title_expression = expression(
    italic("A. deppii") ~ "subgenome -" ~
    italic("A. deppii") ~ "LNG1403"
  ),
  subtitle_text = "All retained sites"
)

save_pdf(
  p_dep2_effect_all,
  "all_sites_deppii_LNG1403_vs_cozumelus_dep_subgenome_effects.pdf",
  width = 7,
  height = 7.5
)


# == REPEAT-FILTERED DEPPII COMPARISON ==
# p_dep_scatter_repeat <- plot_chromosome_scatter(
#   dep_repeat,
#   title_expression = expression(
#     italic("A. deppii") ~
#     "vs." ~
#     italic("A. cozumelus") * "," ~
#     italic("A. deppii") ~
#     "subgenome"
#   ),
#   subtitle_text = "Repeat-filtered sites"
# )

p_dep1_scatter_repeat <- plot_chromosome_scatter(
  dep1_repeat,
  title_expression = expression(
    italic("A. deppii") * " ANMO5280" ~
    "vs." ~
    italic("A. cozumelus") * "," ~
    italic("A. deppii") ~
    "subgenome"
  ),
  subtitle_text = "Repeat-filtered sites"
)

save_pdf(
  p_dep1_scatter_repeat,
  "repeat_filtered_deppii_ANMO5280_vs_cozumelus_dep_subgenome_scatter.pdf",
  width = 6,
  height = 6
)

p_dep2_scatter_repeat <- plot_chromosome_scatter(
  dep2_repeat,
  title_expression = expression(
    italic("A. deppii") * " LNG1403" ~
    "vs." ~
    italic("A. cozumelus") * "," ~
    italic("A. deppii") ~
    "subgenome"
  ),
  subtitle_text = "Repeat-filtered sites"
)

save_pdf(
  p_dep2_scatter_repeat,
  "repeat_filtered_deppii_LNG1403_vs_cozumelus_dep_subgenome_scatter.pdf",
  width = 6,
  height = 6
)

# p_dep_effect_repeat <- plot_chromosome_effects(
#   dep_repeat,
#   title_expression = expression(
#     italic("A. cozumelus") * "," ~
#     italic("A. deppii") ~
#     "subgenome"
#   ),
#   subtitle_text = "Repeat-filtered sites"
# )

p_dep1_effect_repeat <- plot_chromosome_effects(
  dep1_repeat,
  title_expression = expression(
    italic("A. deppii") * " ANMO5280" ~
    "vs." ~
    italic("A. cozumelus") * "," ~
    italic("A. deppii") ~
    "subgenome"
  ),
  subtitle_text = "Repeat-filtered sites"
)

save_pdf(
  p_dep1_effect_repeat,
  "repeat_filtered_deppii_ANMO5280_vs_cozumelus_dep_subgenome_effects.pdf",
  width = 7,
  height = 7.5
)

p_dep2_effect_repeat <- plot_chromosome_effects(
  dep2_repeat,
  title_expression = expression(
    italic("A. deppii") * " LNG1403" ~
    "vs." ~
    italic("A. cozumelus") * "," ~
    italic("A. deppii") ~
    "subgenome"
  ),
  subtitle_text = "Repeat-filtered sites"
)

save_pdf(
  p_dep2_effect_repeat,
  "repeat_filtered_deppii_LNG1403_vs_cozumelus_dep_subgenome_effects.pdf",
  width = 7,
  height = 7.5
)


# CHROMOSOME-LEVEL PARENT/SUBGENOME COMPARISONS
# effects & scatter plots in four panel figure
# panel-specific formatting
# fig4_A <- p_ang_effect_all +
#     labs(
#         title = expression(
#             italic("A. angusticeps") ~ "ancestry"
#         ),
#         subtitle = NULL,
#         x = "Hybrid - parent 5mC (percentage points)"
#     ) +
#     theme(
#         axis.text.y = element_text(
#             size = 8
#         ),
#         axis.text.x = element_text(
#             size = 9
#         ),
#         axis.title.x = element_text(
#             size = 10
#         ),
#         plot.margin = margin(
#             5, 10, 5, 5
#         )
#     )
# 
# fig4_B <- p_dep_effect_all +
#     labs(
#         title = expression(
#             italic("A. deppii") ~ "ancestry"
#         ),
#         subtitle = NULL,
#         x = "Hybrid - parent 5mC (percentage points)"
#     ) +
#     theme(
#         axis.text.y = element_text(
#             size = 8
#         ),
#         axis.text.x = element_text(
#             size = 9
#         ),
#         axis.title.x = element_text(
#             size = 10
#         ),
#         plot.margin = margin(
#             5, 5, 5, 10
#         )
#     )
# 
# fig4_C <- p_ang_scatter_all +
#     labs(
#         title = NULL,
#         subtitle = NULL,
#         # shorter axis titles for combined figure
#         x = "Parental chromosome mean 5mC (%)",
#         y = "Hybrid chromosome mean 5mC (%)"
#     ) +
#     theme(
#         axis.title = element_text(
#             size = 10
#         ),
#         axis.text = element_text(
#             size = 9
#         ),
#         plot.margin = margin(
#             10, 10, 5, 5
#         )
#       )
# 
# fig4_D <- p_dep_scatter_all +
#     labs(
#         title = NULL,
#         subtitle = NULL,
#         x = "Parental chromosome mean 5mC (%)",
#         y = "Hybrid chromosome mean 5mC (%)"
#     ) +
#     theme(
#         axis.title = element_text(
#             size = 10
#         ),
#         axis.text = element_text(
#             size = 9
#         ),
#         plot.margin = margin(
#             10, 5, 5, 10
#         )
#       )
# 
# # four-panel layout
# #   A | B
# #   -----
# #   C | D
# 
# #   A = angusticeps chromosome effects
# #   B = deppii chromosome effects
# #   C = angusticeps parent-vs-hybrid scatter
# #   D = deppii parent-vs-hybrid scatter
# 
# # coord_equal() for fixed aspect ratio
# # nested patchworks to center each scatter plot independently underneath
# # corresponding chromosome-effect plot
# figure4_design <- "
# AAAAAABBBBBB
# AAAAAABBBBBB
# CCCCCCDDDDDD
# CCCCCCDDDDDD
# "
# 
# figure4 <- wrap_plots(
#     A = fig4_A,
#     B = fig4_B,
#     C = fig4_C,
#     D = fig4_D,
#     design = figure4_design
# ) +
#     plot_annotation(
#         tag_levels = "A"
#     ) &
#     theme(
#         plot.tag = element_text(
#             face = "bold",
#             size = 14
#         )
#     )
# 
# ggsave(
#     filename = file.path(
#         OUTPUT_DIR,
#         "Figure4_chromosome_level_results.pdf"
#     ),
#     plot = figure4,
#     width = 14,
#     height = 11,
#     units = "in",
#     device = "pdf"
# )


# == OVERALL MEAN METHYLATION SENSITIVITY ==
all_means <- all_sample_summary[
  ,
  c(
    "sample",
    "mean_pct_5mC"
  )
]

names(all_means)[2] <- "all_retained_sites"

repeat_means <- repeat_sample_summary[
  ,
  c(
    "sample",
    "mean_pct_5mC"
  )
]

names(repeat_means)[2] <- "repeat_filtered"

# match datasets present in both analyses
mean_sensitivity <- merge(
  all_means,
  repeat_means,
  by = "sample",
  all = FALSE
)

# add ordinary character display labels
mean_sensitivity$display_name <-
  sample_display_name(
    mean_sensitivity$sample
  )

# calculate change after repeat filtering
mean_sensitivity$change_after_repeat_filtering <-
  mean_sensitivity$repeat_filtered -
  mean_sensitivity$all_retained_sites

# Save exact values
write.csv(
  mean_sensitivity,
  file.path(
    OUTPUT_DIR,
    "mean_methylation_repeat_filter_sensitivity.csv"
  ),
  row.names = FALSE
)


# == PREPARE LONG-FORM TABLE FOR SENSITIVITY PLOT ==
mean_long <- rbind(
  data.frame(
    display_name =
      mean_sensitivity$display_name,
    run =
      "All retained sites",
    mean_pct_5mC =
      mean_sensitivity$
      all_retained_sites,
    
    stringsAsFactors = FALSE
  ),
  
  data.frame(
    display_name = mean_sensitivity$display_name,
    run = "Repeat filtered",
    mean_pct_5mC = mean_sensitivity$repeat_filtered,
    stringsAsFactors = FALSE
  )
)


# order four datasets consistently
shared_order <- c(
  "A. deppii sample 1",
  "A. deppii sample 2",
  "A. cozumelus, A. angusticeps subgenome",
  "A. cozumelus, A. deppii subgenome"
)

mean_long$display_name <- factor(
  mean_long$display_name,
  levels = shared_order
)

mean_sensitivity$display_name <- factor(
  mean_sensitivity$display_name,
  levels = shared_order
)

mean_long$run <- factor(
  mean_long$run,
  levels = c(
    "All retained sites",
    "Repeat filtered"
  )
)


# == PLOT OVERALL REPEAT-FILTER SENSITIVITY ==
p_mean_sensitivity <- ggplot(
  mean_long,
  aes(
    x = display_name,
    y = mean_pct_5mC,
    color = run
  )
) +
  # connect all-sites and repeat-filtered values for each dataset
  geom_segment(
    data = mean_sensitivity,
    aes(
      x = display_name,
      xend = display_name,
      y = all_retained_sites,
      yend = repeat_filtered
    ),
    inherit.aes = FALSE,
    linewidth = 0.6
  ) +
  geom_point(
    size = 3
  ) +
  scale_color_manual(
    values = RUN_COLORS
  ) +
  scale_x_discrete(
    labels = SAMPLE_AXIS_LABELS
  ) +
  labs(
    x = NULL,
    y = "Mean 5mC (%)",
    color = NULL,
    title = "Effect of repeat filtering on mean methylation",
    subtitle = "Datasets available in both analyses"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x =
      element_text(
        angle = 30,
        hjust = 1
      ),
    legend.position =
      "top"
  )

save_pdf(
  p_mean_sensitivity,
  "sensitivity_overall_mean_all_vs_repeat_filtered.pdf",
  width = 8.5,
  height = 5.5
)


# == DEPPII CHROMOSOME-EFFECT SENSITIVITY ==
# compare chromosome-level hybrid-parent effects before and after repeat filtering
# effect definition: hybrid chromosome mean - parental chromosome mean
# change after repeat filtering: repeat-filtered effect - all-sites effect

# positive change: effect became more positive / less negative after repeat filtering
# negative change: effect became more negative / less positive after repeat filtering

# helper function for building sensitivity table
make_dep_sensitivity <- function(
  all_comp,
  repeat_comp,
  sample_id
) {
  
  # validate comparison objects
  if (
    is.null(all_comp) ||
    is.null(all_comp$matched)
  ) {
    stop(
      paste(
        "ERROR: all-sites comparison is missing matched data for",
        sample_id
      )
    )
  }
  
  if (
    is.null(repeat_comp) ||
    is.null(repeat_comp$matched)
  ) {
    stop(
      paste(
        "ERROR: repeat-filtered comparison is missing matched data for",
        sample_id
      )
    )
  }
  
  # extract chromosome effects from all-retained-sites analysis
  all_effect <- all_comp$matched[
    ,
    c(
      "chrom",
      "mean_diff_hybrid_minus_parent"
    ),
    drop = FALSE
  ]
  
  names(all_effect)[2] <-
    "effect_all_sites"
  
  # extract chromosome effects from repeat-filtered analysis
  repeat_effect <- repeat_comp$matched[
    ,
    c(
      "chrom",
      "mean_diff_hybrid_minus_parent"
    ),
    drop = FALSE
  ]
  
  names(repeat_effect)[2] <-
    "effect_repeat_filtered"
  
  # match chromosomes represented in BOTH analyses
  sensitivity <- merge(
    all_effect,
    repeat_effect,
    by = "chrom",
    all = FALSE
  )
  
  # stop if no chromosomes were shared
  if (nrow(sensitivity) == 0) {
    stop(
      paste(
        "ERROR: No chromosomes were shared between all-sites and",
        "repeat-filtered analyses for",
        sample_id
      )
    )
  }
  
  
  # calculate change caused by repeat filtering
  sensitivity$change_after_repeat_filtering <-
    sensitivity$effect_repeat_filtered -
    sensitivity$effect_all_sites
  
  # add parental sample identity explicitly
  sensitivity$parent_sample <- sample_id
  
  # add numeric chromosome number for sorting
  sensitivity$chrom_number <- as.numeric(
    sub(
      "^Chr_",
      "",
      sensitivity$chrom
    )
  )
  
  if (any(is.na(sensitivity$chrom_number))) {
    warning(
      paste(
        "Some chromosome labels could not be converted",
        "to chromosome numbers for",
        sample_id
      )
    )
  }
  
  # sort chromosomes numerically
  sensitivity <- sensitivity[
    order(
      sensitivity$chrom_number
    ),
    ,
    drop = FALSE
  ]
  
  rownames(sensitivity) <- NULL
  
  return(sensitivity)
}

# build sensitivity tables
dep1_sensitivity <- make_dep_sensitivity(
  all_comp = dep1_all,
  repeat_comp = dep1_repeat,
  sample_id = "deppii_ANMO5280"
)

dep2_sensitivity <- make_dep_sensitivity(
  all_comp = dep2_all,
  repeat_comp = dep2_repeat,
  sample_id = "deppii_LNG1403"
)

# save sensitivity values
write.csv(
  dep1_sensitivity,
  file.path(
    OUTPUT_DIR,
    paste0(
      "deppii_ANMO5280_",
      "chromosome_effect_repeat_filter_sensitivity.csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  dep2_sensitivity,
  file.path(
    OUTPUT_DIR,
    paste0(
      "deppii_LNG1403_",
      "chromosome_effect_repeat_filter_sensitivity.csv"
    )
  ),
  row.names = FALSE
)

# combined CSV for downstream inspection
combined_dep_sensitivity <- rbind(
  dep1_sensitivity,
  dep2_sensitivity
)

write.csv(
  combined_dep_sensitivity,
  file.path(
    OUTPUT_DIR,
    "deppii_individual_chromosome_effect_repeat_filter_sensitivity.csv"
  ),
  row.names = FALSE
)

# calculate descriptive sensitivity correlations
# Do chromosomes showing stronger parent-hybrid effects in the all-sites
#  analysis also tend to show stronger effects after repeats are removed?
calculate_sensitivity_rho <- function(
  sensitivity_table
) {
  
  x <- sensitivity_table$effect_all_sites
  y <- sensitivity_table$effect_repeat_filtered
  
  valid <- is.finite(x) & is.finite(y)
  
  x <- x[valid]
  y <- y[valid]
  
  
  # need at least three observations and variation in both variables
  if (
    length(x) < 3 ||
    sd(x) == 0 ||
    sd(y) == 0
  ) {
    return(NA_real_)
  }
  
  
  cor(
    x,
    y,
    method = "spearman"
  )
}

sensitivity_rho_dep1 <- calculate_sensitivity_rho(
  dep1_sensitivity
)

sensitivity_rho_dep2 <- calculate_sensitivity_rho(
  dep2_sensitivity
)

# helper function for all-sites vs repeat-filtered scatter plots
plot_effect_sensitivity <- function(
  sensitivity_table,
  rho,
  sample_label
) {
  
  p <- ggplot(
    sensitivity_table,
    aes(
      x = effect_all_sites,
      y = effect_repeat_filtered
    )
  ) +
    
  # identity line: identical chromosome effects before and after repeat filtering
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.5
  ) +
    
  geom_point(
    size = 2.7,
    alpha = 0.85
  ) +
    
  # equal x/y scaling (makes displacement from y=x interpretable)
  coord_equal() +
    
  labs(
    x = paste(
      "Hybrid - parent effect",
      "all retained sites",
      "(%)",
      sep = "\n"
    ),
      
    y = paste(
      "Hybrid - parent effect",
      "repeat-filtered sites",
      "(%)",
      sep = "\n"
    ),
      
    title = "Repeat-filter sensitivity of chromosome effects",
      
    subtitle = sample_label
  ) +
    
  theme_classic(
    base_size = 12
  )
  
  
  # add Spearman rho, if calculable
  if (
    length(rho) == 1 &&
    is.finite(rho)
  ) {
    
    p <- p +
      annotate(
        "text",
        x = Inf,
        y = -Inf,
        label = sprintf(
          "Spearman rho = %.3f",
          rho
        ),
        hjust = 1.08,
        vjust = -0.8,
        size = 4
      )
  }
  
  
  # optional chromosome labels
  if (
    LABEL_SCATTER_CHROMOSOMES &&
    HAS_GGREPEL
  ) {
    
    p <- p +
      ggrepel::geom_text_repel(
        aes(
          label = chrom
        ),
        size = 3,
        max.overlaps = Inf
      )
  }
  
  
  return(p)
}

# create effect-concordance plots
p_effect_sensitivity_dep1 <- plot_effect_sensitivity(
  sensitivity_table = dep1_sensitivity,
  rho = sensitivity_rho_dep1,
  
  sample_label = expression(
    italic("A. deppii") ~ "ANMO5280 vs." ~
      italic("A. cozumelus") * "," ~
      italic("A. deppii") ~ "subgenome"
  )
)

p_effect_sensitivity_dep2 <- plot_effect_sensitivity(
  sensitivity_table = dep2_sensitivity,
  rho = sensitivity_rho_dep2,
  
  sample_label = expression(
    italic("A. deppii") ~ "LNG1403 vs." ~
      italic("A. cozumelus") * "," ~
      italic("A. deppii") ~ "subgenome"
  )
)

# save effect-concordance plots
save_pdf(
  p_effect_sensitivity_dep1,
  paste0(
    "sensitivity_deppii_ANMO5280_",
    "chromosome_effects_all_vs_repeat_filtered.pdf"
  ),
  width = 6,
  height = 6
)

save_pdf(
  p_effect_sensitivity_dep2,
  paste0(
    "sensitivity_deppii_LNG1403_",
    "chromosome_effects_all_vs_repeat_filtered.pdf"
  ),
  width = 6,
  height = 6
)

# helper function for chromosome-specific change after repeat filtering
plot_change_by_chromosome <- function(
  sensitivity_table,
  sample_label
) {
  
  x <- sensitivity_table
  
  
  # sort numerically, by chromosome
  x <- x[
    order(
      x$chrom_number
    ),
    ,
    drop = FALSE
  ]
  
  # ggplot places the first factor level at the bottom
  # reverse levels so Chr_1 appears at the top
  x$chrom <- factor(
    x$chrom,
    levels = rev(
      x$chrom
    )
  )
  
  ggplot(
    x,
    aes(
      x = change_after_repeat_filtering,
      y = chrom
    )
  ) +
    
    # zero = repeat filtering did not change effect
    geom_vline(
      xintercept = 0,
      linewidth = 0.5
    ) +
    
    geom_segment(
      aes(
        x = 0,
        xend =
          change_after_repeat_filtering,
        yend = chrom
      ),
      linewidth = 0.65
    ) +
    
    geom_point(
      size = 2.5
    ) +
    
    labs(
      x = paste(
        "Change in chromosome effect after repeat filtering",
        "(percentage points)",
        sep = "\n"
      ),
      
      y = NULL,
      
      title =
        "Chromosome-specific sensitivity to repeat filtering",
      
      subtitle =
        sample_label
    ) +
    
    theme_classic(
      base_size = 12
    )
}


# create chromosome-change plots
p_change_by_chr_dep1 <- plot_change_by_chromosome(
  sensitivity_table = dep1_sensitivity,
  
  sample_label = expression(
    italic("A. deppii") ~ "ANMO5280 vs." ~
      italic("A. cozumelus") * "," ~
      italic("A. deppii") ~ "subgenome"
  )
)


p_change_by_chr_dep2 <- plot_change_by_chromosome(
  sensitivity_table = dep2_sensitivity,
  
  sample_label = expression(
    italic("A. deppii") ~ "LNG1403 vs." ~
      italic("A. cozumelus") * "," ~
      italic("A. deppii") ~ "subgenome"
  )
)

# save chromosome-change plots
save_pdf(
  p_change_by_chr_dep1,
  paste0(
    "sensitivity_deppii_ANMO5280_",
    "change_after_repeat_filtering_by_chromosome.pdf"
  ),
  width = 7,
  height = 7.5
)

save_pdf(
  p_change_by_chr_dep2,
  paste0(
    "sensitivity_deppii_LNG1403_",
    "change_after_repeat_filtering_by_chromosome.pdf"
  ),
  width = 7,
  height = 7.5
)

# combined repeat-filter sensitivity figure
# A: overall mean methylation before/after repeat filtering
# B/C: ANMO5280 chromosome-specific changes and effect concordance
# D/E: LNG1403 chromosome-specific changes and effect concordance

# A
sens_A <- p_mean_sensitivity +
  
  labs(
    title ="Overall mean methylation",
    subtitle = NULL
  ) +
  
  theme(
    plot.title =
      element_text(
        size = 14
      ),
    
    axis.title =
      element_text(
        size = 11
      ),
    
    axis.text =
      element_text(
        size = 10
      ),
    
    legend.position =
      "top"
  )

# B
sens_B <- p_change_by_chr_dep1 +
  
  labs(
    title = "ANMO5280: chromosome-specific changes",
    subtitle = NULL
  ) +
  
  theme(
    plot.title =
      element_text(
        size = 12
      ),
    
    axis.title =
      element_text(
        size = 10
      ),
    
    axis.text =
      element_text(
        size = 8.5
      )
  )

# C
sens_C <- p_effect_sensitivity_dep1 +
  
  labs(
    title = "ANMO5280: chromosome-effect concordance",
    subtitle = NULL,
    
    x = paste(
      "Hybrid - parent effect",
      "all retained sites (percentage points)",
      sep = "\n"
    ),
    
    y = paste(
      "Hybrid - parent effect",
      "repeat-filtered sites (percentage points)",
      sep = "\n"
    )
  ) +
  
  theme(
    plot.title =
      element_text(
        size = 12
      ),
    
    axis.title =
      element_text(
        size = 10
      ),
    
    axis.text =
      element_text(
        size = 9
      )
  )

# D
sens_D <- p_change_by_chr_dep2 +
  
  labs(
    title = "LNG1403: chromosome-specific changes",
    subtitle = NULL
  ) +
  
  theme(
    plot.title =
      element_text(
        size = 12
      ),
    
    axis.title =
      element_text(
        size = 10
      ),
    
    axis.text =
      element_text(
        size = 8.5
      )
  )

# E
sens_E <- p_effect_sensitivity_dep2 +
  
  labs(
    title = "LNG1403: chromosome-effect concordance",
    subtitle = NULL,
    
    x = paste(
      "Hybrid - parent effect",
      "all retained sites (percentage points)",
      sep = "\n"
    ),
    
    y = paste(
      "Hybrid - parent effect",
      "repeat-filtered sites (percentage points)",
      sep = "\n"
    )
  ) +
  
  theme(
    plot.title =
      element_text(
        size = 12
      ),
    
    axis.title =
      element_text(
        size = 10
      ),
    
    axis.text =
      element_text(
        size = 9
      )
  )

# assemble five-panel sensitivity figure
#             A
#
#       B           C
#
#       D           E
#
sensitivity_figure <- (
  
  sens_A
  
) / (
  
  sens_B |
    sens_C
  
) / (
  
  sens_D |
    sens_E
  
) +
  
  plot_layout(
    heights = c(
      0.75,
      1,
      1
    )
  ) +
  
  plot_annotation(
    tag_levels = "A"
  ) &
  
  theme(
    plot.tag =
      element_text(
        face = "bold",
        size = 14
      )
  )

# save combined sensitvity figure
ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "Figure_repeat_filter_sensitivity_by_deppii_individual.pdf"
  ),
  plot = sensitivity_figure,
  width = 13,
  height = 15,
  units = "in",
  device = "pdf"
)

# print sensitivity results to log
cat(
  "\n============================================================\n"
)

cat(
  "REPEAT-FILTER SENSITIVITY SUMMARY\n"
)

cat(
  "============================================================\n"
)

# ANMO5280
cat("\nA. deppii ANMO5280:\n")

cat("Chromosome-effect change after repeat filtering:\n")

print(
  summary(
    dep1_sensitivity$
      change_after_repeat_filtering
  )
)

cat(
  "Spearman correlation of chromosome effects ",
  "between all-sites and repeat-filtered analyses: ",
  if (
    is.finite(
      sensitivity_rho_dep1
    )
  ) {
    sprintf(
      "%.3f",
      sensitivity_rho_dep1
    )
  } else {
    "NA"
  },
  "\n",
  sep = ""
)

# LNG1403
cat("\nA. deppii LNG1403:\n")

cat("Chromosome-effect change after repeat filtering:\n")

print(
  summary(
    dep2_sensitivity$
      change_after_repeat_filtering
  )
)

cat(
  "Spearman correlation of chromosome effects ",
  "between all-sites and repeat-filtered analyses: ",
  if (
    is.finite(
      sensitivity_rho_dep2
    )
  ) {
    sprintf(
      "%.3f",
      sensitivity_rho_dep2
    )
  } else {
    "NA"
  },
  "\n",
  sep = ""
)


# == PARENTAL-PATTERN SUMMARY ==
cat("\nParental-pattern summary -- all retained sites:\n")
print(all_stats$pattern_summary)
cat("\nParental-pattern summary -- repeat filtered:\n")
print(repeat_stats$pattern_summary)

# == OPTIONAL CpG-LEVEL EDA HISTOGRAM FUNCTIONS ==
# load one BED at a time to plot, save, then remove from memory
load_plot_bed <- function(
  filepath
) {
  if (!file.exists(filepath)) {
    stop(
      "ERROR: BED file not found: ",
      filepath
    )
  }

  bed <- read.table(
    filepath,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = ""
  )
  
  # confirm expected columns
  required_cols <- c(
    "chrom",
    "start",
    "end",
    "methyl",
    "strand",
    "prob",
    "cov"
  )
  
  missing_cols <- setdiff(
    required_cols,
    names(bed)
  )
  
  if (length(missing_cols) > 0) {
    stop(
      "ERROR: BED is missing required column(s): ",
      paste(
        missing_cols,
        collapse = ", "
      )
    )
  }
  
  # create clearer internal column names
  bed$pct_5mC <-
    as.numeric(
      bed$prob
    )
  
  bed$valid_cov <-
    as.numeric(
      bed$cov
    )

  # retain valid 5mC rows
  bed <- bed[
    bed$methyl == "m" &
      is.finite(
        bed$pct_5mC
      ) &
      is.finite(
        bed$valid_cov
      ),
    ,
    drop = FALSE
  ]
  return(bed)
}


# == OPTIONAL METHYLATION HISTOGRAM ==
plot_methylation_histogram <- function(
  bed,
  title_expression,
  subtitle_text
) {
  ggplot(
    bed,
    aes(
      x = pct_5mC
    )
  ) +
    geom_histogram(
      binwidth =
        METHYLATION_BINWIDTH,
      boundary = 0,
      closed = "left"
    ) +
    scale_x_continuous(
      limits = c(
        0,
        100
      ),
      breaks = seq(
        0,
        100,
        by = 20
      )
    ) +
    labs(
      x = "5mC (%)",
      y = "CpG sites",
      title =
        title_expression,
      subtitle =
        subtitle_text
    ) +
    theme_classic(
      base_size = 12
    )
}


# == OPTIONAL COVERAGE HISTOGRAM ==
plot_coverage_histogram <- function(
  bed,
  title_expression,
  subtitle_text
) {
  
  # calculate 99.5th percentile
  # extreme coverage values are excluded from the DISPLAY only
  upper <- as.numeric(
    quantile(
      bed$valid_cov,
      0.995,
      na.rm = TRUE
    )
  )
  plot_data <- bed[
    is.finite(
      bed$valid_cov
    ) &
      bed$valid_cov <= upper,
    ,
    drop = FALSE
  ]
  
  ggplot(
    plot_data,
    aes(
      x = valid_cov
    )
  ) +
    geom_histogram(
      binwidth = 1,
      boundary = 0
    ) +
    labs(
      x = "Valid coverage (<=99.5th percentile)",
      y = "CpG sites",
      title = title_expression,
      subtitle = subtitle_text
    ) +
    theme_classic(
      base_size = 12
    )
}


# == BED FILE PATHS ==
BED_PATHS <- list(
  # all retained sites
  all_retained_sites = list(
    angusticeps_ANMO5228 = file.path(
      PROCESSED_DIR,
      "angusticeps_ANMO5228",
      "plot_ready",
      "angusticeps_ANMO5228_final.bed"
    ),
    deppii_ANMO5280 = file.path(
      PROCESSED_DIR,
      "deppii_ANMO5280",
      "plot_ready",
      "deppii_ANMO5280_final.bed"
    ),
    deppii_LNG1403 = file.path(
      PROCESSED_DIR,
      "deppii_LNG1403",
      "plot_ready",
      "deppii_LNG1403_final.bed"
    ),
    cozumelus_ang = file.path(
      PROCESSED_DIR,
      "cozumelus_LNG1401",
      "plot_ready",
      "subgenome_angusticeps",
      "cozumelus_ang_subgenome.bed"
    ),
    cozumelus_dep = file.path(
      PROCESSED_DIR,
      "cozumelus_LNG1401",
      "plot_ready",
      "subgenome_deppii",
      "cozumelus_dep_subgenome.bed"
    )
  ),
  
  # repeat filtered
  repeat_filtered = list(
    # no repeat-filtered A. angusticeps BED exists
    deppii_ANMO5280 = file.path(
      PROCESSED_DIR,
      "deppii_ANMO5280",
      "plot_ready",
      "deppii_ANMO5280_final_norepeat.bed"
    ),
    deppii_LNG1403 = file.path(
      PROCESSED_DIR,
      "deppii_LNG1403",
      "plot_ready",
      "deppii_LNG1403_final_norepeat.bed"
    ),
    cozumelus_ang = file.path(
      PROCESSED_DIR,
      "cozumelus_LNG1401",
      "plot_ready",
      "subgenome_angusticeps",
      "cozumelus_ang_subgenome_norepeat.bed"
    ),
    cozumelus_dep = file.path(
      PROCESSED_DIR,
      "cozumelus_LNG1401",
      "plot_ready",
      "subgenome_deppii",
      "cozumelus_dep_subgenome_norepeat.bed"
    )
  )
)


# == TITLES FOR OPTIONAL HISTOGRAMS ==
sample_title_expression <- function(sample_id) {
  switch(
    sample_id,
    "angusticeps_ANMO5228" = expression(italic("A. angusticeps")),
    "deppii_ANMO5280" = expression(italic("A. deppii") ~ "- sample 1"),
    "deppii_LNG1403" = expression(italic("A. deppii") ~ "- sample 2"),
    "cozumelus_ang" =
      expression(italic("A. cozumelus") * "," ~ italic("A. angusticeps") ~ "subgenome"),
    "cozumelus_dep" =
      expression(italic("A. cozumelus") * "," ~ italic("A. deppii") ~ "subgenome"),
    
    stop(
      "No title expression defined for sample: ",
      sample_id
    )
  )
}


# == GENERATE HISTOGRAMS FOR ONE BED ==
make_one_bed_histograms <- function(
    filepath,
    sample_id,
    run_label
) {
  cat(
    "\nLoading BED:",
    filepath,
    "\n"
  )
  
  # load only this sample
  bed <- load_plot_bed(
    filepath
  )
  
  # determine subtitle
  subtitle_text <- if (
    run_label ==
    "all_retained_sites"
  ) {
    "All retained sites"
  } else {
    "Repeat-filtered sites"
  }

  title_expression <- sample_title_expression(sample_id)
  
  # methylation histogram
  p_meth <-
    plot_methylation_histogram(bed,
      title_expression = title_expression,
      subtitle_text = subtitle_text
    )
  
  # coverage histogram
  p_cov <-
    plot_coverage_histogram(bed,
      title_expression = title_expression,
      subtitle_text = subtitle_text
    )
  
  # filename prefix
  prefix <- if (
    run_label ==
    "all_retained_sites"
  ) {
    "all_sites"
  } else {
    "repeat_filtered"
  }
  
  # save methylation histogram
  save_pdf(
    p_meth,
    paste0(
      prefix,
      "_",
      sample_id,
      "_methylation_distribution.pdf"
    ),
    width = 6,
    height = 4.5
  )
  
  
  # save coverage histogram
  save_pdf(
    p_cov,
    paste0(
      prefix,
      "_",
      sample_id,
      "_coverage_distribution.pdf"
    ),
    width = 6,
    height = 4.5
  )
  
  
  # release memory before next BED
  rm(
    bed,
    p_meth,
    p_cov
  )
  
  invisible(
    gc()
  )
}


# == OPTIONAL HISTOGRAM EXECUTION ==
if (MAKE_CPG_HISTOGRAMS) {
  cat(
    "\nGenerating optional CpG-level histograms ",
    "one BED at a time...\n"
  )
  
  # all retained sites
  for (
    sample_id in
    names(
      BED_PATHS$
      all_retained_sites
    )
  ) {
    make_one_bed_histograms(
      BED_PATHS$
        all_retained_sites[[
          sample_id
        ]],
      sample_id,
      "all_retained_sites"
    )
  }
  
  # repeat filtered
  for (
    sample_id in
    names(
      BED_PATHS$
      repeat_filtered
    )
  ) {
    make_one_bed_histograms(
      BED_PATHS$
        repeat_filtered[[
          sample_id
        ]],
      sample_id,
      "repeat_filtered"
    )
  }
}


# == PRINT FINAL NUMERICAL SUMMARIES ==

cat("\n============================================================\n")
cat("INTERACTIVE FIGURE SCRIPT COMPLETE\n")
cat("============================================================\n")

# original all-sites statistics
cat("\nAll-retained-sites comparison summary:\n")
print(all_stats$comparison_summary)

# original repeat-filtered statistics
cat("\nRepeat-filtered comparison summary:\n")
print(repeat_stats$comparison_summary)

# overall methylation sensitivity
cat("\nOverall mean methylation sensitivity:\n")
print(mean_sensitivity)

# chromosome-effect sensitivity
cat("\nANMO5280 chromosome-effect sensitivity:\n")
print(
  summary(
    dep1_sensitivity$change_after_repeat_filtering
  )
)

cat("\nLNG1403 chromosome-effect sensitivity:\n")
print(
  summary(
    dep2_sensitivity$change_after_repeat_filtering
  )
)

# correlation between all-sites and repeat-filtered chromosome effects
cat(
  "\nANMO5280 repeat-filter sensitivity Spearman rho: ",
  sprintf("%.3f", sensitivity_rho_dep1),
  "\n",
  sep = ""
)

cat(
  "LNG1403 repeat-filter sensitivity Spearman rho: ",
  sprintf("%.3f", sensitivity_rho_dep2),
  "\n",
  sep = ""
)

# output directory
cat("\nFigures and sensitivity CSVs saved to:\n")
cat(OUTPUT_DIR, "\n")

# reminder about optional BED plots
if (!MAKE_CPG_HISTOGRAMS) {
  cat(
    "\nCpG histograms were not generated.\n",
    "Set MAKE_CPG_HISTOGRAMS <- TRUE ",
    "to create them one BED at a time.\n",
    sep = ""
  )
}
