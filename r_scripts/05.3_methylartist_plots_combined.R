# 05.3_methylartist_combined_plots.R
# create combined figures from Methylartist segmeth TSV files
# each row in a segmeth file represents one 100-kb genomic window

# 1. read all segmeth TSV files
# 2. identifiy sample-specific Methylartist columns
# 3. calculate informative methylation calls and 5mC percentage
# 4. apply the same minimum-call threshold to every sample
# 5. retain the appropriate genomic categories for each sample
# 6. combine all retained windows into one plotting table
# 7. create combined strip and violin plots
# 8. write a numerical summary of the plotted data
# 9. save PDF & PNG files

# pct_5mC is a window-level methylation percentage
# each point in the strip plot represents one retained 100-kb window

# install package
library(ggplot2)

# == SETTINGS ==
# input segmeth TSV files
SEGMETH_FILES <- c(
  angusticeps = "/scratch/maelstr1/thesis/02_processed/angusticeps_ANMO5228/methylartist/angusticeps_ANMO5228.segmeth.tsv",
  cozumelus   = "/scratch/maelstr1/thesis/02_processed/cozumelus_LNG1401/methylartist/cozumelus_LNG1401.segmeth.tsv",
  deppii_1    = "/scratch/maelstr1/thesis/02_processed/deppii_ANMO5280/methylartist/deppii_ANMO5280.segmeth.tsv",
  deppii_2    = "/scratch/maelstr1/thesis/02_processed/deppii_LNG1403/methylartist/deppii_LNG1403.segmeth.tsv"
)

# output directory
OUTPUT_DIR <- "/scratch/maelstr1/thesis/02_processed/combined_methylartist_plots"

# filtering:
# minimum number of informative methylation calls per 100-kb window

# informative calls = methylated calls + unmethylated calls
# no-calls are excluded
MIN_CALLS <- 100

# plot settings
POINT_ALPHA  <- 0.25
POINT_SIZE   <- 0.8
VIOLIN_ALPHA <- 0.60

# common biologically interpretable methylation scale
Y_MIN <- 0
Y_MAX <- 100

# colors: Okabe-Ito palette
# pink = angusticeps (maternal)
# blue = deppii (paternal)
ANG_COLOR <- "#CC79A7"
DEP_COLOR <- "#56B4E9"

PLOT_COLORS <- c(
  "angusticeps"   = ANG_COLOR,
  "cozumelus_ang" = ANG_COLOR,
  "cozumelus_dep" = DEP_COLOR,
  "deppii_1"      = DEP_COLOR,
  "deppii_2"      = DEP_COLOR
)


# display labels
PLOT_LABELS <- c(
  "angusticeps"   = expression(italic("A. angusticeps")),
  "cozumelus_ang" = expression(italic("angusticeps") * "-derived"),
  "cozumelus_dep" = expression(italic("deppii") * "-derived"),
  "deppii_1"      = "ANMO5280",
  "deppii_2"      = "LNG1403"
)

# == HELPER FUNCTIONS ==
# FIND ONE METHYLARTIST COLUMN
# Methylartist includes the BAM/sample name in its output column names
# find desired column using the invariant suffix
find_one_column <- function(pattern, column_names, description) {
  
  hit <- grep(
    pattern,
    column_names,
    value = TRUE
  )
  
  if (length(hit) == 0) {
    stop(
      "Could not find ",
      description,
      " column using pattern: ",
      pattern
    )
  }
  
  if (length(hit) > 1) {
    stop(
      "Found more than one possible ",
      description,
      " column:\n  ",
      paste(hit, collapse = "\n  ")
    )
  }
  
  hit
}

# READ & PROCESS ONE SEGMETH FILE
read_segmeth <- function(file, sample_id) {
  
  cat(
    "\n============================================================\n",
    "Processing sample: ", sample_id, "\n",
    "File: ", file, "\n",
    "============================================================\n",
    sep = ""
  )
  
  # check input
  if (!file.exists(file)) {
    stop(
      "Input segmeth file does not exist: ",
      file
    )
  }
  
  # read segmeth table
  seg <- read.delim(
    file,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  cat(
    "Rows/windows loaded: ",
    nrow(seg),
    "\n",
    sep = ""
  )
  
  cat(
    "Columns loaded: ",
    ncol(seg),
    "\n",
    sep = ""
  )
  
  # confirm required segment columns
  required_segment_columns <- c(
    "seg_id",
    "seg_chrom",
    "seg_start",
    "seg_end",
    "seg_name",
    "seg_strand"
  )
  
  missing_segment_columns <- setdiff(
    required_segment_columns,
    names(seg)
  )
  
  if (length(missing_segment_columns) > 0) {
    stop(
      "Missing expected segment columns: ",
      paste(
        missing_segment_columns,
        collapse = ", "
      )
    )
  }
  
  # find sample-specific Methylartist columns
  meth_col <- find_one_column(
    "_m_meth_calls$",
    names(seg),
    "methylated-call"
  )
  
  unmeth_col <- find_one_column(
    "_m_unmeth_calls$",
    names(seg),
    "unmethylated-call"
  )
  
  no_call_col <- find_one_column(
    "_m_no_calls$",
    names(seg),
    "no-call"
  )
  
  methfrac_col <- find_one_column(
    "_m_methfrac$",
    names(seg),
    "methylation-fraction"
  )
  
  readcount_col <- find_one_column(
    "_m_readcount$",
    names(seg),
    "read-count"
  )
  
  motifcount_col <- find_one_column(
    "_m_motifcount$",
    names(seg),
    "motif-count"
  )
  
  cat("\nDetected Methylartist columns:\n")
  cat("  methylated calls:   ", meth_col, "\n", sep = "")
  cat("  unmethylated calls: ", unmeth_col, "\n", sep = "")
  cat("  no calls:           ", no_call_col, "\n", sep = "")
  cat("  methylation frac:   ", methfrac_col, "\n", sep = "")
  cat("  read count:         ", readcount_col, "\n", sep = "")
  cat("  motif count:        ", motifcount_col, "\n", sep = "")
  
  # create readable analysis columns
  seg$meth_calls <- as.numeric(
    seg[[meth_col]]
  )
  
  seg$unmeth_calls <- as.numeric(
    seg[[unmeth_col]]
  )
  
  seg$no_calls <- as.numeric(
    seg[[no_call_col]]
  )
  
  seg$methfrac <- as.numeric(
    seg[[methfrac_col]]
  )
  
  seg$readcount <- as.numeric(
    seg[[readcount_col]]
  )
  
  seg$motifcount <- as.numeric(
    seg[[motifcount_col]]
  )
  
  
  # informative calls exclude no-calls
  seg$total_calls <-
    seg$meth_calls +
    seg$unmeth_calls
  
  # Methylartist reports methfrac on a 0-1 scale
  # convert to percent for plotting
  seg$pct_5mC <-
    seg$methfrac * 100
  
  # select appropriate genomic categories
  if (sample_id == "cozumelus") {
    
    # A. cozumelus: retain the two annotated parental-derived subgenomes
    # unplaced sequence is excluded
    keep_groups <- c(
      "ang_subgenome",
      "dep_subgenome"
    )
    
  } else {
    
    # Sexual species: retain the genome category
    keep_groups <- "genome"
  }
  
  # filter windows
  n_before <- nrow(seg)
  
  seg <- seg[
    !is.na(seg$pct_5mC) &
      seg$total_calls >= MIN_CALLS &
      seg$seg_name %in% keep_groups,
  ]
  
  n_after <- nrow(seg)
  
  cat(
    "\nWindows retained after filtering: ",
    n_after,
    " of ",
    n_before,
    "\n",
    sep = ""
  )
  
  cat(
    "Minimum informative calls: ",
    MIN_CALLS,
    "\n",
    sep = ""
  )
  
  if (n_after == 0) {
    stop(
      "No windows remain after filtering for sample: ",
      sample_id
    )
  }
  
  # add sample identifier
  seg$sample_id <- sample_id
  
  # retain standardized columns only
  seg <- seg[,
    c(
      "seg_id",
      "seg_chrom",
      "seg_start",
      "seg_end",
      "seg_name",
      "seg_strand",
      "meth_calls",
      "unmeth_calls",
      "no_calls",
      "methfrac",
      "readcount",
      "motifcount",
      "total_calls",
      "pct_5mC",
      "sample_id"
      )
    ]


    # Return processed data
    seg
}
  

# == READ ALL SAMPLES ==
processed_samples <- lapply(
  names(SEGMETH_FILES),
  function(sample_id) {
    
    read_segmeth(
      file = SEGMETH_FILES[[sample_id]],
      sample_id = sample_id
    )
  }
)

# combine all retained windows
all_data <- do.call(
  rbind,
  processed_samples
)

rownames(all_data) <- NULL


cat(
  "\n============================================================\n",
  "COMBINED DATA\n",
  "============================================================\n",
  sep = ""
)

cat(
  "Total retained windows across all samples: ",
  nrow(all_data),
  "\n",
  sep = ""
)


# == CREATE BIOLOGICALLY MEANINGFUL PLOTTING GROUPS ==
all_data$plot_group <- NA_character_

# A. angusticeps
all_data$plot_group[
  all_data$sample_id == "angusticeps"
] <- "angusticeps"


# A. cozumelus: angusticeps-derived subgenome
all_data$plot_group[
  all_data$sample_id == "cozumelus" &
    all_data$seg_name == "ang_subgenome"
] <- "cozumelus_ang"


# A. cozumelus: deppii-derived subgenome
all_data$plot_group[
  all_data$sample_id == "cozumelus" &
    all_data$seg_name == "dep_subgenome"
] <- "cozumelus_dep"


# A. deppii Sample 1
all_data$plot_group[
  all_data$sample_id == "deppii_1"
] <- "deppii_1"


# A. deppii Sample 2
all_data$plot_group[
  all_data$sample_id == "deppii_2"
] <- "deppii_2"


# check that every retained window received a plotting group
if (any(is.na(all_data$plot_group))) {
  
  bad_rows <- all_data[
    is.na(all_data$plot_group),
    c("sample_id", "seg_name")
  ]
  
  stop(
    "Some retained windows could not be assigned to a plotting group.\n",
    "Unexpected sample_id / seg_name combinations:\n",
    paste(
      unique(
        paste(
          bad_rows$sample_id,
          bad_rows$seg_name,
          sep = " / "
        )
      ),
      collapse = "\n"
    )
  )
}

# preserve biologically meaningful ordering
PLOT_ORDER <- c(
  "angusticeps",
  "cozumelus_ang",
  "cozumelus_dep",
  "deppii_1",
  "deppii_2"
)

all_data$plot_group <- factor(
  all_data$plot_group,
  levels = PLOT_ORDER
)


# == CHECK WINDOW COUNTS ==
cat("\nRetained windows per plotting group:\n")

print(
  table(all_data$plot_group)
)


# == NUMERICAL SUMMARY ==
summary_list <- lapply(
  split(
    all_data,
    all_data$plot_group
  ),
  function(d) {
    
    data.frame(
      group =
        as.character(
          d$plot_group[1]
        ),
      
      n_windows =
        nrow(d),
      
      mean_5mC =
        mean(
          d$pct_5mC,
          na.rm = TRUE
        ),
      
      median_5mC =
        median(
          d$pct_5mC,
          na.rm = TRUE
        ),
      
      sd_5mC =
        sd(
          d$pct_5mC,
          na.rm = TRUE
        ),
      
      mean_informative_calls =
        mean(
          d$total_calls,
          na.rm = TRUE
        ),
      
      median_informative_calls =
        median(
          d$total_calls,
          na.rm = TRUE
        ),
      
      stringsAsFactors = FALSE
    )
  }
)


combined_summary <- do.call(
  rbind,
  summary_list
)

rownames(combined_summary) <- NULL


cat("\nCombined window-level summary:\n")

print(
  combined_summary,
  digits = 4
)


# == COMBINED STRIP / JITTER PLOT ==
# horizontal jitter separates overlapping windows visually
# jitter does not alter methylation percentage
# white diamond = median methylation for that group
p_strip_combined <- ggplot(
  all_data,
  aes(
    x = plot_group,
    y = pct_5mC,
    color = plot_group
  )
) +
  
  geom_jitter(
    width = 0.15,
    height = 0,
    alpha = POINT_ALPHA,
    size = POINT_SIZE
  ) +
  
  stat_summary(
    fun = median,
    geom = "point",
    shape = 23,
    size = 3,
    fill = "white",
    color = "black"
  ) +
  
  scale_color_manual(
    values = PLOT_COLORS
  ) +
  
  scale_x_discrete(
    labels = PLOT_LABELS
  ) +
  
  scale_y_continuous(
    breaks = seq(Y_MIN, Y_MAX, by = 20),
    expand = expansion(mult = c(0, 0.02))
  ) +
  
  coord_cartesian(
    ylim = c(Y_MIN, Y_MAX),
    clip = "off"
  ) +
  
  annotate(
    "segment",
    x = 2,
    xend = 3,
    y = -7,
    yend = -7,
    linewidth = 0.4
  ) +
  
  annotate(
    "segment",
    x = 4,
    xend = 5,
    y = -7,
    yend = -7,
    linewidth = 0.4
  ) +
  
  annotate(
    "text",
    x = 2.5,
    y = -13,
    label = "A. cozumelus",
    fontface = "italic",
    size = 3.5
  ) +
  
  annotate(
    "text",
    x = 4.5,
    y = -13,
    label = "A. deppii",
    fontface = "italic",
    size = 3.5
  ) +
  
  labs(
    x = NULL,
    y = "5mC (%)",
    title = "Genome-wide DNA methylation",
    subtitle = "5mC across 100-kb genomic windows"
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(
      size = 10,
      margin = margin(t = 5)
    ),
    
    plot.margin = margin(
      t = 5.5,
      r = 5.5,
      b = 30,
      l = 5.5
    )
  )


# == COMBINED VIOLIN PLOT ==
# violin width represents density of retained 100-kb windows at each methylation level
# white box = IQR
# black horizontal line inside the box = median window-level methylation
p_violin_combined <- ggplot(
  all_data,
  aes(
    x = plot_group,
    y = pct_5mC,
    fill = plot_group
  )
) +
  
  geom_violin(
    trim = FALSE,
    alpha = VIOLIN_ALPHA,
    color = "black",
    linewidth = 0.4
  ) +
  
  geom_boxplot(
    width = 0.12,
    fill = "white",
    color = "black",
    linewidth = 0.4,
    outlier.shape = NA
  ) +
  
  scale_fill_manual(
    values = PLOT_COLORS
  ) +
  
  scale_x_discrete(
    labels = PLOT_LABELS
  ) +
  
  scale_y_continuous(
    breaks = seq(Y_MIN, Y_MAX, by = 20),
    expand = expansion(mult = c(0, 0.02))
  ) +
  
  coord_cartesian(
    ylim = c(Y_MIN, Y_MAX),
    clip = "off"
  ) +
  
  annotate(
    "segment",
    x = 2,
    xend = 3,
    y = -7,
    yend = -7,
    linewidth = 0.4
  ) +
  
  annotate(
    "segment",
    x = 4,
    xend = 5,
    y = -7,
    yend = -7,
    linewidth = 0.4
  ) +
  
  annotate(
    "text",
    x = 2.5,
    y = -13,
    label = "A. cozumelus",
    fontface = "italic",
    size = 3.5
  ) +
  
  annotate(
    "text",
    x = 4.5,
    y = -13,
    label = "A. deppii",
    fontface = "italic",
    size = 3.5
  ) +
  
  labs(
    x = NULL,
    y = "5mC (%)",
    title = "Genome-wide DNA methylation",
    subtitle = "5mC across 100-kb genomic windows"
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(
        size = 10,
        margin = margin(t = 5)
    ),

    plot.margin = margin(
        t = 5.5,
        r = 5.5,
        b = 30,
        l = 5.5
    )
  )


# == SAVE OUTPUTS ==
dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# save numerical summary
write.csv(
  combined_summary,
  file.path(
    OUTPUT_DIR,
    "combined_segmeth_window_summary.csv"
  ),
  row.names = FALSE
)

# save combined plotting data
# (reproducibility, downstream statistical analyses)
write.csv(
  all_data[
    ,
    c(
      "sample_id",
      "plot_group",
      "seg_id",
      "seg_chrom",
      "seg_start",
      "seg_end",
      "seg_name",
      "meth_calls",
      "unmeth_calls",
      "total_calls",
      "methfrac",
      "pct_5mC"
    )
  ],
  file.path(
    OUTPUT_DIR,
    "combined_segmeth_plotting_data.csv"
  ),
  row.names = FALSE
)


# strip plot: PDF
ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "combined_segmeth_strip_plot.pdf"
  ),
  plot = p_strip_combined,
  width = 9,
  height = 5.5
)

# strip plot: PNG
ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "combined_segmeth_strip_plot.png"
  ),
  plot = p_strip_combined,
  width = 9,
  height = 5.5,
  dpi = 300
)

# violin plot: PDF
ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "combined_segmeth_violin_plot.pdf"
  ),
  plot = p_violin_combined,
  width = 9,
  height = 5.5
)

# violin plot: PNG
ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "combined_segmeth_violin_plot.png"
  ),
  plot = p_violin_combined,
  width = 9,
  height = 5.5,
  dpi = 300
)


# == FINISHED ==
cat(
  "\n============================================================\n",
  "Finished.\n",
  "Outputs written to:\n  ",
  OUTPUT_DIR,
  "\n",
  "============================================================\n",
  sep = ""
)
