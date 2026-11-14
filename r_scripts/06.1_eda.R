# 06.1_eda.R
# Coordinate-independent exploratory methylation analysis

# OUTPUTS
#  per_sample_methylation_summary.csv
#  per_chromosome_methylation_summary.csv
#  figures/eda/methylation_distributions.pdf
#  figures/eda/coverage_distributions.pdf
#  figures/eda/chromosome_mean_methylation.pdf
#  rds_objects/coordinate_independent_eda.rds


# install package
library(optparse)

# parse command-line arguments
option_list <- list(
  # root containing sample directories
  make_option("--processed_dir", type = "character"),

  # root to write analysis results
  make_option("--output_dir", type = "character"),

  # FALSE = use standard final BEDs; TRUE = use repeat-filtered BEDs
  make_option("--repeat_filtered", type = "character", default = "FALSE")
)
opt <- parse_args(OptionParser(option_list = option_list))

# convert text argument to a logical flag
repeat_filtered <- toupper(opt$repeat_filtered) == "TRUE"
RUN_LABEL <- if (repeat_filtered) "repeat_filtered" else "all_retained_sites"

# ANALYSIS DESIGN
# all retained sites (repeat_filtered = FALSE)
#  angusticeps parent
#  deppii ANMO
#  deppii LNG
#  cozumelus angusticeps-derived subgenome
#  cozumelus deppii-derived subgenome

# repeat filtered (repeat_filtered = TRUE)
#  deppii ANMO
#  deppii LNG
#  cozumelus angusticeps-derived subgenome
#  cozumelus deppii-derived subgenome
#  (no repeat-filtered angusticeps parental BED currently available)

# explicit asymmetry avoids silent substitution of unfiltered file when 
#  repeat-filtered file is unavailable
if (repeat_filtered) {
  PARENT_SUFFIX <- "_final_norepeat.bed"
  COZU_ANG_BED <- "subgenome_angusticeps/cozumelus_ang_subgenome_norepeat.bed"
  COZU_DEP_BED <- "subgenome_deppii/cozumelus_dep_subgenome_norepeat.bed"
  EXPECTED_DATASETS <- c(
    "deppii_ANMO5280",
    "deppii_LNG1403",
    "cozumelus_ang",
    "cozumelus_dep"
  )
    OMITTED_DATASETS <- c("angusticeps_ANMO5228")
} else {
  PARENT_SUFFIX <- "_final.bed"
  COZU_ANG_BED <- "subgenome_angusticeps/cozumelus_ang_subgenome.bed"
  COZU_DEP_BED <- "subgenome_deppii/cozumelus_dep_subgenome.bed"
  EXPECTED_DATASETS <- c(
    "angusticeps_ANMO5228",
    "deppii_ANMO5280",
    "deppii_LNG1403",
    "cozumelus_ang",
    "cozumelus_dep"
  )
  OMITTED_DATASETS <- character(0)
}

# keep all-sites and repeat-filtered results in separate directories
OUTPUT_DIR <- file.path(opt$output_dir, RUN_LABEL)
RDS_DIR <- file.path(OUTPUT_DIR, "rds_objects")
EDA_DIR <- file.path(OUTPUT_DIR, "figures", "eda")
dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(EDA_DIR, recursive = TRUE, showWarnings = FALSE)

cat("COORDINATE-INDEPENDENT EXPLORATORY DATA ANALYSIS\n")
cat("Run type:", RUN_LABEL, "\n\n")

# columns that must be present in every plot-ready BED
required_cols <- c("chrom", "start", "end", "methyl", "strand", "prob", "cov")


# == BED LOADING & VALIDATING FUNCTION ==
# read one BED table
# validate the schema and numeric ranges
# add analysis-friendly metadata/column names
load_bed <- function(filepath, sample_id, ancestry) {

  # fail early with a useful path if a filename in the pipeline is wrong
  if (!file.exists(filepath)) {
    stop(paste("ERROR: File not found:", filepath))
  }

  # check.names=FALSE preserves input header 
  # comment.char="" prevents R from treating any text as a comment
  bed <- read.table(
    filepath,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = ""
  )

  if (nrow(bed) == 0) {
    stop(paste("ERROR: File is empty:", filepath))
  }

  # confirm upstream filtering script produced the expected schema
  missing_cols <- setdiff(required_cols, colnames(bed))
  if (length(missing_cols) > 0) {
    stop(paste(
      "ERROR:", basename(filepath), "is missing columns:",
      paste(missing_cols, collapse = ", ")
    ))
  }

  # ignore unexpected extra columns in Phase 1 and retain known schema
  bed <- bed[, required_cols]

  # explicit conversion prevents character data from silently propagating into
  #  summary statistics
  # start & end are retained for provenance and Phase 2
  bed$start <- as.numeric(bed$start)
  bed$end <- as.numeric(bed$end)
  bed$prob <- as.numeric(bed$prob)
  bed$cov <- as.numeric(bed$cov)

  # validate methylation percentage and coverage before plotting/statistics
  if (any(!is.finite(bed$prob)) || any(!is.finite(bed$cov))) {
    stop(paste(
      "ERROR: Non-numeric methylation percentage or coverage in",
      filepath
    ))
  }
  if (any(bed$prob < 0 | bed$prob > 100)) {
    stop(paste("ERROR: Methylation percentage outside 0-100 in", filepath))
  }
  if (any(bed$cov < 0)) {
    stop(paste("ERROR: Negative coverage in", filepath))
  }

  # plot-ready files are intended to contain 5mC rows only
  # if another modkit modification code survived upstream filtering, warn and remove it
  if (any(bed$methyl != "m")) {
    warning(paste(
      "Non-5mC rows detected in", basename(filepath),
      "- retaining only methyl == 'm'."
    ))
      bed <- bed[bed$methyl == "m", , drop = FALSE]
  }

  # harmonize CHROMOSOME LABELS only
  #   angusticeps_Chr_1 -> Chr_1
  #   deppii_Chr_1      -> Chr_1
  #   Chr_1             -> Chr_1
  # enables chromosome-level matching in 06.2_stats.R
  # does NOT imply numerical base-pair coordinates within those chromosomes are equal
  bed$chrom_parent_label <- sub("^(angusticeps_|deppii_)", "", bed$chrom)

  # attach dataset identity to every row so tables can later be combined safely
  bed$sample <- sample_id
  bed$ancestry <- ancestry

  # rename two historically named columns for clarity
  # (no numerical transformation)
  bed$pct_5mC <- bed$prob
  bed$valid_cov <- bed$cov

  cat("  Loaded:", sample_id, "-", nrow(bed), "5mC sites\n")
  return(bed)
}

# construct input paths valid for this run type
paths <- list(
  deppii_ANMO5280 = file.path(
    opt$processed_dir, "deppii_ANMO5280", "plot_ready",
    paste0("deppii_ANMO5280", PARENT_SUFFIX)
  ),
  deppii_LNG1403 = file.path(
    opt$processed_dir, "deppii_LNG1403", "plot_ready",
    paste0("deppii_LNG1403", PARENT_SUFFIX)
  ),
  cozumelus_ang = file.path(
    opt$processed_dir, "cozumelus_LNG1401", "plot_ready", COZU_ANG_BED
  ),
  cozumelus_dep = file.path(
    opt$processed_dir, "cozumelus_LNG1401", "plot_ready", COZU_DEP_BED
  )
)

# angusticeps is added only to all-retained-sites run
if (!repeat_filtered) {
  paths <- c(
    list(
      angusticeps_ANMO5228 = file.path(
        opt$processed_dir,
        "angusticeps_ANMO5228",
        "plot_ready",
        "angusticeps_ANMO5228_final.bed"
      )
    ),
    paths
  )
}

# confirm the path list matches intended design exactly
# safeguard so if a future edit accidentally adds/removes a dataset,
#  the script stops rather than quietly changing the analysis.
if (!setequal(names(paths), EXPECTED_DATASETS)) {
  stop("ERROR: Internal dataset list does not match the intended run design.")
}

# load available datasets
# ancestry labels are descriptive metadata
# distinguish sexual parents from ancestry-assigned hybrid subgenomes without 
#  changing genomic coordinates
beds <- list()

if (!repeat_filtered) {
  beds$angusticeps_ANMO5228 <- load_bed(
    paths$angusticeps_ANMO5228,
    "angusticeps_ANMO5228",
    "angusticeps_parent"
  )
}

beds$deppii_ANMO5280 <- load_bed(
  paths$deppii_ANMO5280,
  "deppii_ANMO5280",
  "deppii_parent"
)
beds$deppii_LNG1403 <- load_bed(
  paths$deppii_LNG1403,
  "deppii_LNG1403",
  "deppii_parent"
)
beds$cozumelus_ang <- load_bed(
  paths$cozumelus_ang,
  "cozumelus_ang",
  "angusticeps_derived_subgenome"
)
beds$cozumelus_dep <- load_bed(
  paths$cozumelus_dep,
  "cozumelus_dep",
  "deppii_derived_subgenome"
)

cat("Datasets included in this EDA run:", paste(names(beds), collapse = ", "), "\n")
if (length(OMITTED_DATASETS) > 0) {
  cat(
    "Datasets intentionally omitted:",
    paste(OMITTED_DATASETS, collapse = ", "),
    "\nReason: no matched repeat-filtered file is available.\n\n"
  )
}

# per-dataset descriptive methylation/coverage summaries
# statistics are coordinate-independent
sample_summary <- do.call(rbind, lapply(beds, function(d) {
  data.frame(
    sample = d$sample[1],
    ancestry = d$ancestry[1],
    total_sites = nrow(d),

    # central tendency and dispersion of % 5mC across retained CpGs
    mean_pct_5mC = mean(d$pct_5mC, na.rm = TRUE),
    median_pct_5mC = median(d$pct_5mC, na.rm = TRUE),
    sd_pct_5mC = sd(d$pct_5mC, na.rm = TRUE),

    # simple descriptive fractions of highly/weakly methylated sites
    pct_sites_high_80 = mean(d$pct_5mC >= 80, na.rm = TRUE) * 100,
    pct_sites_low_20 = mean(d$pct_5mC <= 20, na.rm = TRUE) * 100,

    # coverage summaries help identify systematic sequencing-depth differences
    mean_valid_cov = mean(d$valid_cov, na.rm = TRUE),
    median_valid_cov = median(d$valid_cov, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
rownames(sample_summary) <- NULL

# per-chromosome descriptive summaries
# each dataset is split by its harmonized ancestry-aware chromosome label
# output gives one row per sample x chromosome
# rows become the input to chromosome-level parent/hybrid matching in 06.2_stats.R
chrom_summary <- do.call(rbind, lapply(beds, function(d) {
  split_d <- split(d, d$chrom_parent_label)

  do.call(rbind, lapply(split_d, function(x) {
    data.frame(
      sample = x$sample[1],
      ancestry = x$ancestry[1],

      # harmonized label used only for chromosome-level correspondence
      chrom = x$chrom_parent_label[1],

      # original assembly-specific label retained for provenance
      original_chrom = x$chrom[1],

      n_sites = nrow(x),
      mean_pct_5mC = mean(x$pct_5mC, na.rm = TRUE),
      median_pct_5mC = median(x$pct_5mC, na.rm = TRUE),
      sd_pct_5mC = sd(x$pct_5mC, na.rm = TRUE),
      mean_valid_cov = mean(x$valid_cov, na.rm = TRUE),
      median_valid_cov = median(x$valid_cov, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
}))
rownames(chrom_summary) <- NULL

# write human-readable CSV outputs before figure generation
# inspect before running downstream statistics
write.csv(
  sample_summary,
  file.path(OUTPUT_DIR, "per_sample_methylation_summary.csv"),
  row.names = FALSE
)
write.csv(
  chrom_summary,
  file.path(OUTPUT_DIR, "per_chromosome_methylation_summary.csv"),
  row.names = FALSE
)

# == EDA FIGURES ==
# methylation distributions
# each histogram shows distribution of % 5mC across retained CpGs in one dataset
pdf(file.path(EDA_DIR, "methylation_distributions.pdf"), width = 11, height = 8.5)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
for (nm in names(beds)) {
  d <- beds[[nm]]
  hist(
    d$pct_5mC,
    breaks = seq(0, 100, by = 5),
    main = nm,
    xlab = "5mC (%)",
    ylab = "CpG sites"
  )
}
par(mfrow = c(1, 1))
dev.off()

# valid-coverage distributions
# restrict each displayed histogram to <= the 99.5th percentile
# prevent a extreme sites from compressing the useful portion of the plot
# affects visualization only; no rows are removed from the stored data or 
#  summary tables
pdf(file.path(EDA_DIR, "coverage_distributions.pdf"), width = 11, height = 8.5)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
for (nm in names(beds)) {
  d <- beds[[nm]]
  upper <- as.numeric(quantile(d$valid_cov, 0.995, na.rm = TRUE))
  hist(
    d$valid_cov[d$valid_cov <= upper],
    breaks = 50,
    main = nm,
    xlab = "Valid coverage (<=99.5th percentile)",
    ylab = "CpG sites"
  )
}
par(mfrow = c(1, 1))
dev.off()

# mean methylation by chromosome within each dataset
# identify broad chromosome-scale differences and possible outlier chromosomes
pdf(file.path(EDA_DIR, "chromosome_mean_methylation.pdf"), width = 12, height = 8)
par(mfrow = c(3, 2), mar = c(7, 4, 3, 1))
for (nm in names(beds)) {
  x <- chrom_summary[chrom_summary$sample == nm, ]
  x <- x[order(x$chrom), ]
  barplot(
    x$mean_pct_5mC,
    names.arg = x$chrom,
    las = 2,
    cex.names = 0.65,
    ylim = c(0, 100),
    ylab = "Mean 5mC (%)",
    main = nm
  )
}
par(mfrow = c(1, 1))
dev.off()

# == SAVE COMPLETE R OBJECT ==
# RDS preserves validated BED tables and summaries
# avoid repeated disk parsing
# ensure stats runs use the same data as EDA figure runs
analysis_data <- list(
  run_label = RUN_LABEL,
  repeat_filtered = repeat_filtered,
  included_datasets = names(beds),
  omitted_datasets = OMITTED_DATASETS,
  paths = paths,
  beds = beds,
  sample_summary = sample_summary,
  chrom_summary = chrom_summary,
  coordinate_note = paste(
    "Chromosome labels are harmonized for chromosome-level summaries only.",
    "Genomic coordinates remain assembly-specific and are not matched across references.",
    if (repeat_filtered) {
      "Repeat-filtered angusticeps is unavailable; cozumelus_ang is descriptive only and is not formally compared with angusticeps in Phase 1."
    } else {
      "All five datasets are available, so both ancestry-matched chromosome-level comparisons can be performed."
    }
  )
)

saveRDS(
  analysis_data,
  file.path(RDS_DIR, "coordinate_independent_eda.rds")
)

cat("\nEDA complete.\n")
cat("Results:", OUTPUT_DIR, "\n")
cat("IMPORTANT: No CpG coordinates were matched across reference assemblies.\n")
