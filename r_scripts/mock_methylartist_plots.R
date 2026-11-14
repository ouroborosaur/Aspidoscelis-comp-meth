# mock_methylartist_plots.R
# simulated conceptual outcomes for segmeth violin plots
# illustrate possible outcomes
# simulated mean 5mC within a 100-kb genomic window

# install package
library(ggplot2)

set.seed(20261107)

# == SETTINGS ==
OUTPUT_DIR <- "/scratch/maelstr1/thesis/02_processed/mock_segmeth_figures"

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

N_WINDOWS <- 5000

LINEAGE_COLORS <- c(
  "Maternal" = "#CC79A7",
  "Paternal" = "#56B4E9"
)

GROUP_LEVELS <- c(
  "A. angusticeps",
  "angusticeps-derived",
  "deppii-derived",
  "A. deppii"
)


# == SIMULATE DISTRIBUTION ==
# mean = approximate center of distribution
# sd = amount of window-to-window variation
# Values are restricted to 0-100%
simulate_windows <- function(
  n,
  mean,
  sd = 8
) {
  
  x <- rnorm(
    n,
    mean = mean,
    sd = sd
  )
  
  x <- pmax(
    x,
    0
  )
  
  x <- pmin(
    x,
    100
  )
  
  return(x)
}


# == CREATE GROUP ==
make_group <- function(
  group,
  lineage,
  mean,
  sd = 8
) {
  
  data.frame(
    group = group,
    lineage = lineage,
    
    methylation = simulate_windows(
      n = N_WINDOWS,
      mean = mean,
      sd = sd
    ),
    
    stringsAsFactors = FALSE
  )
}

# == PLOT ==
plot_mock_segmeth <- function(
    data,
    title_text,
    subtitle_text
) {
  
  data$group <- factor(
    data$group,
    levels = GROUP_LEVELS
  )
  
  data$lineage <- factor(
    data$lineage,
    levels = c(
      "Maternal",
      "Paternal"
    )
  )
  
  
  ggplot(
    data,
    aes(
      x = group,
      y = methylation,
      fill = lineage
    )
  ) +
    
    geom_violin(
      trim = FALSE,
      linewidth = 0.7,
      color = "black"
    ) +
    
    geom_boxplot(
      width = 0.13,
      fill = "white",
      color = "black",
      outlier.shape = NA,
      linewidth = 0.6
    ) +
    
    scale_fill_manual(
      values = LINEAGE_COLORS,
      labels = c(
        "Maternal",
        "Paternal"
      )
    ) +
    
    scale_x_discrete(
    labels = c(
        "A. angusticeps" =
            expression(italic("A. angusticeps")),

        "angusticeps-derived" =
            expression(
                italic("A. angusticeps") * "-derived"
            ),

        "deppii-derived" =
            expression(
                italic("A. deppii") * "-derived"
            ),

        "A. deppii" =
            expression(italic("A. deppii"))
    )
) +
    
    scale_y_continuous(
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
      x = NULL,
      y = "5mC (%)",
      fill = "Lineage",
      title = title_text,
      subtitle = subtitle_text,
      caption = "100-kb genomic windows (simulated)"
    ) +
    
    theme_classic(
      base_size = 14
    ) +
    
    theme(
      legend.position = "top",
      
      axis.text.x = element_text(
        angle = 20,
        hjust = 1
      ),
      
      plot.title = element_text(
        face = "bold",
        size = 17
      ),
      
      plot.subtitle = element_text(
        size = 13
      )
    )
}


# == SCENARIO 1 ==
# both subgenomes resemble parental genomes
# center all datasets at ~70%
scenario_parental <- rbind(
  
  make_group(
    "A. angusticeps",
    "Maternal",
    mean = 70
  ),
  
  make_group(
    "angusticeps-derived",
    "Maternal",
    mean = 70
  ),
  
  make_group(
    "deppii-derived",
    "Paternal",
    mean = 70
  ),
  
  make_group(
    "A. deppii",
    "Paternal",
    mean = 70
  )
)

p_parental <- plot_mock_segmeth(
  scenario_parental,
  title_text = "Parental-like subgenome methylation",
  subtitle_text = "Both hybrid subgenomes retain parental-like methylation"
)


# == SCENARIO 2 ==
# one subgenome changes more than the other
# angusticeps-derived stays near its parent (70% -> 70%)
# deppii-derived shifts strongly downward (70% -> 55%)
scenario_asymmetric <- rbind(
  
  make_group(
    "A. angusticeps",
    "Maternal",
    mean = 70
  ),
  
  make_group(
    "angusticeps-derived",
    "Maternal",
    mean = 70
  ),
  
  make_group(
    "deppii-derived",
    "Paternal",
    mean = 55
  ),
  
  make_group(
    "A. deppii",
    "Paternal",
    mean = 70
  )
)

p_asymmetric <- plot_mock_segmeth(
  scenario_asymmetric,
  title_text = "Ancestry-asymmetric methylation response",
  subtitle_text = "One subgenome changes while the other remains parental-like"
)


# == SCENARIO 3 ==
# both subgenomes diverge from parental genomes
# parental datasets 70%
# subgenomes 85%
scenario_divergent <- rbind(
  
  make_group(
    "A. angusticeps",
    "Maternal",
    mean = 70
  ),
  
  make_group(
    "angusticeps-derived",
    "Maternal",
    mean = 85
  ),
  
  make_group(
    "deppii-derived",
    "Paternal",
    mean = 85
  ),
  
  make_group(
    "A. deppii",
    "Paternal",
    mean = 70
  )
)

p_divergent <- plot_mock_segmeth(
  scenario_divergent,
  title_text = "Divergent subgenome methylation",
  subtitle_text = "Both hybrid subgenomes diverge from parental methylation"
)


# == SAVE PDFs ==
ggsave(
  file.path(
    OUTPUT_DIR,
    "mock_segmeth_subgenomes_resemble_parents.pdf"
  ),
  p_parental,
  width = 10,
  height = 7
)

ggsave(
  file.path(
    OUTPUT_DIR,
    "mock_segmeth_one_subgenome_changes_more.pdf"
  ),
  p_asymmetric,
  width = 10,
  height = 7
)

ggsave(
  file.path(
    OUTPUT_DIR,
    "mock_segmeth_both_subgenomes_diverge.pdf"
  ),
  p_divergent,
  width = 10,
  height = 7
)

# == SAVE PNGs ==
ggsave(
  file.path(
    OUTPUT_DIR,
    "mock_segmeth_subgenomes_resemble_parents.png"
  ),
  p_parental,
  width = 10,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(
    OUTPUT_DIR,
    "mock_segmeth_one_subgenome_changes_more.png"
  ),
  p_asymmetric,
  width = 10,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(
    OUTPUT_DIR,
    "mock_segmeth_both_subgenomes_diverge.png"
  ),
  p_divergent,
  width = 10,
  height = 7,
  dpi = 300
)


# == VERIFY SIMULATION PRODUCED DIFFERENT OUTCOMES ==
print_scenario_means <- function(
    data,
    scenario_name
) {
  
  cat(
    "\n",
    scenario_name,
    "\n",
    sep = ""
  )
  
  x <- aggregate(
    methylation ~ group,
    data = data,
    FUN = mean
  )
  
  x$methylation <- round(
    x$methylation,
    1
  )
  
  print(x)
}


print_scenario_means(
  scenario_parental,
  "PARENTAL-LIKE"
)

print_scenario_means(
  scenario_asymmetric,
  "ASYMMETRIC"
)

print_scenario_means(
  scenario_divergent,
  "BOTH DIVERGE"
)