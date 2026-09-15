options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]]),
                        winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
repo_dir <- normalizePath(file.path(script_dir, ".."),
                          winslash = "/", mustWork = TRUE)
source_data_file <- Sys.getenv(
  "KIRC_BPS_CIBERSORT_SOURCE_DATA",
  unset = file.path(
    repo_dir, "data", "derived", "source_data",
    "30_ESRRB_CIBERSORT_correlation_source_data_187_tumors.csv"
  )
)
output_root <- Sys.getenv(
  "KIRC_BPS_FIGURE_OUTPUT_DIR",
  unset = file.path(repo_dir, "results", "reproduced_ESRRB_CIBERSORT")
)
figure_dir <- file.path(output_root, "figures")
table_dir <- file.path(output_root, "tables")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(source_data_file)) {
  stop("Missing public source-data file: ", source_data_file)
}
source_data <- utils::read.csv(
  source_data_file, check.names = FALSE, stringsAsFactors = FALSE
)
required_columns <- c("profile_id", "ESRRB_log2_CPM")
if (!all(required_columns %in% colnames(source_data))) {
  stop("Source-data file lacks required columns: ",
       paste(setdiff(required_columns, colnames(source_data)), collapse = ", "))
}

immune_cells <- setdiff(
  colnames(source_data),
  c("profile_id", "ESRRB_log2_CPM")
)
if (length(immune_cells) != 22L) {
  stop("Expected 22 CIBERSORT immune-cell columns; found ", length(immune_cells), ".")
}

valid_tumor_ids <- source_data$profile_id
if (length(valid_tumor_ids) < 3L) stop("Fewer than three valid tumor profiles.")

correlation_rows <- lapply(immune_cells, function(cell) {
  x <- source_data$ESRRB_log2_CPM
  y <- source_data[[cell]]
  complete <- is.finite(x) & is.finite(y)
  n_complete <- sum(complete)
  zero_variance <- n_complete < 3L || stats::sd(y[complete]) == 0

  if (zero_variance) {
    return(data.frame(
      gene = "ESRRB", cell_type = cell, n = n_complete,
      rho = NA_real_, p_value = NA_real_, estimable = FALSE
    ))
  }

  test <- suppressWarnings(stats::cor.test(
    x[complete], y[complete], method = "spearman", exact = FALSE
  ))
  data.frame(
    gene = "ESRRB", cell_type = cell, n = n_complete,
    rho = unname(test$estimate), p_value = test$p.value,
    estimable = TRUE
  )
})
correlation_stats <- do.call(rbind, correlation_rows)
correlation_stats$BH_adjusted_p <- stats::p.adjust(correlation_stats$p_value, method = "BH")
correlation_stats$significant_BH_0.05 <- with(
  correlation_stats,
  !is.na(BH_adjusted_p) & BH_adjusted_p < 0.05
)

write.csv(
  correlation_stats,
  file.path(table_dir, "29_ESRRB_CIBERSORT_correlation_statistics_for_figure.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  source_data,
  file.path(
    table_dir,
    sprintf(
      "30_ESRRB_CIBERSORT_correlation_source_data_%d_tumors.csv",
      length(valid_tumor_ids)
    )
  ),
  row.names = FALSE,
  na = ""
)

format_p <- function(p) {
  ifelse(
    is.na(p), "Not estimable",
    ifelse(p < 1e-4, formatC(p, format = "e", digits = 1),
           ifelse(p < 0.01, sprintf("%.4f", p), sprintf("%.3f", p)))
  )
}

palette <- c(
  positive = "#347FA5",
  negative = "#D1783B",
  nonsignificant = "#A7A7A7",
  not_estimable = "#575757"
)

theme_manuscript <- function(base_size = 7.2) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.2, colour = "black"),
      plot.title = element_text(size = base_size + 0.8, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.1, colour = "#404040", hjust = 0),
      plot.caption = element_text(size = base_size - 0.7, colour = "#555555", hjust = 0),
      plot.tag = element_text(size = base_size + 1.2, face = "bold"),
      panel.grid = element_blank(),
      legend.position = "none",
      plot.margin = margin(5, 6, 5, 5)
    )
}

overview_data <- correlation_stats
overview_data$display_group <- with(
  overview_data,
  ifelse(
    !estimable, "not_estimable",
    ifelse(
      significant_BH_0.05 & rho > 0, "positive",
      ifelse(significant_BH_0.05 & rho < 0, "negative", "nonsignificant")
    )
  )
)
overview_data$q_label <- paste0("BH P = ", format_p(overview_data$BH_adjusted_p))
overview_data$q_label[!overview_data$estimable] <- "Zero variance"
overview_data$cell_type <- factor(
  overview_data$cell_type,
  levels = overview_data$cell_type[order(overview_data$rho, na.last = TRUE)]
)

p_overview <- ggplot(overview_data, aes(y = cell_type)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "#555555") +
  geom_segment(
    data = subset(overview_data, estimable),
    aes(x = 0, xend = rho, yend = cell_type),
    linewidth = 0.45,
    colour = "#C7C7C7"
  ) +
  geom_point(
    data = subset(overview_data, estimable),
    aes(x = rho, colour = display_group, size = significant_BH_0.05),
    alpha = 0.95
  ) +
  geom_point(
    data = subset(overview_data, !estimable),
    aes(x = 0),
    shape = 4,
    stroke = 0.65,
    size = 2.1,
    colour = palette[["not_estimable"]]
  ) +
  geom_text(
    aes(x = 0.335, label = q_label),
    hjust = 0,
    size = 2.25,
    colour = "#303030"
  ) +
  scale_colour_manual(values = palette[c("positive", "negative", "nonsignificant")]) +
  scale_size_manual(values = c(`FALSE` = 1.65, `TRUE` = 2.45)) +
  scale_x_continuous(
    limits = c(-0.32, 0.57),
    breaks = c(-0.3, -0.15, 0, 0.15, 0.3),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    title = "ESRRB expression and inferred immune-cell fractions",
    subtitle = paste0(
      "Spearman correlations in ", length(valid_tumor_ids),
      " TCGA-KIRC tumors; BH correction across ", sum(correlation_stats$estimable),
      " estimable tests"
    ),
    x = "Spearman correlation coefficient (rho)",
    y = NULL
  ) +
  theme_manuscript() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 6.6),
    plot.margin = margin(5, 5, 5, 7)
  )

make_scatter <- function(cell, line_colour) {
  stat_row <- correlation_stats[correlation_stats$cell_type == cell, , drop = FALSE]
  plot_data <- data.frame(
    ESRRB_log2_CPM = source_data$ESRRB_log2_CPM,
    immune_fraction_percent = 100 * source_data[[cell]]
  )
  label <- paste0(
    "n = ", stat_row$n,
    "\nSpearman rho = ", sprintf("%.3f", stat_row$rho),
    "\nBH-adjusted P = ", format_p(stat_row$BH_adjusted_p)
  )
  y_upper <- max(plot_data$immune_fraction_percent, na.rm = TRUE) * 1.06

  ggplot(plot_data, aes(x = ESRRB_log2_CPM, y = immune_fraction_percent)) +
    geom_point(size = 1.05, alpha = 0.46, colour = "#4C4C4C") +
    geom_smooth(
      method = "loess", formula = y ~ x, se = TRUE,
      span = 0.85, linewidth = 0.7,
      colour = line_colour, fill = line_colour, alpha = 0.16
    ) +
    annotate(
      "label",
      x = -Inf, y = Inf, label = label,
      hjust = -0.04, vjust = 1.12,
      size = 2.35, lineheight = 1.05,
      linewidth = 0.25, label.padding = grid::unit(0.12, "lines"),
      colour = "#202020", fill = scales::alpha("white", 0.88)
    ) +
    labs(
      title = cell,
      x = "ESRRB expression (log2 CPM)",
      y = "Estimated cell fraction (%)"
    ) +
    coord_cartesian(ylim = c(0, y_upper), clip = "on") +
    theme_manuscript()
}

p_nk <- make_scatter("NK cells resting", palette[["positive"]])
p_neutrophils <- make_scatter("Neutrophils", palette[["negative"]])

combined_figure <- p_overview / (p_nk | p_neutrophils) +
  plot_layout(heights = c(2.25, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 8.4, face = "bold"))

scatter_figure <- p_nk | p_neutrophils

save_publication_figure <- function(plot, stem, width_mm, height_mm, dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  svglite::svglite(paste0(stem, ".svg"), width = width_in, height = height_in)
  print(plot)
  grDevices::dev.off()

  grDevices::cairo_pdf(
    paste0(stem, ".pdf"), width = width_in, height = height_in,
    family = "Arial", onefile = TRUE
  )
  print(plot)
  grDevices::dev.off()

  ragg::agg_tiff(
    paste0(stem, ".tiff"), width = width_in, height = height_in,
    units = "in", res = dpi, compression = "lzw", background = "white"
  )
  print(plot)
  grDevices::dev.off()

  ragg::agg_png(
    paste0(stem, ".png"), width = width_in, height = height_in,
    units = "in", res = 300, background = "white"
  )
  print(plot)
  grDevices::dev.off()
}

save_publication_figure(
  p_overview,
  file.path(figure_dir, "Figure_ESRRB_22_immune_cell_correlation_overview"),
  width_mm = 165, height_mm = 150
)
save_publication_figure(
  scatter_figure,
  file.path(figure_dir, "Figure_ESRRB_significant_immune_correlations_scatter"),
  width_mm = 175, height_mm = 78
)
save_publication_figure(
  combined_figure,
  file.path(figure_dir, "Figure_ESRRB_immune_correlations_combined"),
  width_mm = 183, height_mm = 190
)

legend_text <- paste(
  "Fig. X | Association between ESRRB expression and CIBERSORT-inferred immune-cell fractions in TCGA-KIRC tumors.",
  "a, Spearman correlation coefficients between ESRRB expression and the relative fractions of 22 immune-cell types",
  paste0("among ", length(valid_tumor_ids), " tumors with CIBERSORT deconvolution P < 0.05. "),
  paste0("P values from the ", sum(correlation_stats$estimable), " estimable correlations were adjusted using the Benjamini-Hochberg method. "),
  "The gamma-delta T-cell fraction had zero variance and was therefore not estimable.",
  "BH-significant positive and negative associations are highlighted in blue and orange, respectively.",
  "b,c, Sample-level associations of ESRRB expression with resting NK-cell and neutrophil fractions.",
  "Curves show locally estimated trends with 95% confidence bands and are included for visualization only.",
  "rho, Spearman's rank correlation coefficient; BH, Benjamini-Hochberg.",
  "Source data are provided in the accompanying source-data file.",
  sep = " "
)
writeLines(
  legend_text,
  file.path(figure_dir, "Figure_ESRRB_immune_correlations_legend.txt"),
  useBytes = TRUE
)

qa_notes <- c(
  "Figure contract",
  "Core conclusion: ESRRB expression is associated positively with the resting NK-cell fraction and negatively with the neutrophil fraction in TCGA-KIRC tumors, without implying causality.",
  "Archetype: Quantitative grid with a correlation-overview hero panel and two sample-level supporting panels.",
  "Target/output: Scientific Reports revision; editable SVG/PDF plus 600-dpi TIFF and 300-dpi PNG preview.",
  "Backend: R (ggplot2 and patchwork) only.",
  "Final combined size: 183 x 190 mm.",
  paste0("Replicate unit: individual TCGA-KIRC tumor profile; n = ", length(valid_tumor_ids), "."),
  "Panel a: Spearman rho for 22 CIBERSORT cell fractions; 21 estimable tests; Benjamini-Hochberg adjustment across estimable P values.",
  "Panels b,c: raw sample-level values; LOESS curve and 95% confidence band are descriptive only and are not the inferential test.",
  "Exclusion rule: only tumor profiles with CIBERSORT deconvolution P < 0.05 were included; no additional observations were removed.",
  "Not estimable: gamma-delta T-cell fraction had zero variance.",
  "Reviewer risk: CIBERSORT fractions are compositional estimates from the same bulk RNA-seq profiles as ESRRB expression; associations do not establish immune recruitment, direct regulation, or causality.",
  "Visual QA: all panels inspected at final size; y axes for immune fractions are constrained to non-negative values; labels and confidence bands do not collide with panel boundaries."
)
writeLines(
  qa_notes,
  file.path(figure_dir, "Figure_ESRRB_immune_correlations_QA_notes.txt"),
  useBytes = TRUE
)

cat("Valid tumor profiles:", length(valid_tumor_ids), "\n")
cat("Immune-cell types assessed:", length(immune_cells), "\n")
cat("Estimable correlations:", sum(correlation_stats$estimable), "\n")
cat("BH-significant correlations:", sum(correlation_stats$significant_BH_0.05), "\n")
print(correlation_stats[correlation_stats$significant_BH_0.05, ])
cat("Figures and source-data tables saved under:", output_root, "\n")
