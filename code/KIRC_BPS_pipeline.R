# ==============================================================================
# Corrected TCGA-KIRC / BPS computational reanalysis
# ==============================================================================
# Purpose
#   1. Reconstruct the BPS candidate-target universe from raw database exports
#      with source-specific filters and official human gene-symbol mapping.
#   2. Rebuild TCGA-KIRC STAR-Counts locally with reproducible sample metadata.
#   3. Apply low-expression filtering, TMM normalization, voom-limma modelling,
#      and patient-level blocking for partially paired/repeated observations.
#   4. Assess residual technical variation without blindly removing it.
#   5. Recalculate DEGs, BPS-ccRCC overlap, enrichment, nested-CV machine
#      learning performance, paired GSE53757 validation, optional paired local
#      RNA-seq validation, and CIBERSORT sensitivity analyses.
#
# Interpretation boundary
#   All targets and mechanisms produced here are computational candidates or
#   hypotheses. The workflow does not establish BPS exposure, direct binding,
#   target engagement, causality, or an experimentally validated AOP.
# ==============================================================================

options(stringsAsFactors = FALSE, timeout = 1200, warn = 1)
RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion",
        sample.kind = "Rejection")
set.seed(20260901)
pipeline_version <- "2026-09-14-github-release-candidate-v1"

# ---- Portable repository paths -----------------------------------------------
# Run this file from any working directory with:
#   Rscript code/KIRC_BPS_pipeline.R
# Optional environment variables can override large or restricted inputs:
#   KIRC_BPS_TCGA_DIR, KIRC_BPS_INDEPENDENT_COUNT_FILE,
#   KIRC_BPS_KEGG_SNAPSHOT, and KIRC_BPS_OUTPUT_DIR.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
SCRIPT_DIR <- if (length(script_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]]),
                        winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
REPO_DIR <- normalizePath(file.path(SCRIPT_DIR, ".."),
                          winslash = "/", mustWork = TRUE)
DATA_DIR <- file.path(REPO_DIR, "data")

TCGA_DIR <- Sys.getenv(
  "KIRC_BPS_TCGA_DIR",
  unset = file.path(
    DATA_DIR, "external", "TCGA-KIRC", "Transcriptome_Profiling",
    "Gene_Expression_Quantification"
  )
)
BPS_TARGET_FILE <- file.path(DATA_DIR, "input", "legacy_BPS_210.txt")
TARGET_RAW_DIR <- file.path(DATA_DIR, "input", "target_database_exports")
CHEMBL_TARGET_FILE <- file.path(TARGET_RAW_DIR, "ChEMBL_target.txt")
SEA_TARGET_FILE <- file.path(TARGET_RAW_DIR, "SEA_results.csv")
STITCH_TARGET_FILE <- file.path(TARGET_RAW_DIR, "STITCH_interactions.tsv")
SUPERPRED_TARGET_FILE <- file.path(TARGET_RAW_DIR, "SuperPred_targets.csv")
SWISS_TARGET_FILE <- file.path(TARGET_RAW_DIR, "SwissTargetPrediction.csv")
OLD_DEG_FILE <- file.path(DATA_DIR, "input", "legacy_KIRC_DEG.txt")
GEO_MATRIX_FILE <- file.path(
  DATA_DIR, "external", "GSE53757", "GSE53757_series_matrix.txt.gz"
)
INDEPENDENT_COUNT_FILE <- Sys.getenv(
  "KIRC_BPS_INDEPENDENT_COUNT_FILE",
  unset = file.path(DATA_DIR, "restricted", "independent_count.xlsx")
)

OUTPUT_DIR <- Sys.getenv(
  "KIRC_BPS_OUTPUT_DIR",
  unset = file.path(REPO_DIR, "results", "reproduced")
)
TABLE_DIR <- file.path(OUTPUT_DIR, "tables")
FIGURE_DIR <- file.path(OUTPUT_DIR, "figures")
CACHE_DIR <- file.path(OUTPUT_DIR, "cache")
LOG_DIR <- file.path(OUTPUT_DIR, "logs")
RDATA_FILE <- file.path(CACHE_DIR, "KIRC_BPS_codex.RData")
EXPRESSION_CHECKPOINT_FILE <- file.path(
  CACHE_DIR, "KIRC_BPS_expression_checkpoint.RData"
)
CORE_CHECKPOINT_FILE <- file.path(CACHE_DIR, "KIRC_BPS_core_checkpoint.RData")
CIBERSORT_CHECKPOINT_FILE <- file.path(
  CACHE_DIR, "CIBERSORT_LM22_1000perm_checkpoint.RData"
)
KEGG_SNAPSHOT_FILE <- Sys.getenv(
  "KIRC_BPS_KEGG_SNAPSHOT",
  unset = file.path(DATA_DIR, "restricted", "KEGG_hsa_snapshot_20260906.rds")
)
KEGG_SNAPSHOT_SHA256 <- "ea9aa2198195127197fb4db093da2efc0730793472f47b9aa654b212967b475d"

for (d in c(OUTPUT_DIR, TABLE_DIR, FIGURE_DIR, CACHE_DIR, LOG_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

LOG_FILE <- file.path(LOG_DIR, "KIRC_BPS_codex_analysis.log")
log_con <- file(LOG_FILE, open = "wt", encoding = "UTF-8")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message", append = TRUE)
on.exit({
  while (sink.number(type = "message") > 0) sink(type = "message")
  while (sink.number(type = "output") > 0) sink(type = "output")
  try(close(log_con), silent = TRUE)
}, add = TRUE)

message_stamp <- function(...) {
  txt <- paste0(..., collapse = "")
  message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), txt))
}

message_stamp("Analysis started")

# ---- Packages ----------------------------------------------------------------
required_packages <- c(
  "data.table", "httr", "jsonlite", "edgeR", "limma", "AnnotationDbi",
  "org.Hs.eg.db", "ggplot2", "ggrepel", "patchwork", "pheatmap",
  "glmnet", "e1071", "pROC", "clusterProfiler", "enrichplot",
  "svglite", "ragg", "GEOquery", "Biobase", "hgu133plus2.db",
  "readxl", "CIBERSORT", "digest"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(patchwork)
  library(glmnet)
  library(e1071)
  library(pROC)
})

package_versions <- data.frame(
  package = required_packages,
  version = vapply(
    required_packages,
    function(x) as.character(utils::packageVersion(x)),
    character(1)
  )
)
utils::write.csv(
  package_versions,
  file.path(TABLE_DIR, "00_package_versions.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# Lightweight cache signature. The independently generated SHA-256 manifest
# remains the archival proof of input identity; this signature prevents stale
# checkpoints from being reused after input metadata, package versions, or the
# declared pipeline version change.
signature_tcga_files <- list.files(
  TCGA_DIR,
  pattern = "rna_seq\\.augmented_star_gene_counts\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)
signature_supporting_files <- c(
  BPS_TARGET_FILE, CHEMBL_TARGET_FILE, SEA_TARGET_FILE, STITCH_TARGET_FILE,
  SUPERPRED_TARGET_FILE, SWISS_TARGET_FILE, OLD_DEG_FILE, GEO_MATRIX_FILE,
  INDEPENDENT_COUNT_FILE, KEGG_SNAPSHOT_FILE,
  file.path(DATA_DIR, "input", "manifests",
            "gdc_tcga_kirc_star_counts_metadata.csv")
)
signature_files <- sort(unique(c(
  signature_tcga_files,
  signature_supporting_files[file.exists(signature_supporting_files)]
)))
signature_info <- file.info(signature_files)
project_root_normalized <- normalizePath(
  REPO_DIR, winslash = "/", mustWork = TRUE
)
signature_paths_normalized <- normalizePath(
  signature_files, winslash = "/", mustWork = TRUE
)
tcga_root_normalized <- normalizePath(
  TCGA_DIR, winslash = "/", mustWork = TRUE
)
signature_relative_paths <- ifelse(
  startsWith(signature_paths_normalized, paste0(project_root_normalized, "/")),
  substring(signature_paths_normalized, nchar(project_root_normalized) + 2L),
  ifelse(
    startsWith(signature_paths_normalized, paste0(tcga_root_normalized, "/")),
    paste0(
      "external_tcga/",
      substring(signature_paths_normalized, nchar(tcga_root_normalized) + 2L)
    ),
    basename(signature_paths_normalized)
  )
)
analysis_input_signature <- data.frame(
  relative_path = signature_relative_paths,
  bytes = as.numeric(signature_info$size),
  modified_unix_time = as.numeric(signature_info$mtime),
  stringsAsFactors = FALSE
)
analysis_input_signature <- analysis_input_signature[
  order(analysis_input_signature$relative_path), , drop = FALSE
]
analysis_input_signature_hash <- digest::digest(
  list(
    pipeline_version = pipeline_version,
    R_version = R.version.string,
    package_versions = package_versions,
    inputs = analysis_input_signature
  ),
  algo = "sha256",
  serialize = TRUE
)
write_csv_utf8_safe <- function(x, filename) {
  utils::write.csv(
    x, file.path(TABLE_DIR, filename), row.names = FALSE,
    fileEncoding = "UTF-8", na = ""
  )
}
write_csv_utf8_safe(
  analysis_input_signature,
  "00B_analysis_input_metadata_signature.csv"
)
writeLines(
  c(
    paste("pipeline_version", pipeline_version, sep = "="),
    paste("analysis_input_signature_sha256", analysis_input_signature_hash,
          sep = "=")
  ),
  file.path(LOG_DIR, "analysis_input_signature.txt"),
  useBytes = TRUE
)

# ---- Figure contract ----------------------------------------------------------
# Core conclusion: corrected TCGA processing identifies computational
# BPS-ccRCC candidates and evaluates classification performance without
# patient leakage.
# Archetype: quantitative grid.
# Backend: R only.
# Export: editable SVG/PDF plus 600-dpi TIFF and 300-dpi PNG preview.
# Reviewer risks addressed: n definition, paired/repeated observations,
# low-expression filtering, batch diagnostics, multiplicity, nested CV,
# class imbalance, confidence intervals, and source-data traceability.

palette_contract <- c(
  NT = "#3182BD",
  TP = "#D24B40",
  neutral = "#8A8A8A",
  down = "#3182BD",
  up = "#D24B40",
  accent = "#E28E2C",
  teal = "#33B5A5"
)

theme_pub <- function(base_size = 9, base_family = "sans") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.35, colour = "black"),
      axis.text = ggplot2::element_text(colour = "black"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.key = ggplot2::element_blank()
    )
}
ggplot2::theme_set(theme_pub())

save_gg <- function(plot, stem, width_mm = 183, height_mm = 120) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  svg_file <- file.path(FIGURE_DIR, paste0(stem, ".svg"))
  pdf_file <- file.path(FIGURE_DIR, paste0(stem, ".pdf"))
  tif_file <- file.path(FIGURE_DIR, paste0(stem, ".tiff"))
  png_file <- file.path(FIGURE_DIR, paste0(stem, ".png"))

  svglite::svglite(svg_file, width = width_in, height = height_in)
  print(plot)
  grDevices::dev.off()

  grDevices::cairo_pdf(pdf_file, width = width_in, height = height_in,
                       family = "sans")
  print(plot)
  grDevices::dev.off()

  ragg::agg_tiff(tif_file, width = width_in, height = height_in,
                 units = "in", res = 600, compression = "lzw")
  print(plot)
  grDevices::dev.off()

  ragg::agg_png(png_file, width = width_in, height = height_in,
                units = "in", res = 300)
  print(plot)
  grDevices::dev.off()
  invisible(c(svg_file, pdf_file, tif_file, png_file))
}

write_csv_utf8 <- function(x, filename, row.names = FALSE) {
  utils::write.csv(
    x,
    file.path(TABLE_DIR, filename),
    row.names = row.names,
    fileEncoding = "UTF-8",
    na = ""
  )
}

format_p_number <- function(p) {
  if (length(p) == 0L || !is.finite(p)) return("not estimable")
  if (p == 0) return("<2.2e-308")
  if (p < 0.01) return(formatC(p, format = "e", digits = 2))
  formatC(p, format = "f", digits = 3)
}

make_p_annotation <- function(result_table, primary_column,
                              primary_label, secondary_column = NULL,
                              secondary_label = NULL) {
  if (is.null(result_table) || nrow(result_table) == 0L ||
      !all(c("gene", primary_column) %in% colnames(result_table))) {
    return(data.frame())
  }
  annotation <- data.frame(
    gene = result_table$gene,
    x = 1.5,
    y = Inf,
    label = paste0(
      primary_label, " = ",
      vapply(result_table[[primary_column]], format_p_number, character(1))
    ),
    stringsAsFactors = FALSE
  )
  if (!is.null(secondary_column) &&
      secondary_column %in% colnames(result_table)) {
    annotation$label <- paste0(
      annotation$label, "\n", secondary_label, " = ",
      vapply(result_table[[secondary_column]], format_p_number, character(1))
    )
  }
  annotation
}

# A completed TCGA/enrichment/ML core is cached before external-validation
# modules. If a downstream optional dataset has a format problem, the next run
# resumes from this checkpoint instead of repeating the long patient-blocked
# model. Delete this checkpoint to force a complete raw-count reanalysis.
resume_from_core_checkpoint <- FALSE
if (file.exists(CORE_CHECKPOINT_FILE)) {
  checkpoint_env <- new.env(parent = emptyenv())
  checkpoint_names <- tryCatch(
    load(CORE_CHECKPOINT_FILE, envir = checkpoint_env),
    error = function(e) character()
  )
  required_checkpoint_names <- c(
    "sample_metadata", "dge", "v", "deg_all", "deg_sig", "bps_targets",
    "bps_overlap", "ml_metrics", "hub_genes", "pipeline_version",
    "analysis_input_signature_hash"
  )
  if (all(required_checkpoint_names %in% checkpoint_names) &&
      identical(checkpoint_env$pipeline_version, pipeline_version) &&
      identical(
        checkpoint_env$analysis_input_signature_hash,
        analysis_input_signature_hash
      )) {
    list2env(as.list(checkpoint_env), envir = .GlobalEnv)
    resume_from_core_checkpoint <- TRUE
    message_stamp("Resuming from valid core checkpoint: ", CORE_CHECKPOINT_FILE)
  } else {
    warning("Ignoring incomplete or stale core checkpoint and recomputing the core")
  }
  rm(checkpoint_env, checkpoint_names, required_checkpoint_names)
}

if (!resume_from_core_checkpoint) {

# Cache the expensive raw-count, duplicateCorrelation and limma stage
# separately from target reconstruction and machine learning. This permits
# output-only or ML fixes without repeating the expression model. The cache is
# accepted only when all analysis-defining expression objects are present.
resume_from_expression_checkpoint <- FALSE
if (file.exists(EXPRESSION_CHECKPOINT_FILE)) {
  expression_checkpoint_env <- new.env(parent = emptyenv())
  expression_checkpoint_names <- tryCatch(
    load(EXPRESSION_CHECKPOINT_FILE, envir = expression_checkpoint_env),
    error = function(e) character()
  )
  required_expression_checkpoint_names <- c(
    "sample_metadata_all", "sample_metadata", "sample_summary",
    "gene_annotation_ensembl", "counts_symbol", "group", "patient_id",
    "design", "dge_unfiltered", "keep_gene", "dge", "v", "tmm_logCPM",
    "duplicate_correlation", "consensus_correlation", "deg_all", "deg_sig",
    "filter_summary", "pca_scores", "batch_association",
    "batch_adjustment_applied", "batch_adjustment_reason",
    "pipeline_version", "analysis_input_signature_hash"
  )
  if (all(required_expression_checkpoint_names %in%
          expression_checkpoint_names) &&
      identical(expression_checkpoint_env$pipeline_version, pipeline_version) &&
      identical(
        expression_checkpoint_env$analysis_input_signature_hash,
        analysis_input_signature_hash
      )) {
    list2env(as.list(expression_checkpoint_env), envir = .GlobalEnv)
    resume_from_expression_checkpoint <- TRUE
    message_stamp(
      "Resuming from valid expression checkpoint: ",
      EXPRESSION_CHECKPOINT_FILE
    )
  } else {
    warning("Ignoring incomplete or stale expression checkpoint and recomputing expression analysis")
  }
  rm(
    expression_checkpoint_env, expression_checkpoint_names,
    required_expression_checkpoint_names
  )
}

if (!resume_from_expression_checkpoint) {

# ---- GDC file-to-sample metadata ---------------------------------------------
gdc_cache_file <- file.path(CACHE_DIR, "gdc_tcga_kirc_star_counts_metadata.csv")
legacy_gdc_cache_file <- file.path(
  DATA_DIR, "input", "manifests",
  "gdc_tcga_kirc_star_counts_metadata.csv"
)

first_nested <- function(x, default = NA_character_) {
  if (is.null(x) || length(x) == 0) return(default)
  value <- x[[1]]
  if (is.null(value) || length(value) == 0) return(default)
  as.character(value)
}

query_gdc_metadata <- function() {
  message_stamp("Querying GDC API for TCGA-KIRC STAR-Counts metadata")
  filters <- list(
    op = "and",
    content = list(
      list(
        op = "in",
        content = list(
          field = "cases.project.project_id",
          value = list("TCGA-KIRC")
        )
      ),
      list(
        op = "in",
        content = list(
          field = "data_type",
          value = list("Gene Expression Quantification")
        )
      ),
      list(
        op = "in",
        content = list(
          field = "analysis.workflow_type",
          value = list("STAR - Counts")
        )
      )
    )
  )

  fields <- paste(
    c(
      "file_id", "file_name", "cases.submitter_id",
      "cases.samples.submitter_id", "cases.samples.sample_type",
      "cases.samples.tissue_type",
      "cases.samples.portions.analytes.aliquots.submitter_id"
    ),
    collapse = ","
  )

  response <- httr::POST(
    "https://api.gdc.cancer.gov/files",
    body = list(
      filters = jsonlite::toJSON(filters, auto_unbox = TRUE),
      format = "JSON",
      fields = fields,
      size = "2000"
    ),
    encode = "form",
    httr::timeout(600)
  )
  httr::stop_for_status(response)
  parsed <- httr::content(response, as = "parsed", type = "application/json",
                          simplifyVector = FALSE)
  hits <- parsed$data$hits
  if (length(hits) == 0) stop("GDC returned no TCGA-KIRC STAR-Counts files")

  rows <- lapply(hits, function(h) {
    case <- if (length(h$cases) > 0) h$cases[[1]] else list()
    sample <- if (length(case$samples) > 0) case$samples[[1]] else list()
    portion <- if (length(sample$portions) > 0) sample$portions[[1]] else list()
    analyte <- if (length(portion$analytes) > 0) portion$analytes[[1]] else list()
    aliquot <- if (length(analyte$aliquots) > 0) analyte$aliquots[[1]] else list()
    data.frame(
      file_id = first_nested(h$file_id),
      file_name = first_nested(h$file_name),
      patient_id = first_nested(case$submitter_id),
      sample_barcode = first_nested(sample$submitter_id),
      sample_type = first_nested(sample$sample_type),
      tissue_type = first_nested(sample$tissue_type),
      aliquot_barcode = first_nested(aliquot$submitter_id),
      stringsAsFactors = FALSE
    )
  })
  metadata <- data.table::rbindlist(rows, fill = TRUE)
  metadata <- as.data.frame(metadata)
  metadata
}

available_gdc_cache <- c(legacy_gdc_cache_file, gdc_cache_file)
available_gdc_cache <- available_gdc_cache[file.exists(available_gdc_cache)]
if (length(available_gdc_cache) > 0L) {
  source_gdc_cache <- available_gdc_cache[[1]]
  gdc_metadata <- utils::read.csv(source_gdc_cache, check.names = FALSE)
  message_stamp("Using cached GDC metadata: ", source_gdc_cache)
  if (!identical(normalizePath(source_gdc_cache, mustWork = TRUE),
                 normalizePath(gdc_cache_file, mustWork = FALSE))) {
    utils::write.csv(gdc_metadata, gdc_cache_file, row.names = FALSE,
                     fileEncoding = "UTF-8")
  }
} else {
  gdc_metadata <- query_gdc_metadata()
  utils::write.csv(gdc_metadata, gdc_cache_file, row.names = FALSE,
                   fileEncoding = "UTF-8")
}

if (nrow(gdc_metadata) != 614L) {
  warning("Expected 614 GDC STAR-Counts files, obtained ", nrow(gdc_metadata),
          ". The script will continue only if every local file maps uniquely.")
}

split_barcode <- function(x, index) {
  vapply(strsplit(ifelse(is.na(x), "", x), "-", fixed = TRUE), function(z) {
    if (length(z) >= index) z[[index]] else NA_character_
  }, character(1))
}
gdc_metadata$tss <- split_barcode(gdc_metadata$patient_id, 2)
gdc_metadata$plate <- split_barcode(gdc_metadata$aliquot_barcode, 6)
gdc_metadata$sequencing_center <- split_barcode(gdc_metadata$aliquot_barcode, 7)

local_files <- list.files(
  TCGA_DIR,
  pattern = "rna_seq\\.augmented_star_gene_counts\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)
if (length(local_files) != 614L) {
  stop("Expected 614 local STAR-Counts files, found ", length(local_files))
}

local_index <- data.frame(
  file_id = basename(dirname(local_files)),
  local_file = normalizePath(local_files, winslash = "/", mustWork = TRUE),
  stringsAsFactors = FALSE
)
sample_metadata_all <- merge(
  local_index,
  gdc_metadata,
  by = "file_id",
  all.x = TRUE,
  sort = FALSE
)
if (anyNA(sample_metadata_all$sample_barcode)) {
  missing_ids <- sample_metadata_all$file_id[is.na(sample_metadata_all$sample_barcode)]
  stop("Missing GDC metadata for local file IDs: ", paste(missing_ids, collapse = ", "))
}

sample_metadata_all <- sample_metadata_all[
  match(local_index$file_id, sample_metadata_all$file_id),
]
sample_metadata <- sample_metadata_all[
  sample_metadata_all$sample_type %in% c("Primary Tumor", "Solid Tissue Normal"),
]
sample_metadata$group <- ifelse(
  sample_metadata$sample_type == "Primary Tumor", "TP", "NT"
)
sample_metadata$group <- factor(sample_metadata$group, levels = c("NT", "TP"))
sample_metadata$profile_id <- paste0(
  sample_metadata$sample_barcode, "__", substr(sample_metadata$file_id, 1, 8)
)
rownames(sample_metadata) <- sample_metadata$profile_id

sample_summary <- data.frame(
  metric = c(
    "Local STAR-Counts files", "Profiles retained (TP + NT)",
    "Primary tumor profiles", "Normal profiles", "Unique sample barcodes",
    "Unique patients", "Patients represented by >1 profile",
    "Duplicated sample barcodes", "Excluded other sample types"
  ),
  value = c(
    length(local_files), nrow(sample_metadata),
    sum(sample_metadata$group == "TP"), sum(sample_metadata$group == "NT"),
    length(unique(sample_metadata$sample_barcode)),
    length(unique(sample_metadata$patient_id)),
    sum(table(sample_metadata$patient_id) > 1),
    sum(table(sample_metadata$sample_barcode) > 1),
    nrow(sample_metadata_all) - nrow(sample_metadata)
  )
)
write_csv_utf8(sample_metadata_all, "01_sample_metadata_all_614.csv")
write_csv_utf8(sample_metadata, "02_sample_metadata_analysis.csv")
write_csv_utf8(sample_summary, "03_sample_summary.csv")
message_stamp(
  "Retained profiles: TP=", sum(sample_metadata$group == "TP"),
  ", NT=", sum(sample_metadata$group == "NT"),
  ", unique patients=", length(unique(sample_metadata$patient_id))
)

# ---- Read and aggregate STAR counts ------------------------------------------
read_star_counts <- function(path) {
  x <- data.table::fread(
    path,
    skip = 1,
    select = c("gene_id", "gene_name", "gene_type", "unstranded"),
    data.table = FALSE,
    showProgress = FALSE
  )
  x <- x[!grepl("^N_", x$gene_id), , drop = FALSE]
  x
}

message_stamp("Reading local STAR-Counts files")
first_counts <- read_star_counts(sample_metadata$local_file[[1]])
gene_annotation_ensembl <- first_counts[, c("gene_id", "gene_name", "gene_type")]
counts_ensembl <- matrix(
  0L,
  nrow = nrow(first_counts),
  ncol = nrow(sample_metadata),
  dimnames = list(first_counts$gene_id, sample_metadata$profile_id)
)
counts_ensembl[, 1] <- as.integer(first_counts$unstranded)

for (i in seq.int(2L, nrow(sample_metadata))) {
  current <- read_star_counts(sample_metadata$local_file[[i]])
  if (!identical(current$gene_id, first_counts$gene_id)) {
    idx <- match(first_counts$gene_id, current$gene_id)
    if (anyNA(idx)) stop("Gene IDs differ in file: ", sample_metadata$local_file[[i]])
    current <- current[idx, , drop = FALSE]
  }
  counts_ensembl[, i] <- as.integer(current$unstranded)
  if (i %% 50L == 0L || i == nrow(sample_metadata)) {
    message_stamp("Read ", i, "/", nrow(sample_metadata), " profiles")
  }
}
rm(first_counts, current)
invisible(gc())

valid_symbol <- !is.na(gene_annotation_ensembl$gene_name) &
  nzchar(trimws(gene_annotation_ensembl$gene_name))
counts_symbol <- rowsum(
  counts_ensembl[valid_symbol, , drop = FALSE],
  group = gene_annotation_ensembl$gene_name[valid_symbol],
  reorder = TRUE
)
storage.mode(counts_symbol) <- "integer"
rm(counts_ensembl)
invisible(gc())

# ---- Differential expression -------------------------------------------------
group <- sample_metadata$group
patient_id <- factor(sample_metadata$patient_id)
design <- model.matrix(~ 0 + group)
colnames(design) <- c("NT", "TP")
rownames(design) <- sample_metadata$profile_id

dge_unfiltered <- edgeR::DGEList(counts = counts_symbol, group = group)
keep_gene <- edgeR::filterByExpr(dge_unfiltered, design = design)
dge <- dge_unfiltered[keep_gene, , keep.lib.sizes = FALSE]
dge <- edgeR::calcNormFactors(dge, method = "TMM")

filter_summary <- data.frame(
  metric = c(
    "Genes before low-expression filtering",
    "Genes retained by edgeR::filterByExpr",
    "Genes removed by low-expression filtering"
  ),
  value = c(nrow(dge_unfiltered), sum(keep_gene), sum(!keep_gene))
)
write_csv_utf8(filter_summary, "04_gene_filter_summary.csv")

v_initial <- limma::voom(dge, design, plot = FALSE, save.plot = TRUE)

# Patient-level correlation handles partially paired tumor-normal profiles and
# repeated profiles from the same participant without pretending independence.
duplicate_correlation <- tryCatch(
  limma::duplicateCorrelation(v_initial, design, block = patient_id),
  error = function(e) {
    warning("duplicateCorrelation failed: ", conditionMessage(e))
    list(consensus.correlation = 0)
  }
)
consensus_correlation <- duplicate_correlation$consensus.correlation
if (!is.finite(consensus_correlation)) consensus_correlation <- 0

v <- limma::voom(
  dge,
  design,
  plot = FALSE,
  save.plot = TRUE,
  block = patient_id,
  correlation = consensus_correlation
)
tmm_logCPM <- edgeR::cpm(
  dge, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 0.5
)
fit <- limma::lmFit(
  v,
  design,
  block = patient_id,
  correlation = consensus_correlation
)
contrast_matrix <- limma::makeContrasts(TP - NT, levels = design)
fit2 <- limma::contrasts.fit(fit, contrast_matrix)
fit2 <- limma::eBayes(fit2, robust = TRUE)

deg_all <- limma::topTable(
  fit2,
  number = Inf,
  adjust.method = "BH",
  sort.by = "P",
  confint = 0.95
)
deg_all$gene <- rownames(deg_all)
deg_all$direction <- "Not significant"
deg_all$direction[deg_all$adj.P.Val < 0.05 & deg_all$logFC > 1] <- "Up"
deg_all$direction[deg_all$adj.P.Val < 0.05 & deg_all$logFC < -1] <- "Down"
deg_sig <- deg_all[
  deg_all$adj.P.Val < 0.05 & abs(deg_all$logFC) > 1,
  , drop = FALSE
]

write_csv_utf8(deg_all, "05_DEG_all_filterByExpr_TMM_voom_limma.csv")
write_csv_utf8(deg_sig, "06_DEG_significant_absLogFC1_BH005.csv")

# ---- Batch and unwanted-variation diagnostics --------------------------------
pca <- stats::prcomp(t(v$E), center = TRUE, scale. = FALSE)
pc_var <- 100 * (pca$sdev^2 / sum(pca$sdev^2))
pca_scores <- data.frame(
  profile_id = rownames(pca$x),
  pca$x[, seq_len(min(10L, ncol(pca$x))), drop = FALSE],
  sample_metadata[rownames(pca$x), c(
    "sample_barcode", "patient_id", "group", "tss", "plate",
    "sequencing_center"
  )],
  check.names = FALSE
)

pc_batch_association <- function(scores, variables, n_pc = 5L) {
  out <- list()
  for (variable in variables) {
    z <- scores[[variable]]
    for (j in seq_len(min(n_pc, sum(grepl("^PC", names(scores)))))) {
      pc_name <- paste0("PC", j)
      ok <- is.finite(scores[[pc_name]]) & !is.na(z) & nzchar(as.character(z))
      zz <- droplevels(factor(z[ok]))
      if (sum(ok) < 10L || nlevels(zz) < 2L || nlevels(zz) >= sum(ok)) next
      model <- stats::lm(scores[[pc_name]][ok] ~ zz)
      a <- stats::anova(model)
      ss_total <- sum(a$`Sum Sq`, na.rm = TRUE)
      eta2 <- if (ss_total > 0) a$`Sum Sq`[[1]] / ss_total else NA_real_
      out[[length(out) + 1L]] <- data.frame(
        variable = variable,
        PC = pc_name,
        n = sum(ok),
        levels = nlevels(zz),
        eta_squared = eta2,
        p_value = a$`Pr(>F)`[[1]],
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(out) == 0) return(data.frame())
  result <- do.call(rbind, out)
  result$BH_adjusted_p <- p.adjust(result$p_value, method = "BH")
  result
}

batch_association <- pc_batch_association(
  pca_scores,
  variables = c("group", "tss", "plate", "sequencing_center"),
  n_pc = 5L
)
write_csv_utf8(pca_scores, "07_PCA_scores_with_metadata.csv")
write_csv_utf8(batch_association, "08_PC_batch_association_diagnostics.csv")

# No automatic ComBat/SVA/RUV adjustment is applied. The common GDC processing
# removes pipeline heterogeneity, while PCA/MDS and technical metadata are used
# to diagnose residual variation. Blind correction could remove biology when a
# technical variable is confounded with tissue group. Diagnostics are saved for
# explicit author review.
batch_adjustment_applied <- FALSE
batch_adjustment_reason <- paste(
  "All profiles were generated by the harmonized GDC STAR-Counts workflow.",
  "Residual variation was assessed by PCA/MDS and available TSS, plate and",
  "sequencing-center metadata. No automatic expression-matrix batch removal",
  "was applied because technical factors can be confounded with tissue group."
)

top_tss <- names(sort(table(pca_scores$tss), decreasing = TRUE))[1:min(10, length(unique(pca_scores$tss)))]
pca_scores$tss_plot <- ifelse(pca_scores$tss %in% top_tss, pca_scores$tss, "Other")

p_pca_group <- ggplot(pca_scores, aes(PC1, PC2, colour = group)) +
  geom_point(alpha = 0.75, size = 1.7) +
  scale_colour_manual(values = palette_contract[c("NT", "TP")]) +
  labs(
    title = "PCA by tissue group",
    x = sprintf("PC1 (%.1f%%)", pc_var[[1]]),
    y = sprintf("PC2 (%.1f%%)", pc_var[[2]]),
    colour = "Group"
  )

p_pca_tss <- ggplot(pca_scores, aes(PC1, PC2, colour = tss_plot, shape = group)) +
  geom_point(alpha = 0.75, size = 1.5) +
  labs(
    title = "PCA annotated by tissue source site",
    x = sprintf("PC1 (%.1f%%)", pc_var[[1]]),
    y = sprintf("PC2 (%.1f%%)", pc_var[[2]]),
    colour = "TSS", shape = "Group"
  ) +
  guides(shape = "none") +
  theme(legend.position = "right")

lib_df <- data.frame(
  group = group,
  library_size_million = dge$samples$lib.size / 1e6
)
p_library <- ggplot(lib_df, aes(group, library_size_million, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.65, colour = NA) +
  geom_boxplot(width = 0.16, outlier.size = 0.5) +
  scale_fill_manual(values = palette_contract[c("NT", "TP")]) +
  labs(title = "Library sizes after gene filtering", x = NULL,
       y = "Library size (million reads)") +
  theme(legend.position = "none")

voom_df <- data.frame(
  sqrt_standard_deviation = v$voom.xy$y,
  mean_log2_count = v$voom.xy$x
)
voom_line <- data.frame(
  mean_log2_count = v$voom.line$x,
  sqrt_standard_deviation = v$voom.line$y
)
p_voom <- ggplot(voom_df, aes(mean_log2_count, sqrt_standard_deviation)) +
  geom_point(alpha = 0.18, size = 0.35, colour = palette_contract[["neutral"]]) +
  geom_line(data = voom_line, colour = palette_contract[["accent"]],
            linewidth = 0.7) +
  labs(title = "Voom mean-variance trend", x = "Mean log2 count",
       y = "Square-root residual SD")

qc_figure <- (p_pca_group | p_library) / (p_pca_tss | p_voom) +
  patchwork::plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 11))
save_gg(qc_figure, "Figure_QC_PCA_batch_voom", 183, 150)

mds <- limma::plotMDS(v$E, top = 500, gene.selection = "common", plot = FALSE)
mds_df <- data.frame(
  Dim1 = mds$x,
  Dim2 = mds$y,
  group = group,
  patient_id = sample_metadata$patient_id
)
p_mds <- ggplot(mds_df, aes(Dim1, Dim2, colour = group)) +
  geom_point(alpha = 0.75, size = 1.7) +
  scale_colour_manual(values = palette_contract[c("NT", "TP")]) +
  labs(title = "MDS of the 500 most variable genes", x = "Dimension 1",
       y = "Dimension 2", colour = "Group")
save_gg(p_mds, "Figure_MDS_group", 89, 80)

volcano_df <- deg_all
volcano_df$minus_log10_fdr <- -log10(pmax(volcano_df$adj.P.Val, .Machine$double.xmin))
p_volcano <- ggplot(volcano_df, aes(logFC, minus_log10_fdr, colour = direction)) +
  geom_point(alpha = 0.55, size = 0.7) +
  geom_vline(xintercept = c(-1, 1), linetype = 2, linewidth = 0.35) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.35) +
  scale_colour_manual(values = c(
    "Down" = palette_contract[["down"]],
    "Not significant" = "#BDBDBD",
    "Up" = palette_contract[["up"]]
  )) +
  labs(
    title = "TCGA-KIRC differential expression",
    subtitle = sprintf("BH FDR < 0.05 and |log2FC| > 1; n = %d genes", nrow(deg_sig)),
    x = "log2 fold change (TP - NT)",
    y = "-log10(BH-adjusted P)",
    colour = NULL
  )
save_gg(p_volcano, "Figure_DEG_volcano", 89, 85)

# Heatmap of the top 50 DEGs. Values are voom log2-CPM, row-standardized only
# for visualization; sample labels are hidden because n = 613 profiles.
top_heat_genes <- head(deg_sig$gene, 50L)
heat_mat <- v$E[top_heat_genes, , drop = FALSE]
heat_annotation <- data.frame(Group = group)
rownames(heat_annotation) <- colnames(heat_mat)
heat_colors <- list(Group = palette_contract[c("NT", "TP")])

save_pheatmap <- function(filename, device = c("pdf", "tiff", "png")) {
  device <- match.arg(device)
  path <- file.path(FIGURE_DIR, paste0(filename, ".", ifelse(device == "tiff", "tiff", device)))
  if (device == "pdf") {
    grDevices::cairo_pdf(path, width = 7.2, height = 7.0, family = "sans")
  } else if (device == "tiff") {
    ragg::agg_tiff(path, width = 7.2, height = 7.0, units = "in", res = 600,
                   compression = "lzw")
  } else {
    ragg::agg_png(path, width = 7.2, height = 7.0, units = "in", res = 300)
  }
  pheatmap::pheatmap(
    heat_mat,
    scale = "row",
    show_colnames = FALSE,
    show_rownames = TRUE,
    annotation_col = heat_annotation,
    annotation_colors = heat_colors,
    border_color = NA,
    fontsize = 7,
    fontsize_row = 6,
    clustering_method = "ward.D2"
  )
  grDevices::dev.off()
  invisible(path)
}
save_pheatmap("Figure_DEG_top50_heatmap", "pdf")
save_pheatmap("Figure_DEG_top50_heatmap", "tiff")
save_pheatmap("Figure_DEG_top50_heatmap", "png")

expression_checkpoint_objects <- c(
  "sample_metadata_all", "sample_metadata", "sample_summary",
  "gene_annotation_ensembl", "counts_symbol", "group", "patient_id",
  "design", "dge_unfiltered", "keep_gene", "dge", "v", "tmm_logCPM",
  "duplicate_correlation", "consensus_correlation", "deg_all", "deg_sig",
  "filter_summary", "pca", "pca_scores", "batch_association",
  "batch_adjustment_applied", "batch_adjustment_reason",
  "pipeline_version", "analysis_input_signature",
  "analysis_input_signature_hash"
)
expression_checkpoint_objects <- expression_checkpoint_objects[
  vapply(expression_checkpoint_objects, function(object_name) {
    exists(object_name, envir = .GlobalEnv, inherits = FALSE)
  }, logical(1))
]
save(
  list = expression_checkpoint_objects,
  file = EXPRESSION_CHECKPOINT_FILE,
  compress = "gzip"
)
message_stamp("Saved expression checkpoint: ", EXPRESSION_CHECKPOINT_FILE)
}

# ---- Reproducible BPS target universe and DEG overlap -------------------------
# The primary target universe is rebuilt from the supplied raw database exports.
# Database-specific scores are used only within their source database. ChEMBL,
# STITCH and similarity/model-based predictions are retained as distinct evidence
# classes and are not treated as experimentally equivalent. The original 210-gene
# file is retained only as an audit comparator because it contains legacy UniProt
# entry names (for example, 5NTD and CAH2) rather than current gene symbols.

required_target_files <- c(
  CHEMBL_TARGET_FILE, SEA_TARGET_FILE, STITCH_TARGET_FILE,
  SUPERPRED_TARGET_FILE, SWISS_TARGET_FILE
)
missing_target_files <- required_target_files[!file.exists(required_target_files)]
if (length(missing_target_files) > 0L) {
  stop("Missing raw BPS target files: ", paste(missing_target_files, collapse = "; "))
}

map_uniprot_to_symbol <- function(ids) {
  ids <- unique(trimws(as.character(ids)))
  ids <- ids[nzchar(ids)]
  mapped <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = ids,
    keytype = "UNIPROT",
    column = "SYMBOL",
    multiVals = "first"
  )
  data.frame(
    source_id = names(mapped),
    gene = unname(mapped),
    stringsAsFactors = FALSE
  )
}

chembl_map <- map_uniprot_to_symbol(readLines(
  CHEMBL_TARGET_FILE, warn = FALSE, encoding = "UTF-8"
))
chembl_map <- chembl_map[!is.na(chembl_map$gene) & nzchar(chembl_map$gene), ]
chembl_records <- data.frame(
  database = "ChEMBL",
  gene = chembl_map$gene,
  source_id = chembl_map$source_id,
  evidence_type = "database-annotated association",
  score_name = "not supplied in exported identifier list",
  score_value = NA_real_,
  model_accuracy = NA_real_,
  inclusion_rule = "human UniProt identifiers mapped to official gene symbols",
  high_confidence = FALSE,
  stringsAsFactors = FALSE
)

sea_lines <- readLines(SEA_TARGET_FILE, warn = FALSE, encoding = "UTF-8")
sea_lines <- sub('^"', "", sea_lines)
sea_lines <- sub('"$', "", sea_lines)
sea_raw <- utils::read.csv(
  text = paste(c(sea_lines, ""), collapse = "\n"),
  check.names = FALSE,
  quote = "",
  stringsAsFactors = FALSE
)
sea_keep <- grepl("_HUMAN$", sea_raw[["Target ID"]]) &
  is.finite(sea_raw[["P-Value"]]) & sea_raw[["P-Value"]] <= 0.05 &
  !is.na(sea_raw$Name) & nzchar(trimws(sea_raw$Name))
sea_filtered <- sea_raw[sea_keep, , drop = FALSE]
sea_records <- data.frame(
  database = "SEA",
  gene = trimws(sea_filtered$Name),
  source_id = sea_filtered[["Target ID"]],
  evidence_type = "similarity-based prediction",
  score_name = "P-value",
  score_value = as.numeric(sea_filtered[["P-Value"]]),
  model_accuracy = NA_real_,
  inclusion_rule = "Homo sapiens and SEA P-value <= 0.05",
  high_confidence = TRUE,
  stringsAsFactors = FALSE
)

stitch_raw <- utils::read.delim(
  STITCH_TARGET_FILE, check.names = FALSE, stringsAsFactors = FALSE
)
stitch_direct <-
  (tolower(stitch_raw[["#node1"]]) == "bisphenol s" |
     tolower(stitch_raw$node2) == "bisphenol s") &
  is.finite(stitch_raw$combined_score) & stitch_raw$combined_score >= 0.700
stitch_filtered <- stitch_raw[stitch_direct, , drop = FALSE]
stitch_gene <- ifelse(
  tolower(stitch_filtered[["#node1"]]) == "bisphenol s",
  stitch_filtered$node2,
  stitch_filtered[["#node1"]]
)
stitch_records <- data.frame(
  database = "STITCH",
  gene = trimws(stitch_gene),
  source_id = ifelse(
    tolower(stitch_filtered[["#node1"]]) == "bisphenol s",
    stitch_filtered$node1_external_id,
    stitch_filtered$node2_external_id
  ),
  evidence_type = "direct chemical-protein database edge",
  score_name = "combined score",
  score_value = stitch_filtered$combined_score,
  model_accuracy = NA_real_,
  inclusion_rule = "direct BPS-protein edge and combined score >= 0.700",
  high_confidence = TRUE,
  stringsAsFactors = FALSE
)

superpred_raw <- utils::read.csv(
  SUPERPRED_TARGET_FILE, check.names = FALSE, stringsAsFactors = FALSE
)
superpred_probability <- as.numeric(sub(
  "%", "", superpred_raw$Probability, fixed = TRUE
)) / 100
superpred_model_accuracy <- as.numeric(sub(
  "%", "", superpred_raw$`Model accuracy`, fixed = TRUE
)) / 100
superpred_keep <- is.finite(superpred_probability) & superpred_probability >= 0.50
superpred_map <- map_uniprot_to_symbol(
  superpred_raw$`UniProt ID`[superpred_keep]
)
superpred_meta <- data.frame(
  source_id = superpred_raw$`UniProt ID`[superpred_keep],
  score_value = superpred_probability[superpred_keep],
  model_accuracy = superpred_model_accuracy[superpred_keep],
  stringsAsFactors = FALSE
)
superpred_map <- merge(superpred_map, superpred_meta, by = "source_id", all.x = TRUE)
superpred_map <- superpred_map[!is.na(superpred_map$gene) & nzchar(superpred_map$gene), ]
superpred_records <- data.frame(
  database = "SuperPred",
  gene = superpred_map$gene,
  source_id = superpred_map$source_id,
  evidence_type = "model-based prediction",
  score_name = "probability",
  score_value = superpred_map$score_value,
  model_accuracy = superpred_map$model_accuracy,
  inclusion_rule = "prediction probability >= 0.50",
  high_confidence = superpred_map$model_accuracy >= 0.80,
  stringsAsFactors = FALSE
)

swiss_raw <- utils::read.csv(
  SWISS_TARGET_FILE, check.names = FALSE, stringsAsFactors = FALSE
)
swiss_records <- data.frame(
  database = "SwissTargetPrediction",
  gene = trimws(swiss_raw$`Common name`),
  source_id = swiss_raw$`Uniprot ID`,
  evidence_type = "similarity-based prediction",
  score_name = "probability",
  score_value = as.numeric(swiss_raw$`Probability*`),
  model_accuracy = NA_real_,
  inclusion_rule = "all 100 ranked Homo sapiens predictions retained",
  high_confidence = as.numeric(swiss_raw$`Probability*`) > 0,
  stringsAsFactors = FALSE
)

target_source_records <- unique(rbind(
  chembl_records, sea_records, stitch_records,
  superpred_records, swiss_records
))
target_source_records <- target_source_records[
  !is.na(target_source_records$gene) & nzchar(target_source_records$gene),
]

target_database_counts <- do.call(rbind, lapply(
  split(target_source_records, target_source_records$database),
  function(z) data.frame(
    database = z$database[[1]],
    records_after_database_filter = nrow(z),
    unique_genes_after_database_filter = length(unique(z$gene)),
    unique_genes_in_high_confidence_sensitivity =
      length(unique(z$gene[z$high_confidence])),
    stringsAsFactors = FALSE
  )
))
bps_targets <- sort(unique(target_source_records$gene))
bps_targets_high_confidence <- sort(unique(
  target_source_records$gene[target_source_records$high_confidence]
))
target_database_counts$integrated_unique_genes <- length(bps_targets)

bps_targets_input <- if (file.exists(BPS_TARGET_FILE)) {
  trimws(readLines(BPS_TARGET_FILE, warn = FALSE, encoding = "UTF-8"))
} else {
  character()
}
bps_targets_input <- unique(bps_targets_input[nzchar(bps_targets_input)])
target_legacy_comparison <- data.frame(
  gene_or_entry = union(bps_targets_input, bps_targets),
  in_legacy_210_file = union(bps_targets_input, bps_targets) %in% bps_targets_input,
  in_reconstructed_official_symbol_set = union(bps_targets_input, bps_targets) %in% bps_targets,
  stringsAsFactors = FALSE
)
target_cross_database_frequency <- as.data.frame(table(
  unique(target_source_records[, c("database", "gene")])$gene
), stringsAsFactors = FALSE)
colnames(target_cross_database_frequency) <- c("gene", "database_count")
target_cross_database_frequency <- target_cross_database_frequency[
  order(-target_cross_database_frequency$database_count,
        target_cross_database_frequency$gene),
]

bps_overlap <- deg_sig[deg_sig$gene %in% bps_targets, , drop = FALSE]
bps_overlap <- bps_overlap[order(bps_overlap$adj.P.Val), , drop = FALSE]
bps_overlap_high_confidence <- deg_sig[
  deg_sig$gene %in% bps_targets_high_confidence, , drop = FALSE
]
bps_overlap_high_confidence <- bps_overlap_high_confidence[
  order(bps_overlap_high_confidence$adj.P.Val), , drop = FALSE
]

write_csv_utf8(target_source_records, "09A_BPS_target_source_records.csv")
write_csv_utf8(target_database_counts, "09B_BPS_target_database_counts.csv")
write_csv_utf8(data.frame(gene = bps_targets),
               "09C_BPS_candidate_targets_reconstructed.csv")
write_csv_utf8(data.frame(gene = bps_targets_high_confidence),
               "09D_BPS_candidate_targets_high_confidence_sensitivity.csv")
write_csv_utf8(target_legacy_comparison, "09E_legacy210_vs_reconstructed_targets.csv")
write_csv_utf8(target_cross_database_frequency,
               "09F_BPS_target_cross_database_frequency.csv")
write_csv_utf8(bps_overlap, "10A_BPS_ccRCC_DEG_overlap_reconstructed.csv")
write_csv_utf8(bps_overlap_high_confidence,
               "10B_BPS_ccRCC_DEG_overlap_high_confidence_sensitivity.csv")

old_deg <- if (file.exists(OLD_DEG_FILE)) {
  unique(trimws(readLines(OLD_DEG_FILE, warn = FALSE, encoding = "UTF-8")))
} else character()
old_deg <- old_deg[nzchar(old_deg)]
old_overlap <- intersect(old_deg, bps_targets)
overlap_stability <- data.frame(
  metric = c(
    "Old DEG count", "Corrected DEG count", "Old BPS-DEG overlap count",
    "Corrected BPS-DEG overlap count", "Overlap retained from old analysis",
    "New overlap candidates under corrected analysis"
  ),
  value = c(
    length(old_deg), nrow(deg_sig), length(old_overlap), nrow(bps_overlap),
    length(intersect(old_overlap, bps_overlap$gene)),
    length(setdiff(bps_overlap$gene, old_overlap))
  )
)
write_csv_utf8(overlap_stability, "11_old_vs_corrected_overlap_summary.csv")
write_csv_utf8(
  data.frame(
    gene = union(old_overlap, bps_overlap$gene),
    old_overlap = union(old_overlap, bps_overlap$gene) %in% old_overlap,
    corrected_overlap = union(old_overlap, bps_overlap$gene) %in% bps_overlap$gene
  ),
  "12_old_vs_corrected_overlap_gene_status.csv"
)
message_stamp(
  "Corrected DEGs=", nrow(deg_sig),
  "; BPS-DEG overlap=", nrow(bps_overlap)
)

# ---- GO and KEGG enrichment with tested-gene universe -------------------------
tested_symbols <- rownames(v$E)
tested_entrez <- AnnotationDbi::mapIds(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = tested_symbols,
  keytype = "SYMBOL",
  column = "ENTREZID",
  multiVals = "first"
)
overlap_entrez <- AnnotationDbi::mapIds(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = bps_overlap$gene,
  keytype = "SYMBOL",
  column = "ENTREZID",
  multiVals = "first"
)
tested_entrez <- unique(stats::na.omit(unname(tested_entrez)))
overlap_entrez <- unique(stats::na.omit(unname(overlap_entrez)))

go_result <- NULL
kegg_result <- NULL
go_table <- data.frame()
kegg_table <- data.frame()
kegg_snapshot <- NULL
kegg_snapshot_sha256 <- NA_character_
kegg_snapshot_validation <- data.frame()

if (length(overlap_entrez) >= 5L) {
  go_result <- clusterProfiler::enrichGO(
    gene = overlap_entrez,
    universe = tested_entrez,
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  go_table <- as.data.frame(go_result)
  write_csv_utf8(go_table, "13_GO_enrichment_tested_gene_universe.csv")

  if (!file.exists(KEGG_SNAPSHOT_FILE)) {
    warning(
      "KEGG enrichment was skipped because the licensed frozen snapshot was ",
      "not supplied. Set KIRC_BPS_KEGG_SNAPSHOT to an authorised local copy."
    )
  } else {
  kegg_snapshot_sha256 <- digest::digest(
    file = KEGG_SNAPSHOT_FILE, algo = "sha256"
  )
  if (!identical(kegg_snapshot_sha256, KEGG_SNAPSHOT_SHA256)) {
    stop(
      "KEGG snapshot SHA-256 mismatch. Expected ", KEGG_SNAPSHOT_SHA256,
      "; observed ", kegg_snapshot_sha256
    )
  }
  kegg_snapshot <- readRDS(KEGG_SNAPSHOT_FILE)
  required_kegg_fields <- c(
    "schema_version", "database", "organism", "retrieved_at", "endpoints",
    "source_sha256", "term2gene", "term2name", "term_categories"
  )
  if (!all(required_kegg_fields %in% names(kegg_snapshot)) ||
      !identical(kegg_snapshot$database, "KEGG") ||
      !identical(kegg_snapshot$organism, "hsa")) {
    stop("Invalid or incompatible frozen KEGG snapshot: ", KEGG_SNAPSHOT_FILE)
  }
  kegg_result <- clusterProfiler::enricher(
    gene = overlap_entrez,
    universe = tested_entrez,
    TERM2GENE = kegg_snapshot$term2gene,
    TERM2NAME = kegg_snapshot$term2name,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    minGSSize = 10,
    maxGSSize = 500
  )
  kegg_table <- as.data.frame(kegg_result)
  if (nrow(kegg_table) > 0L) {
    kegg_table <- merge(
      kegg_snapshot$term_categories,
      kegg_table,
      by.x = "term",
      by.y = "ID",
      all.y = TRUE,
      sort = FALSE
    )
    names(kegg_table)[names(kegg_table) == "term"] <- "ID"
    kegg_table <- kegg_table[
      order(kegg_table$p.adjust, kegg_table$pvalue, kegg_table$ID),
      c(
        "category", "subcategory", "ID", "Description", "GeneRatio",
        "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"
      ),
      drop = FALSE
    ]
    rownames(kegg_table) <- NULL
  }
  snapshot_input_matches <-
    !is.null(kegg_snapshot$validation_input) &&
    setequal(overlap_entrez, kegg_snapshot$validation_input$overlap_entrez) &&
    setequal(tested_entrez, kegg_snapshot$validation_input$tested_entrez)
  snapshot_result_matches <- NA
  if (snapshot_input_matches && !is.null(kegg_snapshot$validation_result)) {
    expected_kegg <- kegg_snapshot$validation_result
    shared_numeric <- intersect(
      c("pvalue", "p.adjust", "qvalue", "Count"), names(kegg_table)
    )
    expected_kegg <- expected_kegg[
      match(kegg_table$ID, expected_kegg$ID), names(kegg_table), drop = FALSE
    ]
    snapshot_result_matches <-
      identical(kegg_table$ID, expected_kegg$ID) &&
      isTRUE(all.equal(
        kegg_table[, setdiff(names(kegg_table), shared_numeric), drop = FALSE],
        expected_kegg[, setdiff(names(kegg_table), shared_numeric), drop = FALSE],
        check.attributes = FALSE
      )) &&
      all(vapply(shared_numeric, function(nm) {
        isTRUE(all.equal(
          as.numeric(kegg_table[[nm]]), as.numeric(expected_kegg[[nm]]),
          tolerance = 1e-12, check.attributes = FALSE
        ))
      }, logical(1)))
    if (!isTRUE(snapshot_result_matches)) {
      stop("Frozen KEGG enrichment failed its embedded validation check.")
    }
  }
  kegg_snapshot_validation <- data.frame(
    snapshot_file = normalizePath(
      KEGG_SNAPSHOT_FILE, winslash = "/", mustWork = TRUE
    ),
    snapshot_sha256 = kegg_snapshot_sha256,
    retrieved_at = kegg_snapshot$retrieved_at,
    term2gene_rows = nrow(kegg_snapshot$term2gene),
    term2name_rows = nrow(kegg_snapshot$term2name),
    overlap_entrez_n = length(overlap_entrez),
    tested_entrez_n = length(tested_entrez),
    significant_pathways_n = nrow(kegg_table),
    embedded_input_matches = snapshot_input_matches,
    embedded_result_matches = snapshot_result_matches,
    stringsAsFactors = FALSE
  )
  write_csv_utf8(
    kegg_snapshot_validation,
    "14A_KEGG_frozen_snapshot_validation.csv"
  )
  write_csv_utf8(kegg_table, "14_KEGG_enrichment_tested_gene_universe.csv")
  message_stamp(
    "KEGG enrichment used frozen hsa snapshot retrieved ",
    kegg_snapshot$retrieved_at, "; significant pathways=", nrow(kegg_table)
  )
  }
}

if (nrow(go_table) > 0) {
  go_plot_df <- do.call(rbind, lapply(split(go_table, go_table$ONTOLOGY), function(z) {
    z <- z[order(z$p.adjust), , drop = FALSE]
    head(z, 10L)
  }))
  go_plot_df$Description <- factor(
    go_plot_df$Description,
    levels = rev(unique(go_plot_df$Description[order(go_plot_df$p.adjust)]))
  )
  p_go <- ggplot(go_plot_df, aes(-log10(p.adjust), Description,
                                 size = Count, colour = ONTOLOGY)) +
    geom_point(alpha = 0.85) +
    facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
    labs(title = "GO enrichment of corrected BPS-ccRCC candidates",
         x = "-log10(BH-adjusted P)", y = NULL, size = "Gene count",
         colour = "Ontology") +
    theme(strip.text.y = element_text(angle = 0), legend.position = "right")
  save_gg(p_go, "Figure_GO_enrichment", 183, 150)
}

if (nrow(kegg_table) > 0) {
  kegg_plot_df <- head(kegg_table[order(kegg_table$p.adjust), , drop = FALSE], 15L)
  kegg_plot_df$Description <- factor(
    kegg_plot_df$Description,
    levels = rev(kegg_plot_df$Description)
  )
  p_kegg <- ggplot(kegg_plot_df, aes(-log10(p.adjust), Description, size = Count)) +
    geom_point(colour = palette_contract[["teal"]], alpha = 0.85) +
    labs(title = "KEGG pathway enrichment",
         x = "-log10(BH-adjusted P)", y = NULL, size = "Gene count")
  save_gg(p_kegg, "Figure_KEGG_enrichment", 120, 100)
}

# ---- Patient-grouped nested cross-validation machine learning ----------------
# Candidate features are the reconstructed BPS candidates that pass the
# expression filter; the full-cohort DEG result is not used for nested-CV input.
# All profiles from one patient stay in the same outer and inner fold. SVM-RFE
# ranking, feature-count selection and parameter tuning are repeated inside
# training partitions only.

make_balanced_group_folds <- function(patient, outcome, k = 5L, seed = 1L) {
  patient <- as.character(patient)
  outcome <- factor(outcome, levels = c("NT", "TP"))
  patients <- unique(patient)
  p_counts <- t(vapply(patients, function(id) {
    tab <- table(factor(outcome[patient == id], levels = levels(outcome)))
    as.numeric(tab)
  }, numeric(2)))
  colnames(p_counts) <- levels(outcome)
  rownames(p_counts) <- patients

  set.seed(seed)
  jitter <- runif(nrow(p_counts))
  ordering <- order(rowSums(p_counts > 0), rowSums(p_counts), jitter,
                    decreasing = TRUE)
  fold_totals <- matrix(0, nrow = k, ncol = 2,
                        dimnames = list(seq_len(k), levels(outcome)))
  assignment <- integer(length(patients))
  names(assignment) <- patients
  target <- colSums(p_counts) / k
  target[target == 0] <- 1

  for (idx in ordering) {
    cnt <- p_counts[idx, ]
    scores <- vapply(seq_len(k), function(fold) {
      sum(((fold_totals[fold, ] + cnt) / target)^2) +
        0.02 * sum(fold_totals[fold, ])
    }, numeric(1))
    best <- which(scores == min(scores))
    chosen <- sample(best, 1)
    assignment[patients[[idx]]] <- chosen
    fold_totals[chosen, ] <- fold_totals[chosen, ] + cnt
  }
  unname(assignment[patient])
}

scale_train_test <- function(x_train, x_test) {
  center <- colMeans(x_train)
  scale <- apply(x_train, 2, stats::sd)
  scale[!is.finite(scale) | scale == 0] <- 1
  list(
    train = sweep(sweep(x_train, 2, center, "-"), 2, scale, "/"),
    test = sweep(sweep(x_test, 2, center, "-"), 2, scale, "/"),
    center = center,
    scale = scale
  )
}

class_weights <- function(y) {
  tab <- table(y)
  w <- length(y) / (length(tab) * tab)
  stats::setNames(as.numeric(w), names(tab))
}

fit_linear_svm <- function(x, y, cost = 1, probability = TRUE) {
  e1071::svm(
    x = x,
    y = factor(y, levels = c("NT", "TP")),
    type = "C-classification",
    kernel = "linear",
    cost = cost,
    scale = FALSE,
    class.weights = class_weights(factor(y, levels = c("NT", "TP"))),
    probability = probability
  )
}

predict_tp_probability <- function(model, newx) {
  pred <- stats::predict(model, newx, probability = TRUE)
  probs <- attr(pred, "probabilities")
  if (is.null(probs)) stop("SVM probability output is missing")
  if ("TP" %in% colnames(probs)) return(as.numeric(probs[, "TP"]))
  # e1071 may retain internal class order; choose the column labelled TP only.
  stop("TP probability column is missing from the SVM output")
}

svm_recursive_rank <- function(x, y, ranking_cost = 1) {
  current <- colnames(x)
  removed <- character()
  scaled <- scale_train_test(x, x)$train
  while (length(current) > 1L) {
    model <- fit_linear_svm(scaled[, current, drop = FALSE], y,
                            cost = ranking_cost, probability = FALSE)
    weights <- drop(t(model$coefs) %*% model$SV)
    names(weights) <- current
    remove_n <- max(1L, floor(length(current) * 0.10))
    remove_n <- min(remove_n, length(current) - 1L)
    remove_genes <- names(sort(abs(weights), decreasing = FALSE))[seq_len(remove_n)]
    removed <- c(removed, remove_genes)
    current <- setdiff(current, remove_genes)
  }
  c(current, rev(removed))
}

inner_svm_rfe <- function(x, y, patient, fold_id,
                          cost_grid = c(0.01, 0.1, 1, 10)) {
  p <- ncol(x)
  sizes <- sort(unique(c(1:5, 10, 15, 20, 30, 40, p)))
  sizes <- sizes[sizes <= p]
  performance <- list()
  row_i <- 1L

  # The feature ranking used for each inner validation fold is estimated only
  # from the training partition of that fold. This avoids allowing validation
  # labels to influence SVM-RFE ranks or the selected feature count.
  inner_folds <- sort(unique(fold_id))
  fold_rankings <- lapply(inner_folds, function(fold) {
    tr <- fold_id != fold
    svm_recursive_rank(x[tr, , drop = FALSE], y[tr], ranking_cost = 1)
  })
  names(fold_rankings) <- as.character(inner_folds)

  for (size in sizes) {
    for (cost in cost_grid) {
      fold_auc <- numeric()
      for (fold in inner_folds) {
        tr <- fold_id != fold
        va <- fold_id == fold
        if (length(unique(y[tr])) < 2L || length(unique(y[va])) < 2L) next
        features <- fold_rankings[[as.character(fold)]][seq_len(size)]
        scaled <- scale_train_test(
          x[tr, features, drop = FALSE],
          x[va, features, drop = FALSE]
        )
        model <- fit_linear_svm(scaled$train, y[tr], cost = cost,
                                probability = TRUE)
        prob <- predict_tp_probability(model, scaled$test)
        roc_obj <- pROC::roc(
          response = y[va], predictor = prob,
          levels = c("NT", "TP"), direction = "<", quiet = TRUE
        )
        fold_auc <- c(fold_auc, as.numeric(pROC::auc(roc_obj)))
      }
      performance[[row_i]] <- data.frame(
        size = size,
        cost = cost,
        mean_auc = mean(fold_auc),
        sd_auc = stats::sd(fold_auc),
        se_auc = stats::sd(fold_auc) / sqrt(length(fold_auc)),
        valid_folds = length(fold_auc)
      )
      row_i <- row_i + 1L
    }
  }
  performance <- do.call(rbind, performance)
  performance <- performance[is.finite(performance$mean_auc), , drop = FALSE]
  performance <- performance[order(-performance$mean_auc, performance$size,
                                   performance$cost), , drop = FALSE]
  best <- performance[1, , drop = FALSE]
  one_se_threshold <- best$mean_auc - best$se_auc
  eligible <- performance[
    is.finite(performance$se_auc) & performance$mean_auc >= one_se_threshold,
    , drop = FALSE
  ]
  selected_model <- eligible[
    order(eligible$size, -eligible$mean_auc, eligible$cost),
    , drop = FALSE
  ][1, , drop = FALSE]
  final_ranking <- svm_recursive_rank(x, y, ranking_cost = 1)
  list(
    ranking = final_ranking,
    selected = final_ranking[seq_len(selected_model$size)],
    best_cost = selected_model$cost,
    best_size = selected_model$size,
    best_mean_auc = best$mean_auc,
    best_se_auc = best$se_auc,
    one_se_threshold = one_se_threshold,
    selected_mean_auc = selected_model$mean_auc,
    performance = performance
  )
}

run_lasso_selection <- function(x, y, fold_id) {
  y_num <- as.integer(y == "TP")
  tab <- table(y)
  obs_weights <- ifelse(y == "TP", length(y) / (2 * tab[["TP"]]),
                        length(y) / (2 * tab[["NT"]]))
  cvfit <- glmnet::cv.glmnet(
    x = x,
    y = y_num,
    family = "binomial",
    alpha = 1,
    standardize = TRUE,
    weights = obs_weights,
    foldid = fold_id,
    type.measure = "auc",
    keep = TRUE
  )
  coef_mat <- as.matrix(stats::coef(cvfit, s = "lambda.1se"))
  genes <- setdiff(rownames(coef_mat)[coef_mat[, 1] != 0], "(Intercept)")
  list(cvfit = cvfit, genes = genes)
}

candidate_genes_ml <- intersect(bps_targets, rownames(v$E))
if (length(candidate_genes_ml) < 2L) {
  stop("Fewer than two fixed BPS candidate targets are available for ML")
}

ml_x <- t(v$E[candidate_genes_ml, , drop = FALSE])
ml_y <- factor(group, levels = c("NT", "TP"))
ml_patient <- as.character(sample_metadata$patient_id)
outer_fold_id <- make_balanced_group_folds(
  ml_patient, ml_y, k = 5L, seed = 20260901
)

oof_lasso_probability <- rep(NA_real_, nrow(ml_x))
oof_svm_probability <- rep(NA_real_, nrow(ml_x))
outer_selection <- list()
outer_svm_performance <- list()

message_stamp("Starting patient-grouped nested 5-fold ML evaluation")
for (outer_fold in sort(unique(outer_fold_id))) {
  test_idx <- outer_fold_id == outer_fold
  train_idx <- !test_idx
  x_train <- ml_x[train_idx, , drop = FALSE]
  x_test <- ml_x[test_idx, , drop = FALSE]
  y_train <- ml_y[train_idx]
  y_test <- ml_y[test_idx]
  patient_train <- ml_patient[train_idx]

  inner_fold_id <- make_balanced_group_folds(
    patient_train, y_train, k = 5L, seed = 20260901 + outer_fold
  )
  lasso_fit <- run_lasso_selection(x_train, y_train, inner_fold_id)
  svm_fit <- inner_svm_rfe(
    x_train, y_train, patient_train, inner_fold_id,
    cost_grid = c(0.01, 0.1, 1, 10)
  )
  strict_intersection <- intersect(lasso_fit$genes, svm_fit$selected)

  # Each algorithm is evaluated separately in the untouched outer fold. The
  # LASSO-SVM intersection is used only for candidate prioritization and is not
  # replaced by a union when it is empty.
  oof_lasso_probability[test_idx] <- as.numeric(stats::predict(
    lasso_fit$cvfit,
    newx = x_test,
    s = "lambda.1se",
    type = "response"
  ))
  svm_scaled <- scale_train_test(
    x_train[, svm_fit$selected, drop = FALSE],
    x_test[, svm_fit$selected, drop = FALSE]
  )
  svm_outer_model <- fit_linear_svm(
    svm_scaled$train, y_train,
    cost = svm_fit$best_cost, probability = TRUE
  )
  oof_svm_probability[test_idx] <- predict_tp_probability(
    svm_outer_model, svm_scaled$test
  )

  outer_selection[[outer_fold]] <- data.frame(
    outer_fold = outer_fold,
    lasso_n = length(lasso_fit$genes),
    svm_rfe_n = length(svm_fit$selected),
    intersection_n = length(strict_intersection),
    lasso_genes = paste(lasso_fit$genes, collapse = ";"),
    svm_genes = paste(svm_fit$selected, collapse = ";"),
    intersection_genes = paste(strict_intersection, collapse = ";"),
    svm_cost = svm_fit$best_cost,
    svm_selected_size = svm_fit$best_size,
    svm_one_se_threshold = svm_fit$one_se_threshold,
    train_n = sum(train_idx),
    test_n = sum(test_idx),
    train_NT = sum(y_train == "NT"),
    train_TP = sum(y_train == "TP"),
    test_NT = sum(y_test == "NT"),
    test_TP = sum(y_test == "TP"),
    stringsAsFactors = FALSE
  )
  svm_perf_tmp <- svm_fit$performance
  svm_perf_tmp$outer_fold <- outer_fold
  outer_svm_performance[[outer_fold]] <- svm_perf_tmp
  message_stamp(
    "Completed outer fold ", outer_fold,
    "; LASSO=", length(lasso_fit$genes),
    "; SVM-RFE=", length(svm_fit$selected),
    "; strict intersection=", length(strict_intersection)
  )
}

if (anyNA(oof_lasso_probability) || anyNA(oof_svm_probability)) {
  stop("Nested CV produced missing out-of-fold probabilities")
}
outer_selection_table <- do.call(rbind, outer_selection)
outer_svm_performance_table <- do.call(rbind, outer_svm_performance)

binom_metric <- function(success, total) {
  ci <- stats::binom.test(success, total)$conf.int
  c(estimate = success / total, lower = ci[[1]], upper = ci[[2]])
}

binary_metric_vector <- function(outcome, probability, threshold = 0.5) {
  predicted <- factor(
    ifelse(probability >= threshold, "TP", "NT"),
    levels = c("NT", "TP")
  )
  sensitivity <- mean(predicted[outcome == "TP"] == "TP")
  specificity <- mean(predicted[outcome == "NT"] == "NT")
  c(
    AUC = as.numeric(pROC::auc(pROC::roc(
      response = outcome, predictor = probability,
      levels = c("NT", "TP"), direction = "<", quiet = TRUE
    ))),
    Accuracy = mean(predicted == outcome),
    Sensitivity = sensitivity,
    Specificity = specificity,
    `Balanced accuracy` = mean(c(sensitivity, specificity)),
    `Brier score` = mean((as.integer(outcome == "TP") - probability)^2)
  )
}

cluster_bootstrap_metric_ci <- function(
    outcome, probability, cluster, threshold = 0.5,
    repetitions = 1000L, seed = 1L) {
  cluster <- as.character(cluster)
  cluster_rows <- split(seq_along(cluster), cluster)
  cluster_ids <- names(cluster_rows)
  set.seed(seed)
  boot <- matrix(
    NA_real_, nrow = repetitions, ncol = 6L,
    dimnames = list(NULL, names(binary_metric_vector(outcome, probability, threshold)))
  )
  for (b in seq_len(repetitions)) {
    sampled_ids <- sample(cluster_ids, length(cluster_ids), replace = TRUE)
    idx <- unlist(cluster_rows[sampled_ids], use.names = FALSE)
    if (length(unique(outcome[idx])) < 2L) next
    boot[b, ] <- binary_metric_vector(outcome[idx], probability[idx], threshold)
  }
  t(apply(boot, 2, stats::quantile, probs = c(0.025, 0.975),
          na.rm = TRUE, names = FALSE))
}

evaluate_oof_predictions <- function(
    outcome, probability, patient, algorithm, bootstrap_seed) {
  roc_object <- pROC::roc(
    response = outcome,
    predictor = probability,
    levels = c("NT", "TP"),
    direction = "<",
    quiet = TRUE
  )
  point_metrics <- binary_metric_vector(outcome, probability, threshold = 0.5)
  metric_ci <- cluster_bootstrap_metric_ci(
    outcome, probability, patient,
    threshold = 0.5, repetitions = 1000L, seed = bootstrap_seed
  )
  auc_ci_local <- c(
    metric_ci["AUC", 1], point_metrics[["AUC"]], metric_ci["AUC", 2]
  )
  predicted <- factor(
    ifelse(probability >= 0.5, "TP", "NT"),
    levels = c("NT", "TP")
  )
  brier <- point_metrics[["Brier score"]]
  clipped <- pmin(pmax(probability, 1e-6), 1 - 1e-6)
  calibration <- suppressWarnings(stats::glm(
    I(outcome == "TP") ~ stats::qlogis(clipped),
    family = stats::binomial()
  ))
  calibration_intercept_local <- unname(stats::coef(calibration)[[1]])
  calibration_slope_local <- unname(stats::coef(calibration)[[2]])
  metrics <- data.frame(
    algorithm = algorithm,
    metric = c(
      "AUC", "Accuracy", "Sensitivity", "Specificity", "Balanced accuracy",
      "Brier score", "Calibration intercept", "Calibration slope"
    ),
    estimate = c(
      point_metrics[["AUC"]], point_metrics[["Accuracy"]],
      point_metrics[["Sensitivity"]], point_metrics[["Specificity"]],
      point_metrics[["Balanced accuracy"]],
      brier, calibration_intercept_local, calibration_slope_local
    ),
    lower_95_CI = c(
      metric_ci["AUC", 1], metric_ci["Accuracy", 1],
      metric_ci["Sensitivity", 1], metric_ci["Specificity", 1],
      metric_ci["Balanced accuracy", 1], metric_ci["Brier score", 1],
      NA, NA
    ),
    upper_95_CI = c(
      metric_ci["AUC", 2], metric_ci["Accuracy", 2],
      metric_ci["Sensitivity", 2], metric_ci["Specificity", 2],
      metric_ci["Balanced accuracy", 2], metric_ci["Brier score", 2],
      NA, NA
    ),
    threshold = c(NA, rep(0.5, 4), rep(NA, 3)),
    stringsAsFactors = FALSE
  )
  confusion <- as.data.frame(table(
    algorithm = rep(algorithm, length(outcome)),
    observed = outcome,
    predicted = predicted
  ))
  list(
    roc = roc_object,
    auc_ci = auc_ci_local,
    predicted = predicted,
    metrics = metrics,
    confusion = confusion,
    brier = brier,
    calibration_intercept = calibration_intercept_local,
    calibration_slope = calibration_slope_local
  )
}

lasso_oof_evaluation <- evaluate_oof_predictions(
  ml_y, oof_lasso_probability, ml_patient, "LASSO", 20260921
)
svm_oof_evaluation <- evaluate_oof_predictions(
  ml_y, oof_svm_probability, ml_patient, "SVM-RFE", 20260922
)
ml_metrics <- rbind(
  lasso_oof_evaluation$metrics,
  svm_oof_evaluation$metrics
)
confusion_matrix <- rbind(
  lasso_oof_evaluation$confusion,
  svm_oof_evaluation$confusion
)

oof_predictions <- data.frame(
  profile_id = rownames(ml_x),
  patient_id = ml_patient,
  group = ml_y,
  outer_fold = outer_fold_id,
  LASSO_probability_TP = oof_lasso_probability,
  LASSO_predicted_class = lasso_oof_evaluation$predicted,
  SVM_RFE_probability_TP = oof_svm_probability,
  SVM_RFE_predicted_class = svm_oof_evaluation$predicted
)
write_csv_utf8(outer_selection_table, "15_ML_nestedCV_outer_fold_selections.csv")
write_csv_utf8(outer_svm_performance_table, "16_ML_nestedCV_SVM_RFE_inner_performance.csv")
write_csv_utf8(oof_predictions, "17_ML_nestedCV_out_of_fold_predictions.csv")
write_csv_utf8(ml_metrics, "18_ML_nestedCV_performance_metrics.csv")
write_csv_utf8(confusion_matrix, "19_ML_nestedCV_confusion_matrix.csv")

# Final feature selection on all TCGA profiles. This is for interpretation and
# future external application; its apparent in-sample performance is not used.
# The final mechanism-candidate set additionally requires differential
# expression in the corrected full-cohort analysis. That full-cohort DEG filter
# is never used to define the nested-CV feature space or estimate performance.
full_inner_fold_id <- make_balanced_group_folds(
  ml_patient, ml_y, k = 5L, seed = 20260911
)
lasso_full <- run_lasso_selection(ml_x, ml_y, full_inner_fold_id)
svm_rfe_full <- inner_svm_rfe(
  ml_x, ml_y, ml_patient, full_inner_fold_id,
  cost_grid = c(0.01, 0.1, 1, 10)
)
lasso_genes_full <- lasso_full$genes
svm_genes_full <- svm_rfe_full$selected
ml_intersection_full <- intersect(lasso_genes_full, svm_genes_full)
hub_genes <- intersect(ml_intersection_full, bps_overlap$gene)
hub_genes_primary_retained_high_confidence <- intersect(
  hub_genes, bps_targets_high_confidence
)

# Rerun the same full-cohort feature-selection definitions using the stricter
# target universe. This is a target-definition sensitivity analysis, not an
# additional algorithm and not an external model-performance estimate.
candidate_genes_ml_high_confidence <- intersect(
  bps_targets_high_confidence, rownames(v$E)
)
lasso_high_confidence <- NULL
svm_rfe_high_confidence <- NULL
lasso_genes_high_confidence <- character()
svm_genes_high_confidence <- character()
ml_intersection_high_confidence <- character()
hub_genes_high_confidence_sensitivity <- character()
high_confidence_feature_table <- data.frame()
if (length(candidate_genes_ml_high_confidence) >= 2L) {
  ml_x_high_confidence <- t(v$E[
    candidate_genes_ml_high_confidence, , drop = FALSE
  ])
  lasso_high_confidence <- run_lasso_selection(
    ml_x_high_confidence, ml_y, full_inner_fold_id
  )
  svm_rfe_high_confidence <- inner_svm_rfe(
    ml_x_high_confidence, ml_y, ml_patient, full_inner_fold_id,
    cost_grid = c(0.01, 0.1, 1, 10)
  )
  lasso_genes_high_confidence <- lasso_high_confidence$genes
  svm_genes_high_confidence <- svm_rfe_high_confidence$selected
  ml_intersection_high_confidence <- intersect(
    lasso_genes_high_confidence, svm_genes_high_confidence
  )
  hub_genes_high_confidence_sensitivity <- intersect(
    ml_intersection_high_confidence, bps_overlap_high_confidence$gene
  )
  high_confidence_union <- union(
    lasso_genes_high_confidence, svm_genes_high_confidence
  )
  high_confidence_feature_table <- data.frame(
    gene = high_confidence_union,
    selected_by_LASSO_lambda_1se =
      high_confidence_union %in% lasso_genes_high_confidence,
    selected_by_SVM_RFE_one_SE =
      high_confidence_union %in% svm_genes_high_confidence,
    high_confidence_full_cohort_DEG =
      high_confidence_union %in% bps_overlap_high_confidence$gene,
    final_intersection_hub =
      high_confidence_union %in% hub_genes_high_confidence_sensitivity,
    stringsAsFactors = FALSE
  )
}
final_feature_table <- data.frame(
  gene = union(lasso_genes_full, svm_genes_full),
  selected_by_LASSO_lambda_1se = union(lasso_genes_full, svm_genes_full) %in% lasso_genes_full,
  selected_by_SVM_RFE = union(lasso_genes_full, svm_genes_full) %in% svm_genes_full,
  corrected_full_cohort_DEG = union(lasso_genes_full, svm_genes_full) %in% bps_overlap$gene,
  final_intersection_hub = union(lasso_genes_full, svm_genes_full) %in% hub_genes
)
write_csv_utf8(final_feature_table, "20_ML_final_full_TCGA_feature_selection.csv")
hub_target_evidence <- target_source_records[
  target_source_records$gene %in% hub_genes, , drop = FALSE
]
write_csv_utf8(hub_target_evidence, "20B_final_hub_target_evidence.csv")
write_csv_utf8(
  data.frame(
    gene = hub_genes,
    retained_in_high_confidence_target_sensitivity =
      hub_genes %in% bps_targets_high_confidence,
    stringsAsFactors = FALSE
  ),
  "20C_final_hub_high_confidence_target_sensitivity.csv"
)
write_csv_utf8(
  high_confidence_feature_table,
  "20D_high_confidence_target_sensitivity_feature_selection.csv"
)

selection_frequency_for <- function(gene_strings, method) {
  split_genes <- strsplit(gene_strings, ";", fixed = TRUE)
  genes <- sort(unique(unlist(split_genes)))
  genes <- genes[nzchar(genes)]
  if (length(genes) == 0L) {
    return(data.frame(
      method = character(), gene = character(), selected_outer_folds = integer(),
      total_outer_folds = integer(), frequency = numeric()
    ))
  }
  selected_counts <- vapply(genes, function(gene) {
    sum(vapply(split_genes, function(z) gene %in% z, logical(1)))
  }, integer(1))
  data.frame(
    method = method,
    gene = genes,
    selected_outer_folds = selected_counts,
    total_outer_folds = length(split_genes),
    frequency = selected_counts / length(split_genes),
    stringsAsFactors = FALSE
  )
}
selection_frequency <- rbind(
  selection_frequency_for(outer_selection_table$lasso_genes, "LASSO"),
  selection_frequency_for(outer_selection_table$svm_genes, "SVM-RFE"),
  selection_frequency_for(
    outer_selection_table$intersection_genes,
    "Strict LASSO-SVM intersection"
  )
)
selection_frequency <- selection_frequency[
  order(selection_frequency$method, -selection_frequency$frequency,
        selection_frequency$gene),
]
write_csv_utf8(selection_frequency, "21_ML_outer_fold_selection_frequency.csv")

roc_df <- rbind(
  data.frame(
    false_positive_rate = 1 - lasso_oof_evaluation$roc$specificities,
    true_positive_rate = lasso_oof_evaluation$roc$sensitivities,
    algorithm = sprintf(
      "LASSO, AUC %.3f (95%% CI %.3f-%.3f)",
      as.numeric(pROC::auc(lasso_oof_evaluation$roc)),
      lasso_oof_evaluation$auc_ci[[1]], lasso_oof_evaluation$auc_ci[[3]]
    )
  ),
  data.frame(
    false_positive_rate = 1 - svm_oof_evaluation$roc$specificities,
    true_positive_rate = svm_oof_evaluation$roc$sensitivities,
    algorithm = sprintf(
      "SVM-RFE, AUC %.3f (95%% CI %.3f-%.3f)",
      as.numeric(pROC::auc(svm_oof_evaluation$roc)),
      svm_oof_evaluation$auc_ci[[1]], svm_oof_evaluation$auc_ci[[3]]
    )
  )
)
p_roc <- ggplot(
  roc_df,
  aes(false_positive_rate, true_positive_rate, colour = algorithm)
) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#888888") +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = c(palette_contract[["NT"]], palette_contract[["up"]])) +
  coord_equal() +
  labs(
    title = "Nested-CV ROC",
    subtitle = sprintf(
      "OOF AUC: LASSO %.3f; SVM-RFE %.3f",
      as.numeric(pROC::auc(lasso_oof_evaluation$roc)),
      as.numeric(pROC::auc(svm_oof_evaluation$roc))
    ),
    x = "False-positive rate", y = "True-positive rate", colour = NULL
  ) +
  theme(legend.position = "none")

calibration_df <- rbind(
  data.frame(
    probability = oof_lasso_probability,
    observed = as.integer(ml_y == "TP"),
    algorithm = "LASSO"
  ),
  data.frame(
    probability = oof_svm_probability,
    observed = as.integer(ml_y == "TP"),
    algorithm = "SVM-RFE"
  )
)
calibration_summary <- do.call(rbind, lapply(
  split(calibration_df, calibration_df$algorithm),
  function(z) {
    z$bin <- cut(
      z$probability,
      breaks = unique(stats::quantile(
        z$probability, probs = seq(0, 1, 0.1), na.rm = TRUE
      )),
      include.lowest = TRUE
    )
    out <- aggregate(cbind(probability, observed) ~ bin, data = z, FUN = mean)
    out$algorithm <- z$algorithm[[1]]
    out
  }
))
p_calibration <- ggplot(
  calibration_summary,
  aes(probability, observed, colour = algorithm)
) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#888888") +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.7) +
  scale_colour_manual(values = c(palette_contract[["NT"]], palette_contract[["up"]])) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "OOF calibration",
    x = "Mean predicted probability", y = "Observed tumor fraction",
    colour = NULL
  ) +
  theme(legend.position = "bottom")

strict_frequency_plot <- selection_frequency[
  selection_frequency$method == "Strict LASSO-SVM intersection",
]
p_frequency <- ggplot(
  head(strict_frequency_plot, 20L),
  aes(reorder(gene, frequency), frequency)
) +
  geom_col(fill = palette_contract[["accent"]], width = 0.72) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(title = "Strict-intersection stability",
       subtitle = "Patient-grouped outer folds (n = 5)",
       x = NULL, y = "Selection frequency")

lasso_cv_df <- data.frame(
  log_lambda = log(lasso_full$cvfit$lambda),
  mean_auc = lasso_full$cvfit$cvm,
  se_auc = lasso_full$cvfit$cvsd
)
p_lasso <- ggplot(lasso_cv_df, aes(log_lambda, mean_auc)) +
  geom_ribbon(aes(ymin = mean_auc - se_auc, ymax = mean_auc + se_auc),
              fill = "#D9D9D9", alpha = 0.7) +
  geom_line(colour = palette_contract[["NT"]], linewidth = 0.7) +
  geom_vline(xintercept = log(lasso_full$cvfit$lambda.1se),
             linetype = 2, colour = palette_contract[["up"]]) +
  labs(title = "Full-cohort LASSO inner CV",
       subtitle = "Dashed line: lambda.1se",
       x = "log(lambda)", y = "Cross-validated AUC")

svm_cv_df <- svm_rfe_full$performance[
  svm_rfe_full$performance$cost == svm_rfe_full$best_cost,
]
p_svm <- ggplot(svm_cv_df, aes(size, mean_auc)) +
  geom_ribbon(
    aes(ymin = pmax(0, mean_auc - se_auc), ymax = pmin(1, mean_auc + se_auc)),
    fill = "#D9D9D9", alpha = 0.7
  ) +
  geom_line(colour = palette_contract[["teal"]], linewidth = 0.7) +
  geom_point(colour = palette_contract[["teal"]], size = 1.2) +
  geom_hline(yintercept = svm_rfe_full$one_se_threshold, linetype = 3,
             colour = palette_contract[["neutral"]]) +
  geom_vline(xintercept = svm_rfe_full$best_size, linetype = 2,
             colour = palette_contract[["up"]]) +
  scale_x_log10(breaks = unique(svm_cv_df$size)) +
  labs(
    title = "Full-cohort SVM-RFE inner CV",
    subtitle = sprintf("One-SE subset: %d genes; cost %.2g",
                       svm_rfe_full$best_size, svm_rfe_full$best_cost),
    x = "Number of features (log scale)", y = "Cross-validated AUC"
  )

ml_figure <- (p_roc | p_calibration) / (p_lasso | p_svm) +
  patchwork::plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 11))
save_gg(ml_figure, "Figure_ML_nestedCV_performance", 183, 145)
save_gg(p_frequency, "Figure_ML_strict_intersection_stability", 89, 75)

message_stamp(
  "Nested-CV LASSO AUC=",
  sprintf("%.3f", as.numeric(pROC::auc(lasso_oof_evaluation$roc))),
  "; nested-CV SVM-RFE AUC=",
  sprintf("%.3f", as.numeric(pROC::auc(svm_oof_evaluation$roc))),
  "; final LASSO genes=", length(lasso_genes_full),
  "; final SVM-RFE genes=", length(svm_genes_full),
  "; final intersection hubs=", length(hub_genes)
)

# ---- TCGA full-cohort and matched-pair hub validation -------------------------
tcga_hub_results <- deg_all[
  match(hub_genes, deg_all$gene),
  , drop = FALSE
]
tcga_hub_results <- tcga_hub_results[!is.na(tcga_hub_results$gene), ]
write_csv_utf8(tcga_hub_results, "21A_TCGA_full_cohort_hub_validation.csv")

p_tcga_full <- NULL
p_tcga_paired <- NULL
tcga_pair_results <- data.frame()
tcga_pair_metadata <- data.frame()
tcga_pair_v <- NULL

if (length(hub_genes) > 0L) {
  tcga_long <- do.call(rbind, lapply(hub_genes, function(gene) {
    data.frame(
      gene = gene,
      expression = as.numeric(v$E[gene, ]),
      group = group,
      patient_id = sample_metadata$patient_id,
      stringsAsFactors = FALSE
    )
  }))
  tcga_full_annotation <- make_p_annotation(
    tcga_hub_results, "adj.P.Val", "Patient-blocked limma BH P"
  )
  p_tcga_full <- ggplot(tcga_long, aes(group, expression, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.55, colour = NA) +
    geom_boxplot(width = 0.16, outlier.shape = NA, linewidth = 0.35) +
    geom_point(
      position = position_jitter(
        width = 0.10, height = 0, seed = 20260961
      ),
      size = 0.45, alpha = 0.25
    ) +
    geom_text(
      data = tcga_full_annotation,
      aes(x = x, y = y, label = label), inherit.aes = FALSE,
      vjust = 1.15, size = 2.35, lineheight = 0.92
    ) +
    facet_wrap(~ gene, scales = "free_y") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
    scale_fill_manual(values = palette_contract[c("NT", "TP")]) +
    labs(
      title = "TCGA-KIRC: full cohort",
      subtitle = sprintf("%d normal; %d tumor (patient-blocked)",
                         sum(group == "NT"), sum(group == "TP")),
      x = NULL, y = "Voom log2-CPM", fill = "Group"
    ) +
    theme(legend.position = "none")
  save_gg(p_tcga_full, "Figure_TCGA_hub_expression_full_cohort", 89, 78)

  patient_group_profile <- split(
    seq_len(nrow(sample_metadata)),
    sample_metadata$patient_id
  )
  paired_patients <- names(Filter(function(idx) {
    all(c("NT", "TP") %in% as.character(sample_metadata$group[idx]))
  }, patient_group_profile))

  paired_rows <- do.call(rbind, lapply(paired_patients, function(id) {
    idx <- patient_group_profile[[id]]
    data.frame(
      patient_id = id,
      group = c("NT", "TP"),
      profile_id = c(
        sort(sample_metadata$profile_id[idx][sample_metadata$group[idx] == "NT"])[1],
        sort(sample_metadata$profile_id[idx][sample_metadata$group[idx] == "TP"])[1]
      ),
      stringsAsFactors = FALSE
    )
  }))
  tcga_pair_metadata <- paired_rows
  rownames(tcga_pair_metadata) <- tcga_pair_metadata$profile_id
  pair_counts <- counts_symbol[, tcga_pair_metadata$profile_id, drop = FALSE]
  pair_group <- factor(tcga_pair_metadata$group, levels = c("NT", "TP"))
  pair_id <- factor(tcga_pair_metadata$patient_id)
  pair_design <- model.matrix(~ pair_id + pair_group)
  pair_dge <- edgeR::DGEList(pair_counts, group = pair_group)
  pair_keep <- edgeR::filterByExpr(pair_dge, design = pair_design)
  pair_dge <- pair_dge[pair_keep, , keep.lib.sizes = FALSE]
  pair_dge <- edgeR::calcNormFactors(pair_dge, method = "TMM")
  tcga_pair_v <- limma::voom(pair_dge, pair_design, plot = FALSE)
  pair_fit <- limma::lmFit(tcga_pair_v, pair_design)
  pair_fit <- limma::eBayes(pair_fit, robust = TRUE)
  pair_deg <- limma::topTable(
    pair_fit,
    coef = "pair_groupTP",
    number = Inf,
    adjust.method = "BH",
    sort.by = "none",
    confint = 0.95
  )
  pair_deg$gene <- rownames(pair_deg)
  tcga_pair_results <- pair_deg[
    pair_deg$gene %in% hub_genes, , drop = FALSE
  ]

  tcga_pair_wilcox <- do.call(rbind, lapply(hub_genes, function(gene) {
    pair_expr <- data.frame(
      patient_id = tcga_pair_metadata$patient_id,
      group = tcga_pair_metadata$group,
      expression = as.numeric(tcga_pair_v$E[gene, tcga_pair_metadata$profile_id])
    )
    tumor <- pair_expr$expression[pair_expr$group == "TP"]
    normal <- pair_expr$expression[pair_expr$group == "NT"]
    wt <- suppressWarnings(stats::wilcox.test(
      tumor, normal, paired = TRUE, exact = FALSE, conf.int = TRUE
    ))
    data.frame(
      gene = gene,
      matched_pairs = length(paired_patients),
      median_paired_difference = stats::median(tumor - normal),
      paired_wilcoxon_p = wt$p.value,
      pseudomedian_difference = if (!is.null(wt$estimate)) unname(wt$estimate) else NA_real_,
      pseudomedian_CI_lower = if (!is.null(wt$conf.int)) wt$conf.int[[1]] else NA_real_,
      pseudomedian_CI_upper = if (!is.null(wt$conf.int)) wt$conf.int[[2]] else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  tcga_pair_wilcox$paired_wilcoxon_BH <- p.adjust(
    tcga_pair_wilcox$paired_wilcoxon_p, method = "BH"
  )
  tcga_pair_results <- merge(
    tcga_pair_results, tcga_pair_wilcox, by = "gene", all = TRUE
  )
  write_csv_utf8(tcga_pair_metadata, "21B_TCGA_matched_pair_metadata.csv")
  write_csv_utf8(tcga_pair_results, "21C_TCGA_matched_pair_hub_sensitivity.csv")

  tcga_pair_long <- do.call(rbind, lapply(hub_genes, function(gene) {
    data.frame(
      gene = gene,
      expression = as.numeric(tcga_pair_v$E[gene, tcga_pair_metadata$profile_id]),
      group = factor(tcga_pair_metadata$group, levels = c("NT", "TP")),
      patient_id = tcga_pair_metadata$patient_id,
      stringsAsFactors = FALSE
    )
  }))
  tcga_pair_annotation <- make_p_annotation(
    tcga_pair_results, "adj.P.Val", "Paired limma BH P",
    "paired_wilcoxon_BH", "Paired Wilcoxon BH P"
  )
  p_tcga_paired <- ggplot(
    tcga_pair_long,
    aes(group, expression, group = patient_id)
  ) +
    geom_line(colour = "#BDBDBD", linewidth = 0.25, alpha = 0.45) +
    geom_point(aes(colour = group), size = 0.9, alpha = 0.75) +
    geom_text(
      data = tcga_pair_annotation,
      aes(x = x, y = y, label = label), inherit.aes = FALSE,
      vjust = 1.15, size = 2.35, lineheight = 0.92
    ) +
    facet_wrap(~ gene, scales = "free_y") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    scale_colour_manual(values = palette_contract[c("NT", "TP")]) +
    labs(
      title = "TCGA-KIRC: matched pairs",
      subtitle = sprintf("%d tumor-normal pairs", length(paired_patients)),
      x = NULL, y = "Voom log2-CPM", colour = "Group"
    ) +
    theme(legend.position = "top")
  save_gg(p_tcga_paired, "Figure_TCGA_hub_expression_matched_pairs", 89, 78)
}

core_checkpoint_objects <- c(
  "sample_metadata_all", "sample_metadata", "sample_summary",
  "gene_annotation_ensembl", "counts_symbol", "group", "patient_id",
  "design", "dge_unfiltered", "keep_gene", "dge", "v", "tmm_logCPM",
  "duplicate_correlation", "consensus_correlation", "deg_all", "deg_sig",
  "pca", "pca_scores", "batch_association", "batch_adjustment_applied",
  "batch_adjustment_reason", "target_source_records", "target_database_counts",
  "target_legacy_comparison", "target_cross_database_frequency",
  "bps_targets_input", "bps_targets", "bps_targets_high_confidence",
  "bps_overlap", "bps_overlap_high_confidence", "overlap_stability",
  "go_result", "go_table", "kegg_result", "kegg_table", "kegg_snapshot",
  "kegg_snapshot_sha256", "kegg_snapshot_validation",
  "candidate_genes_ml", "ml_x", "ml_y",
  "ml_patient", "outer_fold_id", "outer_selection_table",
  "outer_svm_performance_table", "oof_lasso_probability",
  "oof_svm_probability", "oof_predictions", "lasso_oof_evaluation",
  "svm_oof_evaluation", "ml_metrics", "confusion_matrix",
  "lasso_full", "svm_rfe_full", "lasso_genes_full", "svm_genes_full",
  "ml_intersection_full", "hub_genes",
  "hub_genes_primary_retained_high_confidence",
  "candidate_genes_ml_high_confidence", "ml_x_high_confidence",
  "lasso_high_confidence", "svm_rfe_high_confidence",
  "lasso_genes_high_confidence", "svm_genes_high_confidence",
  "ml_intersection_high_confidence",
  "hub_genes_high_confidence_sensitivity", "high_confidence_feature_table",
  "hub_target_evidence",
  "final_feature_table", "selection_frequency", "tcga_hub_results",
  "tcga_pair_metadata", "tcga_pair_v", "tcga_pair_results",
  "package_versions", "pipeline_version", "analysis_input_signature",
  "analysis_input_signature_hash"
)
core_checkpoint_objects <- core_checkpoint_objects[
  vapply(core_checkpoint_objects, function(object_name) {
    exists(object_name, envir = .GlobalEnv, inherits = FALSE)
  }, logical(1))
]
save(list = core_checkpoint_objects, file = CORE_CHECKPOINT_FILE,
     compress = "gzip")
message_stamp("Saved core checkpoint: ", CORE_CHECKPOINT_FILE)
}

# ---- Paired external-discrimination helper ------------------------------------
paired_expression_auc <- function(
    expression, outcome, pair_id, expected_logFC,
    repetitions = 2000L, seed = 1L) {
  outcome <- factor(outcome, levels = c("NT", "TP"))
  direction_multiplier <- ifelse(expected_logFC >= 0, 1, -1)
  score <- as.numeric(expression) * direction_multiplier
  roc_object <- pROC::roc(
    outcome, score, levels = c("NT", "TP"), direction = "<", quiet = TRUE
  )
  pair_rows <- split(seq_along(pair_id), as.character(pair_id))
  pair_ids <- names(pair_rows)
  set.seed(seed)
  boot_auc <- rep(NA_real_, repetitions)
  for (b in seq_len(repetitions)) {
    sampled_pairs <- sample(pair_ids, length(pair_ids), replace = TRUE)
    idx <- unlist(pair_rows[sampled_pairs], use.names = FALSE)
    if (length(unique(outcome[idx])) < 2L) next
    boot_auc[[b]] <- as.numeric(pROC::auc(pROC::roc(
      outcome[idx], score[idx],
      levels = c("NT", "TP"), direction = "<", quiet = TRUE
    )))
  }
  ci <- stats::quantile(
    boot_auc, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE
  )
  c(
    AUC = as.numeric(pROC::auc(roc_object)),
    AUC_CI_lower = ci[[1]],
    AUC_CI_upper = ci[[2]]
  )
}

# ---- GSE53757 paired external expression validation ---------------------------
geo_validation <- NULL
geo_expression_symbol <- NULL
geo_metadata <- NULL
geo_pair_results <- data.frame()
p_geo <- NULL

if (file.exists(GEO_MATRIX_FILE)) {
  message_stamp("Starting paired GSE53757 validation")
  geo_eset <- GEOquery::getGEO(
    filename = GEO_MATRIX_FILE,
    GSEMatrix = TRUE,
    getGPL = FALSE
  )
  geo_raw <- Biobase::exprs(geo_eset)
  geo_metadata <- Biobase::pData(geo_eset)
  if (max(geo_raw, na.rm = TRUE) > 100) geo_raw <- log2(geo_raw + 1)

  geo_tissue <- if ("tissue:ch1" %in% colnames(geo_metadata)) {
    geo_metadata[["tissue:ch1"]]
  } else {
    geo_metadata[["characteristics_ch1"]]
  }
  geo_group <- ifelse(grepl("normal kidney", geo_tissue, ignore.case = TRUE),
                      "NT", "TP")
  geo_group <- factor(geo_group, levels = c("NT", "TP"))

  # GEO records are 72 consecutive pairs. Some blocks are ordered TP-NT and
  # others NT-TP, so pairing is checked without assuming which member is first.
  pair_members <- split(
    seq_len(ncol(geo_raw)),
    rep(seq_len(ncol(geo_raw) / 2L), each = 2L)
  )
  pair_pattern_ok <- ncol(geo_raw) == 144L && all(vapply(
    pair_members,
    function(idx) identical(sort(as.character(geo_group[idx])), c("NT", "TP")),
    logical(1)
  ))
  if (!pair_pattern_ok) {
    stop("GSE53757 does not contain 72 consecutive one-TP/one-NT pairs")
  }
  geo_pair_id <- factor(rep(seq_len(72L), each = 2L))

  probe_symbol <- AnnotationDbi::mapIds(
    hgu133plus2.db::hgu133plus2.db,
    keys = rownames(geo_raw),
    keytype = "PROBEID",
    column = "SYMBOL",
    multiVals = "first"
  )
  valid_geo <- !is.na(probe_symbol) & nzchar(probe_symbol)
  geo_expression_symbol <- rowsum(
    geo_raw[valid_geo, , drop = FALSE],
    group = unname(probe_symbol[valid_geo]),
    reorder = TRUE
  )
  # Averaging probes after log2 transformation avoids favoring one probe solely
  # because of a higher mean signal.
  probe_counts <- table(unname(probe_symbol[valid_geo]))
  geo_expression_symbol <- sweep(
    geo_expression_symbol,
    1,
    as.numeric(probe_counts[rownames(geo_expression_symbol)]),
    "/"
  )

  geo_metadata$group <- geo_group
  geo_metadata$pair_id <- geo_pair_id
  geo_metadata$sample_order <- seq_len(nrow(geo_metadata))

  validation_genes <- intersect(hub_genes, rownames(geo_expression_symbol))
  if (length(validation_genes) > 0L) {
    geo_design <- model.matrix(~ geo_pair_id + geo_group)
    geo_fit <- limma::lmFit(geo_expression_symbol[validation_genes, , drop = FALSE],
                            geo_design)
    geo_fit <- limma::eBayes(geo_fit, robust = TRUE)
    geo_pair_results <- limma::topTable(
      geo_fit,
      coef = "geo_groupTP",
      number = Inf,
      adjust.method = "BH",
      sort.by = "none",
      confint = 0.95
    )
    # sort.by = "none" preserves the validation_genes input order. Assigning
    # explicitly also handles the single-gene topTable case, where row names can
    # otherwise be replaced by the numeric label "1".
    geo_pair_results$gene <- validation_genes

    wilcox_rows <- lapply(validation_genes, function(gene) {
      tumor <- vapply(pair_members, function(idx) {
        as.numeric(geo_expression_symbol[gene, idx[geo_group[idx] == "TP"]])
      }, numeric(1))
      normal <- vapply(pair_members, function(idx) {
        as.numeric(geo_expression_symbol[gene, idx[geo_group[idx] == "NT"]])
      }, numeric(1))
      wt <- suppressWarnings(stats::wilcox.test(
        tumor, normal, paired = TRUE, exact = FALSE, conf.int = TRUE
      ))
      data.frame(
        gene = gene,
        median_paired_difference = stats::median(tumor - normal),
        paired_wilcoxon_p = wt$p.value,
        pseudomedian_difference = if (!is.null(wt$estimate)) unname(wt$estimate) else NA_real_,
        pseudomedian_CI_lower = if (!is.null(wt$conf.int)) wt$conf.int[[1]] else NA_real_,
        pseudomedian_CI_upper = if (!is.null(wt$conf.int)) wt$conf.int[[2]] else NA_real_
      )
    })
    geo_wilcox <- do.call(rbind, wilcox_rows)
    geo_wilcox$paired_wilcoxon_BH <- p.adjust(geo_wilcox$paired_wilcoxon_p,
                                             method = "BH")
    geo_pair_results <- merge(geo_pair_results, geo_wilcox, by = "gene", all = TRUE)
    geo_auc <- do.call(rbind, lapply(seq_along(validation_genes), function(i) {
      gene <- validation_genes[[i]]
      expected_logFC <- tcga_hub_results$logFC[match(gene, tcga_hub_results$gene)]
      auc_result <- paired_expression_auc(
        expression = geo_expression_symbol[gene, ],
        outcome = geo_group,
        pair_id = geo_pair_id,
        expected_logFC = expected_logFC,
        repetitions = 2000L,
        seed = 20260930 + i
      )
      data.frame(gene = gene, t(auc_result), check.names = FALSE)
    }))
    geo_pair_results <- merge(geo_pair_results, geo_auc, by = "gene", all = TRUE)
    write_csv_utf8(geo_pair_results, "22_GSE53757_paired_hub_validation.csv")
    write_csv_utf8(geo_metadata, "23_GSE53757_sample_metadata_pairing.csv")

    geo_long <- do.call(rbind, lapply(validation_genes, function(gene) {
      data.frame(
        gene = gene,
        expression = as.numeric(geo_expression_symbol[gene, ]),
        group = geo_group,
        pair_id = geo_pair_id
      )
    }))
    geo_annotation <- make_p_annotation(
      geo_pair_results, "adj.P.Val", "Paired limma BH P",
      "paired_wilcoxon_BH", "Paired Wilcoxon BH P"
    )
    p_geo <- ggplot(geo_long, aes(group, expression, group = pair_id)) +
      geom_line(colour = "#BDBDBD", linewidth = 0.25, alpha = 0.45) +
      geom_point(aes(colour = group), size = 1.0, alpha = 0.8) +
      geom_text(
        data = geo_annotation,
        aes(x = x, y = y, label = label), inherit.aes = FALSE,
        vjust = 1.15, size = 2.35, lineheight = 0.92
      ) +
      facet_wrap(~ gene, scales = "free_y", ncol = 3) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
      scale_colour_manual(values = palette_contract[c("NT", "TP")]) +
      labs(title = "GSE53757: matched pairs",
           subtitle = "72 pairs; paired limma and Wilcoxon",
           x = NULL, y = "log2 MAS5 expression", colour = "Group") +
      theme(legend.position = "top")
    save_gg(p_geo, "Figure_GSE53757_paired_validation", 183, 120)
  } else {
    warning("None of the final hub genes mapped to GSE53757")
  }
}

# ---- Optional paired independent RNA-seq validation ---------------------------
independent_validation <- NULL
independent_pair_results <- data.frame()
independent_v <- NULL
p_independent <- NULL

if (file.exists(INDEPENDENT_COUNT_FILE) && length(hub_genes) > 0L) {
  message_stamp("Starting paired independent RNA-seq validation")
  independent_raw <- as.data.frame(readxl::read_excel(INDEPENDENT_COUNT_FILE))
  required_columns <- c(
    "gene_name", paste0("N", 1:10), paste0("T", 1:10)
  )
  if (all(required_columns %in% colnames(independent_raw))) {
    independent_counts_df <- independent_raw[, required_columns]
    independent_counts_df <- independent_counts_df[
      !is.na(independent_counts_df$gene_name) &
        nzchar(trimws(independent_counts_df$gene_name)),
    ]
    independent_count_matrix <- as.matrix(
      independent_counts_df[, setdiff(required_columns, "gene_name")]
    )
    storage.mode(independent_count_matrix) <- "numeric"
    independent_count_matrix <- rowsum(
      independent_count_matrix,
      group = independent_counts_df$gene_name,
      reorder = TRUE
    )
    independent_group <- factor(
      c(rep("NT", 10), rep("TP", 10)),
      levels = c("NT", "TP")
    )
    independent_pair <- factor(rep(1:10, times = 2))
    independent_design <- model.matrix(~ independent_pair + independent_group)
    independent_dge <- edgeR::DGEList(independent_count_matrix,
                                      group = independent_group)
    independent_keep <- edgeR::filterByExpr(
      independent_dge, design = independent_design
    )
    independent_dge <- independent_dge[independent_keep, , keep.lib.sizes = FALSE]
    independent_dge <- edgeR::calcNormFactors(independent_dge, method = "TMM")
    independent_v <- limma::voom(independent_dge, independent_design, plot = FALSE)
    independent_fit <- limma::lmFit(independent_v, independent_design)
    independent_fit <- limma::eBayes(independent_fit, robust = TRUE)
    independent_deg <- limma::topTable(
      independent_fit,
      coef = "independent_groupTP",
      number = Inf,
      adjust.method = "BH",
      sort.by = "none",
      confint = 0.95
    )
    independent_deg$gene <- rownames(independent_deg)
    independent_pair_results <- independent_deg[
      independent_deg$gene %in% hub_genes,
      , drop = FALSE
    ]
    if (nrow(independent_pair_results) > 0L) {
      independent_wilcox <- do.call(rbind, lapply(
        independent_pair_results$gene,
        function(gene) {
          pair_expr <- data.frame(
            pair_id = independent_pair,
            group = independent_group,
            expression = as.numeric(independent_v$E[gene, ])
          )
          pair_levels <- levels(independent_pair)
          tumor <- vapply(pair_levels, function(id) {
            pair_expr$expression[
              pair_expr$pair_id == id & pair_expr$group == "TP"
            ][1]
          }, numeric(1))
          normal <- vapply(pair_levels, function(id) {
            pair_expr$expression[
              pair_expr$pair_id == id & pair_expr$group == "NT"
            ][1]
          }, numeric(1))
          wt <- suppressWarnings(stats::wilcox.test(
            tumor, normal, paired = TRUE, exact = FALSE, conf.int = TRUE
          ))
          data.frame(
            gene = gene,
            matched_pairs = length(pair_levels),
            median_paired_difference = stats::median(tumor - normal),
            paired_wilcoxon_p = wt$p.value,
            pseudomedian_difference = if (!is.null(wt$estimate)) unname(wt$estimate) else NA_real_,
            pseudomedian_CI_lower = if (!is.null(wt$conf.int)) wt$conf.int[[1]] else NA_real_,
            pseudomedian_CI_upper = if (!is.null(wt$conf.int)) wt$conf.int[[2]] else NA_real_,
            stringsAsFactors = FALSE
          )
        }
      ))
      independent_wilcox$paired_wilcoxon_BH <- p.adjust(
        independent_wilcox$paired_wilcoxon_p, method = "BH"
      )
      independent_pair_results <- merge(
        independent_pair_results, independent_wilcox, by = "gene", all = TRUE
      )
      independent_auc <- do.call(rbind, lapply(
        seq_len(nrow(independent_pair_results)),
        function(i) {
          gene <- independent_pair_results$gene[[i]]
          expected_logFC <- tcga_hub_results$logFC[
            match(gene, tcga_hub_results$gene)
          ]
          auc_result <- paired_expression_auc(
            expression = independent_v$E[gene, ],
            outcome = independent_group,
            pair_id = independent_pair,
            expected_logFC = expected_logFC,
            repetitions = 2000L,
            seed = 20260940 + i
          )
          data.frame(gene = gene, t(auc_result), check.names = FALSE)
        }
      ))
      independent_pair_results <- merge(
        independent_pair_results, independent_auc, by = "gene", all = TRUE
      )
    }
    write_csv_utf8(
      independent_pair_results,
      "24_independent_RNAseq_paired_hub_validation.csv"
    )

    independent_genes <- intersect(hub_genes, rownames(independent_v$E))
    if (length(independent_genes) > 0L) {
      independent_long <- do.call(rbind, lapply(independent_genes, function(gene) {
        data.frame(
          gene = gene,
          expression = as.numeric(independent_v$E[gene, ]),
          group = independent_group,
          pair_id = independent_pair
        )
      }))
      independent_annotation <- make_p_annotation(
        independent_pair_results, "adj.P.Val", "Paired limma BH P",
        "paired_wilcoxon_BH", "Paired Wilcoxon BH P"
      )
      p_independent <- ggplot(
        independent_long,
        aes(group, expression, group = pair_id)
      ) +
        geom_line(colour = "#BDBDBD", linewidth = 0.35, alpha = 0.65) +
        geom_point(aes(colour = group), size = 1.6) +
        geom_text(
          data = independent_annotation,
          aes(x = x, y = y, label = label), inherit.aes = FALSE,
          vjust = 1.15, size = 2.35, lineheight = 0.92
        ) +
        facet_wrap(~ gene, scales = "free_y", ncol = 3) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
        scale_colour_manual(values = palette_contract[c("NT", "TP")]) +
        labs(title = "Independent RNA-seq",
             subtitle = "10 matched pairs; paired limma and Wilcoxon",
             x = NULL, y = "Voom log2-CPM", colour = "Group") +
        theme(legend.position = "top")
      save_gg(p_independent, "Figure_independent_RNAseq_paired_validation",
              183, 120)
    }
  } else {
    warning("Independent count file lacks required N1-N10/T1-T10 columns")
  }
}

# ---- Integrated single-gene validation figure --------------------------------
multicohort_validation_summary <- data.frame()
if (length(hub_genes) > 0L) {
  validation_rows <- list()
  if (nrow(tcga_hub_results) > 0L) {
    validation_rows[["TCGA_full"]] <- data.frame(
      dataset = "TCGA-KIRC full cohort",
      gene = tcga_hub_results$gene,
      logFC = tcga_hub_results$logFC,
      CI_lower = tcga_hub_results$CI.L,
      CI_upper = tcga_hub_results$CI.R,
      adjusted_p = tcga_hub_results$adj.P.Val,
      pairs = NA_integer_,
      external_AUC = NA_real_,
      AUC_CI_lower = NA_real_,
      AUC_CI_upper = NA_real_,
      stringsAsFactors = FALSE
    )
  }
  if (nrow(tcga_pair_results) > 0L) {
    validation_rows[["TCGA_pairs"]] <- data.frame(
      dataset = "TCGA-KIRC matched pairs",
      gene = tcga_pair_results$gene,
      logFC = tcga_pair_results$logFC,
      CI_lower = tcga_pair_results$CI.L,
      CI_upper = tcga_pair_results$CI.R,
      adjusted_p = tcga_pair_results$adj.P.Val,
      pairs = tcga_pair_results$matched_pairs,
      external_AUC = NA_real_,
      AUC_CI_lower = NA_real_,
      AUC_CI_upper = NA_real_,
      stringsAsFactors = FALSE
    )
  }
  if (nrow(geo_pair_results) > 0L) {
    validation_rows[["GEO"]] <- data.frame(
      dataset = "GSE53757 matched pairs",
      gene = geo_pair_results$gene,
      logFC = geo_pair_results$logFC,
      CI_lower = geo_pair_results$CI.L,
      CI_upper = geo_pair_results$CI.R,
      adjusted_p = geo_pair_results$adj.P.Val,
      pairs = 72L,
      external_AUC = geo_pair_results$AUC,
      AUC_CI_lower = geo_pair_results$AUC_CI_lower,
      AUC_CI_upper = geo_pair_results$AUC_CI_upper,
      stringsAsFactors = FALSE
    )
  }
  if (nrow(independent_pair_results) > 0L) {
    validation_rows[["independent"]] <- data.frame(
      dataset = "Independent RNA-seq matched pairs",
      gene = independent_pair_results$gene,
      logFC = independent_pair_results$logFC,
      CI_lower = independent_pair_results$CI.L,
      CI_upper = independent_pair_results$CI.R,
      adjusted_p = independent_pair_results$adj.P.Val,
      pairs = 10L,
      external_AUC = independent_pair_results$AUC,
      AUC_CI_lower = independent_pair_results$AUC_CI_lower,
      AUC_CI_upper = independent_pair_results$AUC_CI_upper,
      stringsAsFactors = FALSE
    )
  }
  multicohort_validation_summary <- do.call(rbind, validation_rows)
  rownames(multicohort_validation_summary) <- NULL
  write_csv_utf8(
    multicohort_validation_summary,
    "24B_multicohort_single_hub_validation_summary.csv"
  )

  forest_df <- multicohort_validation_summary
  forest_df$dataset <- factor(
    forest_df$dataset,
    levels = rev(c(
      "TCGA-KIRC full cohort", "TCGA-KIRC matched pairs",
      "GSE53757 matched pairs", "Independent RNA-seq matched pairs"
    ))
  )
  p_forest <- ggplot(forest_df, aes(logFC, dataset)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "#888888") +
    geom_errorbar(
      aes(xmin = CI_lower, xmax = CI_upper),
      orientation = "y", width = 0.16, linewidth = 0.45,
      colour = palette_contract[["neutral"]]
    ) +
    geom_point(size = 2.0, colour = palette_contract[["up"]]) +
    labs(
      title = paste0(hub_genes[[1]], " expression across cohorts"),
      subtitle = "Tumor vs normal; log2FC (95% CI)",
      x = "Log2 fold change", y = NULL
    )
  save_gg(p_forest, "Figure_multicohort_hub_effect_forest", 120, 82)

  # Plot objects are deliberately not required for the numeric core
  # checkpoint. Reconstruct the two TCGA panels when resuming downstream work.
  if (!exists("p_tcga_full", inherits = FALSE) || is.null(p_tcga_full)) {
    tcga_long <- do.call(rbind, lapply(hub_genes, function(gene) {
      data.frame(
        gene = gene,
        expression = as.numeric(v$E[gene, ]),
        group = group,
        patient_id = sample_metadata$patient_id,
        stringsAsFactors = FALSE
      )
    }))
    tcga_full_annotation <- make_p_annotation(
      tcga_hub_results, "adj.P.Val", "Patient-blocked limma BH P"
    )
    p_tcga_full <- ggplot(tcga_long, aes(group, expression, fill = group)) +
      geom_violin(trim = FALSE, alpha = 0.55, colour = NA) +
      geom_boxplot(width = 0.16, outlier.shape = NA, linewidth = 0.35) +
      geom_point(
        position = position_jitter(
          width = 0.10, height = 0, seed = 20260961
        ),
        size = 0.45, alpha = 0.25
      ) +
      geom_text(
        data = tcga_full_annotation,
        aes(x = x, y = y, label = label), inherit.aes = FALSE,
        vjust = 1.15, size = 2.35, lineheight = 0.92
      ) +
      facet_wrap(~ gene, scales = "free_y") +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
      scale_fill_manual(values = palette_contract[c("NT", "TP")]) +
      labs(title = "TCGA-KIRC: full cohort",
           subtitle = sprintf(
             "%d normal; %d tumor (patient-blocked)",
             sum(group == "NT"), sum(group == "TP")
           ), x = NULL, y = "Voom log2-CPM", fill = "Group") +
      theme(legend.position = "none")
  }
  if ((!exists("p_tcga_paired", inherits = FALSE) ||
       is.null(p_tcga_paired)) &&
      !is.null(tcga_pair_v) && nrow(tcga_pair_metadata) > 0L) {
    tcga_pair_long <- do.call(rbind, lapply(hub_genes, function(gene) {
      data.frame(
        gene = gene,
        expression = as.numeric(
          tcga_pair_v$E[gene, tcga_pair_metadata$profile_id]
        ),
        group = factor(tcga_pair_metadata$group, levels = c("NT", "TP")),
        patient_id = tcga_pair_metadata$patient_id,
        stringsAsFactors = FALSE
      )
    }))
    tcga_pair_annotation <- make_p_annotation(
      tcga_pair_results, "adj.P.Val", "Paired limma BH P",
      "paired_wilcoxon_BH", "Paired Wilcoxon BH P"
    )
    p_tcga_paired <- ggplot(
      tcga_pair_long, aes(group, expression, group = patient_id)
    ) +
      geom_line(colour = "#BDBDBD", linewidth = 0.25, alpha = 0.45) +
      geom_point(aes(colour = group), size = 0.9, alpha = 0.75) +
      geom_text(
        data = tcga_pair_annotation,
        aes(x = x, y = y, label = label), inherit.aes = FALSE,
        vjust = 1.15, size = 2.35, lineheight = 0.92
      ) +
      facet_wrap(~ gene, scales = "free_y") +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
      scale_colour_manual(values = palette_contract[c("NT", "TP")]) +
      labs(title = "TCGA-KIRC: matched pairs",
           subtitle = sprintf(
             "%d tumor-normal pairs",
             length(unique(tcga_pair_metadata$patient_id))
           ), x = NULL, y = "Voom log2-CPM", colour = "Group") +
      theme(legend.position = "top")
  }

  # Refresh the standalone TCGA panels as well when this run resumed from the
  # numeric core checkpoint and therefore reconstructed the plot objects here.
  if (!is.null(p_tcga_full)) {
    save_gg(p_tcga_full, "Figure_TCGA_hub_expression_full_cohort", 89, 78)
  }
  if (!is.null(p_tcga_paired)) {
    save_gg(p_tcga_paired, "Figure_TCGA_hub_expression_matched_pairs", 89, 78)
  }

  combined_panel_names <- c(
    "p_tcga_full", "p_tcga_paired", "p_geo", "p_independent"
  )
  combined_panels_available <- all(vapply(
    combined_panel_names,
    function(plot_name) {
      exists(plot_name, envir = .GlobalEnv, inherits = FALSE) &&
        !is.null(get(plot_name, envir = .GlobalEnv, inherits = FALSE))
    },
    logical(1)
  ))
  if (combined_panels_available) {
    multicohort_figure <- (p_tcga_full | p_tcga_paired) /
      (p_geo | p_independent) +
      patchwork::plot_annotation(tag_levels = "a") &
      theme(plot.tag = element_text(face = "bold", size = 11))
    save_gg(
      multicohort_figure,
      "Figure_ESRRB_multicohort_validation",
      183, 145
    )
  }
}

# ---- CIBERSORT sensitivity analysis -------------------------------------------
cibersort_results <- NULL
cibersort_qc <- NULL
immune_group_results <- data.frame()
immune_wilcox_results <- data.frame()
hub_immune_correlations <- data.frame()
immune_consensus_correlation <- NA_real_
cibersort_completed <- FALSE

if (nrow(v$E) > 0L) {
  message_stamp("Starting CIBERSORT with RNA-seq quantile normalization disabled")
  tryCatch({
    data("LM22", package = "CIBERSORT", envir = environment())
    # CIBERSORT expects non-log, non-negative expression. TMM-normalized CPM is
    # therefore used as the RNA-seq mixture matrix; QN remains disabled.
    immune_cpm <- edgeR::cpm(dge, normalized.lib.sizes = TRUE, log = FALSE)
    common_immune_genes <- intersect(rownames(LM22), rownames(immune_cpm))
    immune_input <- immune_cpm[common_immune_genes, , drop = FALSE]
    cibersort_signature <- list(
      signature_matrix = "LM22",
      permutations = 1000L,
      quantile_normalization = FALSE,
      RNG_seed = 20260950L,
      RNG_kind = RNGkind(),
      genes = rownames(immune_input),
      profiles = colnames(immune_input),
      mixture_sha256 = digest::digest(
        immune_input, algo = "sha256", serialize = TRUE
      ),
      LM22_sha256 = digest::digest(
        LM22, algo = "sha256", serialize = TRUE
      ),
      CIBERSORT_version = as.character(utils::packageVersion("CIBERSORT")),
      R_version = R.version.string
    )
    cibersort_cache_valid <- FALSE
    if (file.exists(CIBERSORT_CHECKPOINT_FILE)) {
      cibersort_cache_env <- new.env(parent = emptyenv())
      cibersort_cache_names <- tryCatch(
        load(CIBERSORT_CHECKPOINT_FILE, envir = cibersort_cache_env),
        error = function(e) character()
      )
      if (all(c("cibersort_results", "cibersort_signature") %in%
              cibersort_cache_names) &&
          identical(
            cibersort_cache_env$cibersort_signature,
            cibersort_signature
          )) {
        cibersort_results <- cibersort_cache_env$cibersort_results
        cibersort_cache_valid <- TRUE
        message_stamp(
          "Using validated CIBERSORT checkpoint: ",
          CIBERSORT_CHECKPOINT_FILE
        )
      }
      rm(cibersort_cache_env, cibersort_cache_names)
    }
    if (!cibersort_cache_valid) {
      set.seed(20260950)
      cibersort_results <- CIBERSORT::cibersort(
        sig_matrix = LM22,
        mixture_file = immune_input,
        perm = 1000,
        QN = FALSE
      )
      save(
        cibersort_results, cibersort_signature,
        file = CIBERSORT_CHECKPOINT_FILE,
        compress = "gzip"
      )
      message_stamp("Saved CIBERSORT checkpoint: ", CIBERSORT_CHECKPOINT_FILE)
    }
    cibersort_qc <- data.frame(
      profile_id = rownames(cibersort_results),
      sample_metadata[rownames(cibersort_results), c("patient_id", "group")],
      cibersort_results[, setdiff(colnames(cibersort_results), colnames(LM22)),
                        drop = FALSE],
      check.names = FALSE
    )
    write_csv_utf8(cibersort_results, "25_CIBERSORT_all_profiles.csv",
                   row.names = TRUE)
    write_csv_utf8(cibersort_qc, "26_CIBERSORT_QC_metrics.csv")

    p_col <- grep("P.value|P-value|P.value", colnames(cibersort_results),
                  ignore.case = TRUE, value = TRUE)[1]
    valid_profiles <- if (!is.na(p_col) && length(p_col) > 0L) {
      rownames(cibersort_results)[
        as.numeric(cibersort_results[, p_col]) < 0.05
      ]
    } else rownames(cibersort_results)

    fractions <- as.matrix(cibersort_results[valid_profiles, colnames(LM22),
                                             drop = FALSE])
    immune_estimable_cells <- colnames(fractions)[vapply(
      as.data.frame(fractions),
      function(z) stats::sd(z, na.rm = TRUE) > 0,
      logical(1)
    )]
    immune_non_estimable_cells <- setdiff(
      colnames(fractions), immune_estimable_cells
    )
    if (length(immune_estimable_cells) < 2L) {
      stop("Fewer than two non-constant CIBERSORT cell fractions")
    }
    # Exclude all-zero/constant components before zero replacement and CLR.
    # Otherwise closure can create an apparent CLR difference for a component
    # whose raw estimated fraction is identically zero.
    fractions_clr_input <- pmax(
      fractions[, immune_estimable_cells, drop = FALSE], 1e-6
    )
    clr_fractions <- log(fractions_clr_input) -
      rowMeans(log(fractions_clr_input))
    immune_meta <- sample_metadata[rownames(clr_fractions), , drop = FALSE]
    immune_group <- factor(immune_meta$group, levels = c("NT", "TP"))
    immune_patient <- factor(immune_meta$patient_id)
    immune_design <- model.matrix(~ 0 + immune_group)
    colnames(immune_design) <- c("NT", "TP")
    immune_v <- t(clr_fractions)
    immune_corfit <- tryCatch(
      limma::duplicateCorrelation(immune_v, immune_design, block = immune_patient),
      error = function(e) list(consensus.correlation = NA_real_)
    )
    immune_consensus_correlation <- immune_corfit$consensus.correlation
    if (!is.finite(immune_consensus_correlation)) {
      immune_consensus_correlation <- 0
      warning(
        "CIBERSORT-valid profiles contained insufficient repeated-patient ",
        "information; immune block correlation was set to 0"
      )
    }
    immune_fit <- limma::lmFit(
      immune_v, immune_design, block = immune_patient,
      correlation = immune_consensus_correlation
    )
    immune_fit <- limma::contrasts.fit(
      immune_fit,
      limma::makeContrasts(TP - NT, levels = immune_design)
    )
    immune_fit <- limma::eBayes(immune_fit, robust = TRUE)
    immune_group_results <- limma::topTable(
      immune_fit, number = Inf, adjust.method = "BH", sort.by = "P"
    )
    immune_group_results$cell_type <- rownames(immune_group_results)
    write_csv_utf8(
      immune_group_results,
      "27_CIBERSORT_group_comparison_CLR_patient_blocked.csv"
    )

    immune_long <- as.data.frame(fractions)
    immune_long$profile_id <- rownames(immune_long)
    immune_long$group <- immune_group
    immune_long <- data.table::melt(
      data.table::as.data.table(immune_long),
      id.vars = c("profile_id", "group"),
      variable.name = "cell_type",
      value.name = "fraction"
    )
    immune_long <- as.data.frame(immune_long)
    immune_wilcox_results <- do.call(rbind, lapply(
      colnames(fractions),
      function(cell) {
        tp <- fractions[immune_group == "TP", cell]
        nt <- fractions[immune_group == "NT", cell]
        estimable <- stats::sd(c(tp, nt), na.rm = TRUE) > 0
        wt <- if (estimable) {
          suppressWarnings(stats::wilcox.test(
            tp, nt, paired = FALSE, exact = FALSE, conf.int = FALSE
          ))
        } else NULL
        u_statistic <- if (estimable) unname(wt$statistic) else NA_real_
        data.frame(
          cell_type = cell,
          n_NT = length(nt),
          n_TP = length(tp),
          median_NT = stats::median(nt),
          median_TP = stats::median(tp),
          median_difference_TP_minus_NT = stats::median(tp) - stats::median(nt),
          Wilcoxon_W = u_statistic,
          rank_biserial_TP_vs_NT =
            if (estimable) {
              2 * u_statistic / (length(tp) * length(nt)) - 1
            } else NA_real_,
          p_value = if (estimable) wt$p.value else NA_real_,
          estimable = estimable,
          stringsAsFactors = FALSE
        )
      }
    ))
    immune_wilcox_results$BH_adjusted_p <- p.adjust(
      immune_wilcox_results$p_value, method = "BH"
    )
    immune_wilcox_results <- immune_wilcox_results[
      order(immune_wilcox_results$BH_adjusted_p), , drop = FALSE
    ]
    write_csv_utf8(
      immune_wilcox_results,
      "27B_CIBERSORT_group_comparison_Wilcoxon_BH.csv"
    )

    immune_annotation <- immune_wilcox_results[, c(
      "cell_type", "BH_adjusted_p"
    )]
    immune_annotation$x <- 1.5
    immune_annotation$y <- Inf
    immune_annotation$label <- paste0(
      "BH P = ",
      vapply(immune_annotation$BH_adjusted_p, format_p_number, character(1))
    )
    immune_long$cell_type <- factor(
      immune_long$cell_type, levels = colnames(fractions)
    )
    immune_annotation$cell_type <- factor(
      immune_annotation$cell_type, levels = colnames(fractions)
    )
    p_immune <- ggplot(immune_long, aes(group, fraction, fill = group)) +
      geom_boxplot(outlier.shape = NA, linewidth = 0.3, width = 0.56) +
      geom_point(
        position = position_jitter(
          width = 0.10, height = 0, seed = 20260962
        ),
        size = 0.35, alpha = 0.22
      ) +
      geom_text(
        data = immune_annotation,
        aes(x = x, y = y, label = label), inherit.aes = FALSE,
        vjust = 1.15, size = 2.15
      ) +
      scale_y_continuous(
        expand = expansion(mult = c(0.04, 0.18), add = c(0, 0.002))
      ) +
      facet_wrap(
        ~ cell_type, ncol = 4, scales = "free_y",
        labeller = ggplot2::labeller(
          cell_type = ggplot2::label_wrap_gen(width = 20)
        )
      ) +
      scale_fill_manual(values = palette_contract[c("NT", "TP")]) +
      labs(
        title = "CIBERSORT-estimated relative immune-cell fractions",
        subtitle = sprintf(
          paste0(
            "Valid deconvolution profiles (P < 0.05): NT=%d, TP=%d\n",
            "Two-sided Wilcoxon rank-sum tests; 22 assessed, %d estimable; BH correction across estimable tests"
          ),
          sum(immune_group == "NT"), sum(immune_group == "TP"),
          sum(immune_wilcox_results$estimable)
        ),
        x = NULL, y = "Estimated relative fraction", fill = "Group"
      ) +
      theme(
        legend.position = "none",
        strip.text = element_text(size = 7),
        axis.text = element_text(size = 7),
        plot.subtitle = element_text(size = 8)
      )
    save_gg(p_immune, "Figure_CIBERSORT_group_comparison", 183, 188)

    tumor_profiles <- rownames(immune_meta)[immune_group == "TP"]
    correlation_genes <- intersect(hub_genes, rownames(v$E))
    if (length(correlation_genes) > 0L && length(tumor_profiles) > 20L) {
      cor_rows <- list()
      idx <- 1L
      for (gene in correlation_genes) {
        for (cell in colnames(fractions)) {
          expression_values <- as.numeric(v$E[gene, tumor_profiles])
          fraction_values <- as.numeric(fractions[tumor_profiles, cell])
          complete <- is.finite(expression_values) & is.finite(fraction_values)
          estimable <- sum(complete) >= 3L &&
            stats::sd(fraction_values[complete]) > 0
          ct <- if (estimable) {
            suppressWarnings(stats::cor.test(
              expression_values[complete], fraction_values[complete],
              method = "spearman", exact = FALSE
            ))
          } else NULL
          cor_rows[[idx]] <- data.frame(
            gene = gene, cell_type = cell,
            n = sum(complete),
            rho = if (estimable) unname(ct$estimate) else NA_real_,
            p_value = if (estimable) ct$p.value else NA_real_,
            estimable = estimable
          )
          idx <- idx + 1L
        }
      }
      hub_immune_correlations <- do.call(rbind, cor_rows)
      hub_immune_correlations$BH_adjusted_p <- p.adjust(
        hub_immune_correlations$p_value, method = "BH"
      )
      write_csv_utf8(
        hub_immune_correlations,
        "28_hub_CIBERSORT_correlations_tumor_only_BH.csv"
      )
    }
    cibersort_completed <- TRUE
  }, error = function(e) {
    warning("CIBERSORT module failed but core analysis was retained: ",
            conditionMessage(e))
  })
}

# ---- Analysis summary and saved workspace ------------------------------------
lasso_auc_value <- as.numeric(pROC::auc(lasso_oof_evaluation$roc))
svm_auc_value <- as.numeric(pROC::auc(svm_oof_evaluation$roc))
strict_hub_frequency <- if (length(hub_genes) > 0L) {
  selection_frequency$frequency[
    selection_frequency$method == "Strict LASSO-SVM intersection" &
      selection_frequency$gene == hub_genes[[1]]
  ][1]
} else NA_real_
analysis_summary <- data.frame(
  item = c(
    "TCGA tumor profiles", "TCGA normal profiles", "TCGA unique patients",
    "Genes before filtering", "Genes after filterByExpr",
    "Patient-block consensus correlation", "Corrected DEGs",
    "Legacy author-supplied BPS entries", "Reconstructed BPS candidate targets",
    "High-confidence BPS sensitivity targets", "Corrected BPS-DEG overlap",
    "High-confidence BPS-DEG sensitivity overlap",
    "Nested-CV LASSO AUC", "Nested-CV LASSO AUC lower 95% CI",
    "Nested-CV LASSO AUC upper 95% CI",
    "Nested-CV SVM-RFE AUC", "Nested-CV SVM-RFE AUC lower 95% CI",
    "Nested-CV SVM-RFE AUC upper 95% CI", "Final LASSO features",
    "Final SVM-RFE features", "Final intersection hub genes",
    "Primary final hubs present in high-confidence target set",
    "High-confidence target sensitivity hub genes",
    "Strict outer-fold intersection frequency of final hub",
    "TCGA matched tumor-normal pairs",
    "Batch adjustment applied", "GSE53757 paired validation available",
    "Independent paired RNA-seq validation available",
    "CIBERSORT completed"
  ),
  value = c(
    sum(group == "TP"), sum(group == "NT"), length(unique(patient_id)),
    nrow(dge_unfiltered), nrow(dge), consensus_correlation, nrow(deg_sig),
    length(bps_targets_input), length(bps_targets),
    length(bps_targets_high_confidence), nrow(bps_overlap),
    nrow(bps_overlap_high_confidence),
    lasso_auc_value, lasso_oof_evaluation$auc_ci[[1]],
    lasso_oof_evaluation$auc_ci[[3]],
    svm_auc_value, svm_oof_evaluation$auc_ci[[1]],
    svm_oof_evaluation$auc_ci[[3]],
    length(lasso_genes_full), length(svm_genes_full), length(hub_genes),
    length(hub_genes_primary_retained_high_confidence),
    length(hub_genes_high_confidence_sensitivity), strict_hub_frequency,
    nrow(tcga_pair_metadata) / 2,
    batch_adjustment_applied, nrow(geo_pair_results) > 0,
    nrow(independent_pair_results) > 0, cibersort_completed
  )
)
write_csv_utf8(analysis_summary, "29_analysis_summary.csv")
writeLines(
  c(
    "Corrected TCGA-KIRC/BPS computational reanalysis",
    paste("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste("Pipeline version:", pipeline_version),
    paste("Input signature SHA-256:", analysis_input_signature_hash),
    paste("R version:", R.version.string),
    paste("Tumor profiles:", sum(group == "TP")),
    paste("Normal profiles:", sum(group == "NT")),
    paste("Unique patients:", length(unique(patient_id))),
    paste("Genes retained after filterByExpr:", nrow(dge)),
    paste("Corrected DEGs:", nrow(deg_sig)),
    paste("Reconstructed BPS candidate targets:", length(bps_targets)),
    paste("High-confidence target sensitivity set:",
          length(bps_targets_high_confidence)),
    paste("Corrected BPS-DEG overlap:", nrow(bps_overlap)),
    paste("Nested-CV LASSO AUC:", sprintf("%.3f", lasso_auc_value)),
    paste("Nested-CV SVM-RFE AUC:", sprintf("%.3f", svm_auc_value)),
    paste("Final hub genes:", paste(hub_genes, collapse = ", ")),
    paste("Strict outer-fold intersection frequency:",
          sprintf("%.3f", strict_hub_frequency)),
    paste("Hub retained in high-confidence target sensitivity:",
          paste(hub_genes_primary_retained_high_confidence, collapse = ", ")),
    paste("High-confidence sensitivity hub genes:",
          paste(hub_genes_high_confidence_sensitivity, collapse = ", ")),
    "Interpretation: computational candidate identification and hypothesis generation only."
  ),
  file.path(OUTPUT_DIR, "README_results_summary.txt"),
  useBytes = TRUE
)

session_info <- utils::sessionInfo()
capture.output(session_info, file = file.path(LOG_DIR, "sessionInfo.txt"))

objects_to_save <- c(
  "sample_metadata_all", "sample_metadata", "sample_summary",
  "gene_annotation_ensembl", "counts_symbol", "group", "patient_id",
  "design", "dge_unfiltered", "keep_gene", "dge", "v", "tmm_logCPM",
  "duplicate_correlation",
  "consensus_correlation", "deg_all", "deg_sig", "pca", "pca_scores",
  "batch_association", "batch_adjustment_applied", "batch_adjustment_reason",
  "target_source_records", "target_database_counts", "target_legacy_comparison",
  "target_cross_database_frequency", "bps_targets_input", "bps_targets",
  "bps_targets_high_confidence", "bps_overlap", "bps_overlap_high_confidence",
  "overlap_stability",
  "go_result", "go_table", "kegg_result", "kegg_table", "kegg_snapshot",
  "kegg_snapshot_sha256", "kegg_snapshot_validation",
  "candidate_genes_ml", "ml_x", "ml_y", "outer_fold_id",
  "outer_selection_table", "outer_svm_performance_table",
  "oof_lasso_probability", "oof_svm_probability", "oof_predictions",
  "lasso_oof_evaluation", "svm_oof_evaluation", "ml_metrics",
  "confusion_matrix", "lasso_full", "svm_rfe_full", "lasso_genes_full",
  "svm_genes_full", "ml_intersection_full", "hub_genes",
  "hub_genes_primary_retained_high_confidence",
  "candidate_genes_ml_high_confidence", "ml_x_high_confidence",
  "lasso_high_confidence", "svm_rfe_high_confidence",
  "lasso_genes_high_confidence", "svm_genes_high_confidence",
  "ml_intersection_high_confidence",
  "hub_genes_high_confidence_sensitivity", "high_confidence_feature_table",
  "hub_target_evidence",
  "final_feature_table", "selection_frequency", "tcga_hub_results",
  "tcga_pair_metadata", "tcga_pair_v", "tcga_pair_results",
  "geo_expression_symbol", "geo_metadata",
  "geo_pair_results", "independent_v", "independent_pair_results",
  "multicohort_validation_summary",
  "cibersort_results", "cibersort_qc", "immune_group_results",
  "immune_wilcox_results", "immune_estimable_cells",
  "immune_non_estimable_cells",
  "immune_consensus_correlation", "cibersort_completed",
  "hub_immune_correlations", "analysis_summary", "lasso_auc_value",
  "svm_auc_value", "strict_hub_frequency", "package_versions",
  "session_info", "pipeline_version", "analysis_input_signature",
  "analysis_input_signature_hash"
)
objects_to_save <- objects_to_save[vapply(objects_to_save, function(object_name) {
  exists(object_name, envir = .GlobalEnv, inherits = FALSE)
}, logical(1))]
save(list = objects_to_save, file = RDATA_FILE, compress = "xz")

message_stamp("Saved R workspace: ", RDATA_FILE)
message_stamp("Analysis completed successfully")
