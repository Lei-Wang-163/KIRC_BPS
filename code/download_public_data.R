options(stringsAsFactors = FALSE, timeout = 3600)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]]),
                        winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
repo_dir <- normalizePath(file.path(script_dir, ".."),
                          winslash = "/", mustWork = TRUE)
args <- commandArgs(trailingOnly = TRUE)

if (!any(args %in% c("--geo", "--tcga"))) {
  cat("Usage:\n")
  cat("  Rscript code/download_public_data.R --geo\n")
  cat("  GDC_CLIENT=/path/to/gdc-client Rscript code/download_public_data.R --tcga\n")
  quit(status = 0L)
}

if ("--geo" %in% args) {
  geo_dir <- file.path(repo_dir, "data", "external", "GSE53757")
  dir.create(geo_dir, recursive = TRUE, showWarnings = FALSE)
  geo_file <- file.path(geo_dir, "GSE53757_series_matrix.txt.gz")
  geo_url <- paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE53nnn/GSE53757/matrix/",
    "GSE53757_series_matrix.txt.gz"
  )
  if (!file.exists(geo_file)) {
    message("Downloading GSE53757 series matrix from NCBI GEO")
    utils::download.file(geo_url, geo_file, mode = "wb", quiet = FALSE)
  } else {
    message("GSE53757 file already exists: ", geo_file)
  }
}

if ("--tcga" %in% args) {
  gdc_client <- Sys.getenv("GDC_CLIENT", unset = "")
  if (!nzchar(gdc_client) || !file.exists(gdc_client)) {
    stop("Set GDC_CLIENT to the installed gdc-client executable.")
  }
  manifest <- file.path(
    repo_dir, "data", "input", "manifests", "gdc_manifest_tcga_kirc.txt"
  )
  destination <- file.path(
    repo_dir, "data", "external", "TCGA-KIRC", "Transcriptome_Profiling",
    "Gene_Expression_Quantification"
  )
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  status <- system2(
    gdc_client,
    args = c("download", "-m", shQuote(manifest), "-d", shQuote(destination))
  )
  if (!identical(status, 0L)) stop("gdc-client exited with status ", status)
}

