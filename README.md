# BPS–ccRCC computational analysis code

This minimal repository contains the author-written R code underlying the
revised manuscript:

> Multi-cohort transcriptomic validation combined with machine learning and
> network toxicology reveals potential mechanisms of bisphenol S in ccRCC

The workflow identifies computational candidate targets and generates a
mechanistic hypothesis. It does not establish BPS exposure, direct target
engagement, causality, or clinical utility.

## Repository contents

```text
code/
  KIRC_BPS_pipeline.R
  KIRC_BPS_ESRRB_CIBERSORT_correlation_figure.R
  download_public_data.R
data/input/manifests/
  gdc_manifest_tcga_kirc.txt
renv.lock
README.md
.gitignore
```

No independent-cohort data or results, expression matrices, figure source
data, database exports, docking files, MD files, KEGG content, R workspaces, or
analysis caches are included in this minimal code repository.

## Public datasets

- TCGA-KIRC STAR-count data: NCI Genomic Data Commons, project `TCGA-KIRC`.
  The exact public file manifest used in the analysis is included. It contains
  614 files; the code retains 541 primary-tumour and 72 solid-tissue-normal
  profiles after sample-type filtering.
- External microarray validation: NCBI GEO accession `GSE53757` (72 matched
  tumour–normal pairs).
- Structural model: Protein Data Bank entry `6LN4`.

The independent 10-pair RNA-seq cohort is not part of this GitHub repository.
Its repository/accession and access conditions must be reported separately in
the manuscript's Data Availability statement.

## R environment

The recorded analysis used R 4.3.2 and Bioconductor 3.18. Restore the package
environment with:

```r
install.packages("renv")
renv::restore()
```

The CIBERSORT package is pinned in `renv.lock` to GitHub commit
`cf2b173f8806a7a93a9e48a176531eaeb8db91f5`.

## Download public expression data

Install the GDC Data Transfer Tool and run from the repository root:

```powershell
Rscript code/download_public_data.R --geo
$env:GDC_CLIENT = "C:/path/to/gdc-client.exe"
Rscript code/download_public_data.R --tcga
```

If TCGA files are stored elsewhere:

```powershell
$env:KIRC_BPS_TCGA_DIR = "D:/data/TCGA-KIRC/Transcriptome_Profiling/Gene_Expression_Quantification"
```

## Required target-database inputs

The database exports are not redistributed in this minimal repository. To run
the target-reconstruction stage, place the exact exports below in
`data/input/target_database_exports/`:

```text
ChEMBL_target.txt
SEA_results.csv
STITCH_interactions.tsv
SuperPred_targets.csv
SwissTargetPrediction.csv
```

The implemented filters are:

- ChEMBL: supplied human UniProt identifiers mapped to official gene symbols;
- SEA: Homo sapiens and P value <= 0.05;
- STITCH: direct BPS–protein edges and combined score >= 0.700;
- SuperPred: probability >= 0.50; model accuracy >= 0.80 defines the
  high-confidence sensitivity flag;
- SwissTargetPrediction: all 100 ranked Homo sapiens predictions; probability
  > 0 defines the high-confidence sensitivity flag.

The reported reconstruction produced 206 unique BPS candidate genes and a
122-gene high-confidence sensitivity set. Because the database exports are not
included, an exact target-set rerun additionally requires the original exports
or their separately archived record.

## Run

```powershell
Rscript code/KIRC_BPS_pipeline.R
```

Outputs are created under `results/reproduced/`. Fixed random seeds are encoded
in the script. The independent-cohort and KEGG modules are skipped when their
separately controlled inputs are absent.

The ESRRB–CIBERSORT figure script requires its separately archived public
source-data CSV. Supply its location with:

```powershell
$env:KIRC_BPS_CIBERSORT_SOURCE_DATA = "D:/source_data/ESRRB_CIBERSORT.csv"
Rscript code/KIRC_BPS_ESRRB_CIBERSORT_correlation_figure.R
```

## Availability and reproducibility boundary

This GitHub repository is a **code archive**, not the sole data repository.
Public dataset accessions, the independent-cohort access route, supplementary
source data, and docking/MD supporting data could be described 
separately in the manuscript and/or a recognised data repository.



