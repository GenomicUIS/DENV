#!/usr/bin/env Rscript
# 03_report.R — verify all outputs and build the methodology + results report.

project_root <- local({
  env <- Sys.getenv("DENGUE_SURAMERICA_ROOT", "")
  if (nzchar(env)) return(normalizePath(env, winslash = "/", mustWork = FALSE))
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]),
                                         winslash = "/", mustWork = FALSE))))
  }
  probe <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  while (!file.exists(file.path(probe, "DESCRIPTION"))) {
    parent <- dirname(probe)
    if (identical(parent, probe)) stop("Cannot find project root.")
    probe <- parent
  }
  probe
})
source(file.path(project_root, "R", "utils.R"))
load_project_helpers(project_root)
load_pipeline_packages(project_root)

resumen_dir <- file.path(project_root, "results", "summary")
qc_dir      <- file.path(project_root, "results", "sequence_curation", "summary")
tiers_root  <- file.path(project_root, "results", "sequence_curation")

paths <- list(
  metadata_csv  = file.path(resumen_dir, "metadata_DENV_south_america_up_to_2026.csv"),
  stats_csv     = file.path(resumen_dir, "genome_statistics_DENV_south_america.csv"),
  missing_txt   = file.path(resumen_dir, "countries_without_records.txt"),
  qc_csv        = file.path(qc_dir, "qc_coding_region_sequences.csv"),
  tiers_csv     = file.path(qc_dir, "depuration_tiers_statistics.csv"),
  review_csv    = file.path(qc_dir, "manual_review.csv"),
  report_md     = file.path(resumen_dir, "methodology_and_results_DENGUE_SOUTH_AMERICA.md"),
  report_txt    = file.path(resumen_dir, "methodology_and_results_DENGUE_SOUTH_AMERICA.txt")
)
inputs <- paths[1:6]
missing <- inputs[!vapply(inputs, file.exists, logical(1L))]
if (length(missing) > 0L) {
  stop("Missing inputs: ", paste(unlist(missing), collapse = ", "))
}

# ---- load -----------------------------------------------------------------
metadata   <- read.csv(paths$metadata_csv, stringsAsFactors = FALSE,
                       check.names = FALSE, fileEncoding = "UTF-8")
stats      <- read.csv(paths$stats_csv,    stringsAsFactors = FALSE,
                       check.names = FALSE, fileEncoding = "UTF-8")
qc         <- read.csv(paths$qc_csv,       stringsAsFactors = FALSE,
                       check.names = FALSE, fileEncoding = "UTF-8")
tier_stats <- read.csv(paths$tiers_csv,    stringsAsFactors = FALSE,
                       check.names = FALSE, fileEncoding = "UTF-8")
review     <- read.csv(paths$review_csv,   stringsAsFactors = FALSE,
                       check.names = FALSE, fileEncoding = "UTF-8")
missing_countries <- readLines(paths$missing_txt, warn = FALSE, encoding = "UTF-8")

# ---- verification ---------------------------------------------------------
verification <- list()
verification$metadata_rows <- nrow(metadata)
verification$stats_sum     <- sum(stats$total_valid_sequences, na.rm = TRUE)
verification$counts_match  <- verification$metadata_rows == verification$stats_sum
verification$combined_fasta_exists <- vapply(
  c("DENV1","DENV2","DENV3","DENV4"), function(s) {
    file.exists(file.path(resumen_dir, "fasta_consolidados_por_serotipo",
                          sprintf("%s_suramerica_hasta_2026.fasta", s)))
  }, logical(1L))

if (!verification$counts_match) {
  stop(sprintf("Verification failed: metadata rows (%d) != stats sum (%d)",
               verification$metadata_rows, verification$stats_sum))
}

# ---- summary tables -------------------------------------------------------
md_table <- function(df) {
  if (is.null(df) || nrow(df) == 0L) return("_(no rows)_")
  df[] <- lapply(df, function(x) gsub("\\|", "/", as.character(x)))
  header <- paste0("| ", paste(names(df), collapse = " | "), " |")
  sep <- paste0("| ", paste(rep("---", ncol(df)), collapse = " | "), " |")
  body <- apply(df, 1L, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  paste(c(header, sep, body), collapse = "\n")
}

serotype_summary <- aggregate(
  cbind(total_valid_sequences, total_complete_genomes, total_complete_cds,
        duplicate_accessions, duplicate_exact_sequences) ~ serotype,
  data = stats, FUN = sum)
serotype_summary$countries_with_records <-
  as.integer(table(metadata$serotype)[serotype_summary$serotype])

country_summary <- as.data.frame(table(metadata$country), stringsAsFactors = FALSE)
names(country_summary) <- c("country", "total_sequences")
country_summary <- country_summary[order(-country_summary$total_sequences), ,
                                   drop = FALSE]

qc_status <- as.data.frame(table(qc$cds_status), stringsAsFactors = FALSE)
names(qc_status) <- c("cds_status", "n")

complete_genome_trimmed <- sum(qc$completeness == "complete_genome" &
                                 qc$length_delta > 0L, na.rm = TRUE)
short_cds_n <- sum(qc$coding_length < 9000L, na.rm = TRUE)

tier_counts <- tier_stats

# ---- interpretation numbers (computed) ------------------------------------
n_metadata_total <- nrow(metadata)
n_level3 <- tier_stats$total_sequences[tier_stats$tier_id == "nivel_3_ACGT_puro"]
n_level4 <- tier_stats$total_sequences[tier_stats$tier_id == "nivel_4_auto_filtered"]
n_level5 <- tier_stats$total_sequences[tier_stats$tier_id == "nivel_5_manual_curated"]
n_review <- nrow(review)

# ---- write markdown -------------------------------------------------------
lines <- c(
  "# DENGUE_SURAMERICA — Methodology and results",
  "",
  sprintf("Date: %s", format(Sys.Date(), "%Y-%m-%d")),
  sprintf("Project root: `%s`", project_root),
  "",
  "## 1. Objective",
  "",
  "Assemble a reproducible set of dengue virus sequences from South America,",
  "organized by serotype and country, with FASTA and full GenBank files,",
  "consolidated metadata, download statistics, and QC tiers suitable for",
  "phylogenetic analysis.",
  "",
  "Time range: 1900-2026 inclusive. Countries: Argentina, Bolivia, Brasil,",
  "Chile, Colombia, Ecuador, Guyana, Paraguay, Peru, Surinam, Uruguay,",
  "Venezuela. Serotypes: DENV1-DENV4.",
  "",
  "## 2. Scripts",
  "",
  "- `scripts/00_setup.R` — install required R packages into `.Rlib/`.",
  "- `scripts/01_download.R` — query NCBI, download, filter, deduplicate.",
  "- `scripts/02_qc_tiers.R` — extract CDS, QC, and build depuration tiers.",
  "- `scripts/03_report.R` — verify outputs and write this report.",
  "- `scripts/run_all.R` — runs 00 -> 03 in order.",
  "",
  "## 3. Methods",
  "",
  "### 3.1 Query",
  "",
  "Each (serotype, country) pair is queried independently via NCBI E-utilities",
  "(`esearch` -> `efetch`). The query restricts to complete genomes / complete",
  "CDS / complete coding sequences, bounds `PDAT` to 1900-2026, and excludes",
  "records whose titles indicate partial, fragment, amplicon, incomplete,",
  "synthetic, artificial, construct, or patent content.",
  "",
  "### 3.2 Record-level filters",
  "",
  "A record is retained only when:",
  "",
  "- `geo_loc_name` / `country` parses to the target South American country.",
  "- The definition/organism matches dengue (`Dengue virus`, `DENV`, or",
  "  `Orthoflavivirus denguei`).",
  "- If a serotype can be inferred from the header, it agrees with the query.",
  "- The definition does not indicate partial / fragment / synthetic.",
  "- The accession has not been seen before.",
  "- The sequence does not duplicate a previously accepted sequence.",
  "",
  "### 3.3 File layout",
  "",
  "Per accepted record:",
  "",
  "- `<SEROTYPE>/fasta/<accession>.fasta`",
  "- `<SEROTYPE>/genbank_full/<accession>.gb`",
  "- `results/by_country/<country>/<serotype>/{fasta,genbank_full,metadata}/`",
  "",
  "### 3.4 CDS extraction and QC",
  "",
  "For each accepted record the GenBank `CDS` annotation is parsed and the",
  "coding sequence is extracted (5'UTR and 3'UTR removed). When multiple CDS",
  "features exist, the polyprotein CDS is preferred, then the longest.",
  "",
  "Sequences are then scanned for characters: A, C, G, T, N, gaps, IUPAC",
  "ambiguity codes, and other non-ACGT characters.",
  "",
  "### 3.5 Depuration tiers",
  "",
  "| Tier | Rule |",
  "|---|---|",
  "| `nivel_1_metadata_complete_or_cds` | Metadata-accepted, CDS-trimmed |",
  "| `nivel_2_sin_N_sin_gaps` | Level 1 minus N and gaps |",
  "| `nivel_3_ACGT_puro` | Level 2 minus any non-ACGT character |",
  "| `nivel_4_auto_filtered` | Level 3 minus sequences auto-flagged as short CDS or non-dengue |",
  "| `nivel_5_manual_curated` | Level 3 minus sequences the researcher marked `exclude` in `manual_review.csv` |",
  "",
  "### 3.6 Manual review workflow",
  "",
  "`scripts/02_qc_tiers.R` writes `results/sequence_curation/summary/manual_review.csv`",
  "with every auto-flagged sequence, defaulting to `decision = exclude`. The",
  "researcher flips individual rows to `keep` and re-runs the script. Level 5",
  "is then rebuilt from those decisions. Level 4 (automatic filter) is left",
  "unchanged, so the two tiers are always comparable.",
  "",
  "## 4. Results — download",
  "",
  sprintf("Total accepted sequences in metadata: **%d**.", n_metadata_total),
  "",
  "### 4.1 Per serotype",
  "",
  md_table(serotype_summary),
  "",
  "### 4.2 Per country",
  "",
  md_table(country_summary),
  "",
  "### 4.3 Countries with zero records per serotype",
  "",
  paste(sprintf("- %s", missing_countries), collapse = "\n"),
  "",
  "## 5. Results — CDS extraction and QC",
  "",
  sprintf("Records evaluated in QC: %d.", nrow(qc)),
  sprintf("Records with a CDS extracted: %d.",
          sum(qc$cds_status == "cds_extracted", na.rm = TRUE)),
  sprintf("Complete genomes with UTR removed: %d.", complete_genome_trimmed),
  sprintf("Sequences with coding length < 9000 nt: %d.", short_cds_n),
  "",
  "### 5.1 CDS extraction status",
  "",
  md_table(qc_status),
  "",
  "### 5.2 Tier counts",
  "",
  md_table(tier_counts),
  "",
  "## 6. Sequences under manual review",
  "",
  sprintf("Sequences flagged for manual review: %d.", n_review),
  sprintf("Marked `keep` by the researcher: %d.",
          sum(tolower(trimws(review$decision)) == "keep", na.rm = TRUE)),
  sprintf("Marked `exclude` (or default): %d.",
          sum(tolower(trimws(review$decision)) != "keep", na.rm = TRUE)),
  "",
  "Auto-flagged sequences are kept under",
  "`depuracion_secuencias/revision_auto/` for offline inspection.",
  "",
  "## 7. Interpretation",
  "",
  sprintf(paste0("The pipeline accepted **%d** sequences by GenBank metadata. ",
                 "After CDS extraction and character filtering, the strict ",
                 "ACGT-only tier contains **%d** sequences. The automatic ",
                 "filter (level 4) retains **%d**, and the manual-curated ",
                 "tier (level 5) contains **%d** sequences."),
          n_metadata_total, n_level3, n_level4, n_level5),
  "",
  "For alignment (MAFFT) and tree inference (IQ-TREE), the recommended",
  "starting point is `nivel_4_auto_filtered`, or `nivel_5_manual_curated` once",
  "the review CSV has been curated.",
  "",
  "## 8. Limitations",
  "",
  "- Serotype classification for each record is derived from the header and",
  "  organism fields; no independent lineage validation is performed here.",
  "- Subnational geographic coordinates are not reliable for most GenBank",
  "  records.",
  "- Level 3 controls characters, but does not inspect alignment quality,",
  "  coverage, or recombination. Confirm lineage / serotype with an",
  "  independent tool (e.g., Genome Detective) before publication.",
  "",
  "## 9. Verification",
  "",
  sprintf("- Metadata rows: %d", verification$metadata_rows),
  sprintf("- Sum of per-query stats: %d", verification$stats_sum),
  sprintf("- Match: %s", verification$counts_match),
  sprintf("- Consolidated FASTA files present: %s",
          paste(names(verification$combined_fasta_exists)[
            verification$combined_fasta_exists], collapse = ", ")),
  "",
  "## 10. Main output files",
  "",
  "- `resumen/metadata_DENV_suramerica_hasta_2026.csv`",
  "- `resumen/estadisticas_genomas_DENV_suramerica.csv`",
  "- `results/sequence_curation/summary/qc_secuencias_coding_region.csv`",
  "- `results/sequence_curation/summary/estadistic_level_depuration.csv`",
  "- `results/sequence_curation/summary/manual_review.csv`",
  "- `depuracion_secuencias/revision_auto/`",
  "- `depuracion_secuencias/nivel_{1..5}_*/`"
)

writeLines(lines, paths$report_md, useBytes = TRUE)

plain <- lines
plain <- gsub("`", "", plain, fixed = TRUE)
plain <- sub("^### ", "", plain)
plain <- sub("^## ",  "", plain)
plain <- sub("^# ",   "", plain)
writeLines(plain, paths$report_txt, useBytes = TRUE)

cat("Report written:\n")
cat("  -", paths$report_md, "\n")
cat("  -", paths$report_txt, "\n")