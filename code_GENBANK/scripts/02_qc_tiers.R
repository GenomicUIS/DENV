#!/usr/bin/env Rscript
# 02_qc_tiers.R - extract CDS, apply QC rules, build depuration tiers.
#
# Tiers produced:
#   nivel_1_metadata_complete_or_cds  - accepted by metadata, CDS-trimmed
#   nivel_2_sin_N_sin_gaps            - level 1 minus N and gaps
#   nivel_3_ACGT_puro                 - level 2 minus any non-ACGT char
#   nivel_4_auto_filtered             - level 3 minus auto-flagged records
#   nivel_5_manual_curated            - level 3 minus records the researcher
#                                       marked "exclude" in manual_review.csv
#
# Auto-flagged records (short CDS, non-dengue taxonomy) are copied to
#   depuracion_secuencias/revision_auto/ for inspection, and listed in
#   depuracion_secuencias/resumen/manual_review.csv. The researcher edits the
#   `decision` column (keep / exclude) and re-runs this script.

# ---- bootstrap ------------------------------------------------------------
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

# ---- configuration --------------------------------------------------------
tiers_root  <- file.path(project_root, "results", "sequence_curation")
summary_dir <- file.path(tiers_root, "summary")
review_dir  <- file.path(tiers_root, "auto_review")
metadata_csv <- file.path(project_root, "results", "summary",
                          "metadata_DENV_south_america_up_to_2026.csv")
min_coding_length <- as.integer(Sys.getenv("DENV_MIN_CODING_LENGTH", "9000"))
if (is.na(min_coding_length) || min_coding_length <= 0L) min_coding_length <- 9000L

log <- new_logger(file.path(summary_dir, "log_qc_tiers.txt"))

# ---- optional cleanup -----------------------------------------------------
clean <- tolower(Sys.getenv("DENV_CLEAN_OUTPUT", "false")) %in% c("1", "true", "yes")
if (clean && dir.exists(tiers_root)) {
  resolved <- normalizePath(tiers_root, winslash = "/", mustWork = FALSE)
  if (startsWith(resolved, paste0(project_root, "/")) &&
      basename(resolved) == "depuracion_secuencias") {
    log$warn("Removing %s", resolved)
    unlink(resolved, recursive = TRUE, force = TRUE)
  }
}
ensure_dir(summary_dir)

if (!file.exists(metadata_csv)) stop("Metadata file not found: ", metadata_csv)
metadata <- read.csv(metadata_csv, stringsAsFactors = FALSE, check.names = FALSE,
                     fileEncoding = "UTF-8")
log$info("Loaded %d metadata rows", nrow(metadata))

# ---- per-record QC --------------------------------------------------------
qc_rows <- vector("list", nrow(metadata))
coding_sequences <- list()   # accession -> CDS/ORF sequence

for (i in seq_len(nrow(metadata))) {
  row <- metadata[i, , drop = FALSE]
  fasta_path <- file.path(project_root, row$fasta_relative_path[[1L]])
  gb_path    <- file.path(project_root, row$genbank_relative_path[[1L]])
  
  if (!file.exists(fasta_path) || !file.exists(gb_path)) {
    log$warn("Missing files for %s - skipping", row$accession[[1L]])
    next
  }
  
  original <- read_fasta(fasta_path)
  gb_text  <- paste(readLines(gb_path, warn = FALSE, encoding = "UTF-8"),
                    collapse = "\n")
  record   <- parse_genbank_record(gb_text)
  cds      <- extract_cds_sequence(record, original$sequence)
  coding   <- toupper(cds$sequence)
  qc       <- sequence_qc(coding)
  
  coding_sequences[[row$accession[[1L]]]] <- coding
  
  utr5 <- if (!is.na(cds$start)) max(cds$start - 1L, 0L) else NA_integer_
  utr3 <- if (!is.na(cds$end))   max(nchar(original$sequence) - cds$end, 0L) else NA_integer_
  
  auto_reasons <- character(0)
  if (!is_dengue_record(record$definition, record$organism, record$source)) {
    auto_reasons <- c(auto_reasons, "taxonomy_not_dengue")
  }
  if (!is.na(qc$coding_length) && qc$coding_length < min_coding_length) {
    auto_reasons <- c(auto_reasons,
                      sprintf("coding_length_below_%d_nt", min_coding_length))
  }
  if (identical(cds$status, "cds_parse_failed")) {
    auto_reasons <- c(auto_reasons, "cds_parse_failed")
  }
  
  qc_row <- data.frame(
    accession                = row$accession[[1L]],
    serotype                 = row$serotype[[1L]],
    country                  = row$country[[1L]],
    year                     = row$year[[1L]],
    completeness             = row$completeness[[1L]],
    definition               = record$definition,
    organism                 = record$organism,
    original_fasta_length    = nchar(original$sequence),
    sequence_used            = cds$source,
    cds_status               = cds$status,
    cds_location             = cds$location,
    cds_product              = cds$product,
    cds_start                = cds$start,
    cds_end                  = cds$end,
    cds_fuzzy                = cds$fuzzy,
    utr5_removed             = utr5,
    utr3_removed             = utr3,
    length_delta             = nchar(original$sequence) - nchar(coding),
    coding_length            = qc$coding_length,
    n_bases                  = qc$N,
    gaps                     = qc$gaps,
    iupac_ambiguous          = qc$iupac_ambiguous,
    other_non_acgt           = qc$other_non_acgt,
    non_acgt_total           = qc$non_acgt_total,
    pct_N                    = qc$pct_N,
    pct_gaps                 = qc$pct_gaps,
    pct_non_acgt             = qc$pct_non_acgt,
    level1_pass              = qc$coding_length > 0L,
    level2_pass              = qc$coding_length > 0L && qc$N == 0L && qc$gaps == 0L,
    level3_pass              = qc$coding_length > 0L && qc$non_acgt_total == 0L,
    auto_flag_required       = length(auto_reasons) > 0L,
    auto_flag_reason         = paste(auto_reasons, collapse = ";"),
    stringsAsFactors = FALSE
  )
  qc_rows[[i]] <- qc_row
  
  if (i %% 100L == 0L || i == nrow(metadata)) {
    log$info("QC processed %d/%d", i, nrow(metadata))
  }
}

# Some rows may have been skipped, leaving NULL entries. Drop them before rbind.
qc_rows <- qc_rows[!vapply(qc_rows, is.null, logical(1L))]
if (length(qc_rows) == 0L) {
  stop("No records survived QC. Check that 01_download.R ran successfully.")
}
qc_df <- do.call(rbind, qc_rows)
write.csv(qc_df, file.path(summary_dir, "qc_coding_region_sequences.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# ---- tier writer ----------------------------------------------------------
tier_ids <- c("nivel_1_metadata_complete_or_cds",
              "nivel_2_sin_N_sin_gaps",
              "nivel_3_ACGT_puro",
              "nivel_4_auto_filtered",
              "nivel_5_manual_curated")

tier_labels <- c(
  nivel_1_metadata_complete_or_cds = "Level 1: metadata-accepted, CDS-trimmed",
  nivel_2_sin_N_sin_gaps           = "Level 2: no N, no gaps",
  nivel_3_ACGT_puro                = "Level 3: pure ACGT",
  nivel_4_auto_filtered            = "Level 4: auto-filtered (short CDS or non-dengue)",
  nivel_5_manual_curated           = "Level 5: manual curated (researcher decisions)"
)

tier_paths <- function(tier_dir, acc, ser, ctry) {
  list(
    fasta         = file.path(tier_dir, ser, "fasta",        paste0(acc, ".fasta")),
    gb            = file.path(tier_dir, ser, "genbank_full", paste0(acc, ".gb")),
    country_fasta = file.path(tier_dir, "by_country", ctry, ser, "fasta",
                              paste0(acc, ".fasta")),
    country_gb    = file.path(tier_dir, "by_country", ctry, ser, "genbank_full",
                              paste0(acc, ".gb"))
  )
}

write_tier <- function(tier_id, accessions) {
  tier_dir <- file.path(tiers_root, tier_id)
  if (dir.exists(tier_dir)) unlink(tier_dir, recursive = TRUE, force = TRUE)
  ensure_dir(file.path(tier_dir, "resumen"))
  
  subset <- qc_df[qc_df$accession %in% accessions, , drop = FALSE]
  if (nrow(subset) == 0L) {
    log$warn("Tier %s has zero sequences", tier_id)
    write.csv(subset,
              file.path(tier_dir, "resumen",
                        sprintf("metadata_%s.csv", tier_id)),
              row.names = FALSE)
    return(invisible(subset))
  }
  
  meta_rows <- vector("list", nrow(subset))
  for (i in seq_len(nrow(subset))) {
    row <- subset[i, , drop = FALSE]
    acc <- row$accession[[1L]]
    seq <- coding_sequences[[acc]]
    if (is.null(seq) || !nzchar(seq)) next
    
    header <- sprintf("%s|serotype=%s|country=%s|year=%s|tier=%s|length=%d",
                      clean_header_value(acc),
                      clean_header_value(row$serotype[[1L]]),
                      clean_header_value(row$country[[1L]]),
                      clean_header_value(row$year[[1L]]),
                      tier_id, nchar(seq))
    paths <- tier_paths(tier_dir, acc, row$serotype[[1L]], row$country[[1L]])
    for (p in paths) ensure_dir(dirname(p))
    
    write_fasta(seq, header, paths$fasta)
    write_fasta(seq, header, paths$country_fasta)
    
    gb_match <- match(acc, metadata$accession)
    if (!is.na(gb_match)) {
      src_gb <- file.path(project_root, metadata$genbank_relative_path[gb_match])
      if (file.exists(src_gb)) {
        file.copy(src_gb, paths$gb,         overwrite = TRUE)
        file.copy(src_gb, paths$country_gb, overwrite = TRUE)
      } else {
        log$warn("GenBank source missing for %s", acc)
      }
    } else {
      log$warn("Accession %s not in metadata index", acc)
    }
    
    meta_rows[[i]] <- cbind(row, data.frame(
      tier_id = tier_id,
      tier_label = tier_labels[[tier_id]],
      tier_fasta_relative_path = rel_path(paths$fasta, tier_dir),
      tier_genbank_relative_path = rel_path(paths$gb, tier_dir),
      stringsAsFactors = FALSE))
  }
  
  meta_rows <- meta_rows[!vapply(meta_rows, is.null, logical(1L))]
  meta_df <- if (length(meta_rows) > 0L) do.call(rbind, meta_rows) else subset
  write.csv(meta_df,
            file.path(tier_dir, "resumen", sprintf("metadata_%s.csv", tier_id)),
            row.names = FALSE, fileEncoding = "UTF-8")
  invisible(meta_df)
}

# ---- build levels 1-3 -----------------------------------------------------
level1_acc <- qc_df$accession[qc_df$level1_pass]
level2_acc <- qc_df$accession[qc_df$level2_pass]
level3_acc <- qc_df$accession[qc_df$level3_pass]

log$info("Level 1: %d | Level 2: %d | Level 3: %d",
         length(level1_acc), length(level2_acc), length(level3_acc))

write_tier("nivel_1_metadata_complete_or_cds", level1_acc)
write_tier("nivel_2_sin_N_sin_gaps",           level2_acc)
write_tier("nivel_3_ACGT_puro",                level3_acc)

# ---- auto-flagged review package ------------------------------------------
auto_flagged <- qc_df[qc_df$auto_flag_required & qc_df$accession %in% level3_acc,
                      , drop = FALSE]

if (nrow(auto_flagged) > 0L) {
  if (dir.exists(review_dir)) unlink(review_dir, recursive = TRUE, force = TRUE)
  for (i in seq_len(nrow(auto_flagged))) {
    row <- auto_flagged[i, , drop = FALSE]
    acc <- row$accession[[1L]]
    seq <- coding_sequences[[acc]]
    header <- sprintf("%s|serotype=%s|country=%s|auto_flag=%s",
                      clean_header_value(acc),
                      clean_header_value(row$serotype[[1L]]),
                      clean_header_value(row$country[[1L]]),
                      clean_header_value(row$auto_flag_reason[[1L]]))
    write_fasta(seq, header,
                file.path(review_dir, row$serotype[[1L]], "fasta",
                          paste0(acc, ".fasta")))
    gb_match <- match(acc, metadata$accession)
    if (!is.na(gb_match)) {
      src_gb <- file.path(project_root, metadata$genbank_relative_path[gb_match])
      if (file.exists(src_gb)) {
        tgt <- file.path(review_dir, row$serotype[[1L]], "genbank_full",
                         paste0(acc, ".gb"))
        ensure_dir(dirname(tgt))
        file.copy(src_gb, tgt, overwrite = TRUE)
      }
    }
  }
  write.csv(auto_flagged,
            file.path(review_dir, "auto_flagged.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
}

# ---- manual review CSV (persistent) ---------------------------------------
manual_review_path <- file.path(summary_dir, "manual_review.csv")
template <- data.frame(
  accession        = auto_flagged$accession,
  serotype         = auto_flagged$serotype,
  country          = auto_flagged$country,
  coding_length    = auto_flagged$coding_length,
  auto_flag_reason = auto_flagged$auto_flag_reason,
  definition       = auto_flagged$definition,
  decision         = rep("exclude", nrow(auto_flagged)),
  notes            = rep("", nrow(auto_flagged)),
  stringsAsFactors = FALSE
)

if (file.exists(manual_review_path)) {
  previous <- read.csv(manual_review_path, stringsAsFactors = FALSE)
  for (j in seq_len(nrow(template))) {
    k <- match(template$accession[[j]], previous$accession)
    if (!is.na(k)) {
      if ("decision" %in% names(previous)) {
        template$decision[[j]] <- as.character(previous$decision[[k]])
      }
      if ("notes" %in% names(previous)) {
        template$notes[[j]] <- as.character(previous$notes[[k]])
      }
    }
  }
}
write.csv(template, manual_review_path, row.names = FALSE, fileEncoding = "UTF-8")

# ---- level 4 (auto-filtered) ----------------------------------------------
level4_acc <- setdiff(level3_acc, auto_flagged$accession)
write_tier("nivel_4_auto_filtered", level4_acc)

# ---- level 5 (manual curated) ---------------------------------------------
# Default: same as level 4. If the researcher marked `keep` for some flagged
# sequences, those are added back into level 5.
keep_acc <- template$accession[tolower(trimws(template$decision)) == "keep"]
level5_acc <- union(level4_acc, keep_acc)
write_tier("nivel_5_manual_curated", level5_acc)

# ---- statistics -----------------------------------------------------------
tier_counts <- data.frame(
  tier_id = tier_ids,
  tier_label = unname(tier_labels[tier_ids]),
  total_sequences = c(length(level1_acc), length(level2_acc),
                      length(level3_acc), length(level4_acc), length(level5_acc)),
  stringsAsFactors = FALSE
)
write.csv(tier_counts, file.path(summary_dir, "depuration_tiers_statistics.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write_excel(file.path(summary_dir, "depuration_tiers_statistics.xlsx"),
            list(by_tier = tier_counts, qc = qc_df))

# ---- README ---------------------------------------------------------------
writeLines(c(
  "# Sequence depuration tiers",
  "",
  "Levels 1-3 are automatic. Levels 4-5 filter further:",
  "",
  "- **nivel_4_auto_filtered**: level 3 minus sequences flagged automatically",
  sprintf("  (coding length < %d nt or non-dengue taxonomy).", min_coding_length),
  "- **nivel_5_manual_curated**: level 3 minus the sequences the researcher",
  "  marked `exclude` in `resumen/manual_review.csv`.",
  "",
  "## Manual review workflow",
  "",
  "1. Run this script once. It writes `resumen/manual_review.csv` listing every",
  "   auto-flagged sequence with a `decision` column set to `exclude`.",
  "2. Open the CSV, review each row, and change `decision` to `keep` for",
  "   sequences you want to retain. Use `notes` to justify your choice.",
  "3. Re-run this script. Level 5 is rebuilt from your decisions. Level 4 is",
  "   unchanged (it reflects the automatic rules only).",
  "",
  "FASTA/GenBank files for auto-flagged sequences are also kept under",
  "`revision_auto/` for offline inspection."
), file.path(tiers_root, "README.md"), useBytes = TRUE)

log$info("Done. Tier counts: L1=%d L2=%d L3=%d L4=%d L5=%d",
         length(level1_acc), length(level2_acc), length(level3_acc),
         length(level4_acc), length(level5_acc))
cat("\nTier summary:\n")
print(tier_counts, row.names = FALSE)
cat("\nManual review file:", manual_review_path, "\n")