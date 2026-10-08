#!/usr/bin/env Rscript
# 01_download.R — download dengue records from NCBI GenBank for South America.
#
# Outputs:
#   resumen/metadata_DENV_suramerica_hasta_<year>.{csv,xlsx}
#   resumen/estadisticas_genomas_DENV_suramerica.{csv,xlsx}
#   resumen/paises_sin_serotipos_encontrados.txt
#   resumen/log_descarga.txt
#   resumen/log_filtros_registros.csv
#   resumen/fasta_consolidados_por_serotipo/<SEROTYPE>_suramerica_hasta_<year>.fasta
#   <SEROTYPE>/{fasta,genbank_full,metadata}/
#   por_pais/<country>/<serotype>/{fasta,genbank_full,metadata}/

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
cfg <- list(
  root         = project_root,
  resumen_dir  = file.path(project_root, "results", "summary"),
  country_view = file.path(project_root, "results", "by_country"),
  max_year     = as.integer(Sys.getenv("DENV_MAX_YEAR", "2026")),
  batch_size   = as.integer(Sys.getenv("DENV_BATCH_SIZE", "100")),
  retmax       = as.integer(Sys.getenv("DENV_RETMAX", "500")),
  max_ids      = as.integer(Sys.getenv("DENV_MAX_IDS_PER_QUERY", "0")),
  dedupe_seq   = tolower(Sys.getenv("DENV_DEDUPE_SEQUENCE", "true")) %in%
    c("1", "true", "yes"),
  clean_output = tolower(Sys.getenv("DENV_CLEAN_OUTPUT", "false")) %in%
    c("1", "true", "yes")
)
if (is.na(cfg$max_year))   cfg$max_year   <- 2026L
if (is.na(cfg$batch_size)) cfg$batch_size <- 100L
if (is.na(cfg$retmax))     cfg$retmax     <- 500L
if (is.na(cfg$max_ids) || cfg$max_ids <= 0L) cfg$max_ids <- NULL

ensure_dir(cfg$resumen_dir)
ensure_dir(file.path(cfg$resumen_dir, "consolidated_fasta_by_serotype"))

log <- new_logger(file.path(cfg$resumen_dir, "download_log.txt"))
log$info("Project root: %s", project_root)
log$info("Max year: %d", cfg$max_year)
log$info("Batch size: %d", cfg$batch_size)
log$info("Dedupe exact sequences: %s", cfg$dedupe_seq)
log$info("Clean output: %s", cfg$clean_output)

# ---- optional cleanup -----------------------------------------------------
if (cfg$clean_output) {
  targets <- c(file.path(project_root, "results", c("DENV1", "DENV2", "DENV3", "DENV4")),
               cfg$country_view)
  for (t in targets) {
    if (!dir.exists(t)) next
    resolved <- normalizePath(t, winslash = "/", mustWork = FALSE)
    if (startsWith(resolved, paste0(project_root, "/"))) {
      log$warn("Removing %s", resolved)
      unlink(resolved, recursive = TRUE, force = TRUE)
    }
  }
  summary_files <- c("metadata_DENV_south_america_up_to_2026.csv",
                     "metadata_DENV_south_america_up_to_2026.xlsx",
                     "genome_statistics_DENV_south_america.csv",
                     "genome_statistics_DENV_south_america.xlsx",
                     "countries_without_records.txt",
                     "records_filter_log.csv",
                     "verification_statistics_vs_metadata.csv")
  for (f in summary_files) {
    p <- file.path(cfg$resumen_dir, f)
    if (file.exists(p)) unlink(p, force = TRUE)
  }
}

# ---- query inputs ---------------------------------------------------------
countries_df <- data.frame(
  country = c("Argentina", "Bolivia", "Brazil", "Chile", "Colombia",
              "Ecuador", "Guyana", "Paraguay", "Peru", "Suriname",
              "Uruguay", "Venezuela"),
  query_alias = c(
    '"Argentina"[All Fields]',
    '"Bolivia"[All Fields]',
    '("Brazil"[All Fields] OR "Brasil"[All Fields])',
    '"Chile"[All Fields]',
    '"Colombia"[All Fields]',
    '"Ecuador"[All Fields]',
    '"Guyana"[All Fields]',
    '"Paraguay"[All Fields]',
    '"Peru"[All Fields]',
    '("Suriname"[All Fields] OR "Surinam"[All Fields])',
    '"Uruguay"[All Fields]',
    '"Venezuela"[All Fields]'
  ),
  stringsAsFactors = FALSE
)

serotypes <- list(
  DENV1 = c('"Dengue virus type 1"[Organism]', '"Dengue virus 1"[Organism]',
            'DENV-1[All Fields]', 'DENV1[All Fields]'),
  DENV2 = c('"Dengue virus type 2"[Organism]', '"Dengue virus 2"[Organism]',
            'DENV-2[All Fields]', 'DENV2[All Fields]'),
  DENV3 = c('"Dengue virus type 3"[Organism]', '"Dengue virus 3"[Organism]',
            'DENV-3[All Fields]', 'DENV3[All Fields]'),
  DENV4 = c('"Dengue virus type 4"[Organism]', '"Dengue virus 4"[Organism]',
            'DENV-4[All Fields]', 'DENV4[All Fields]')
)

sel_countries <- trimws(strsplit(Sys.getenv("DENV_COUNTRIES", ""), ",",
                                 fixed = TRUE)[[1L]])
sel_countries <- sel_countries[nzchar(sel_countries)]
if (length(sel_countries) > 0L) {
  countries_df <- countries_df[countries_df$country %in% sel_countries, ,
                               drop = FALSE]
}
sel_serotypes <- trimws(strsplit(Sys.getenv("DENV_SEROTYPES", ""), ",",
                                 fixed = TRUE)[[1L]])
sel_serotypes <- sel_serotypes[nzchar(sel_serotypes)]
if (length(sel_serotypes) > 0L) {
  serotypes <- serotypes[names(serotypes) %in% sel_serotypes]
}
if (nrow(countries_df) == 0L) stop("No countries selected.")
if (length(serotypes) == 0L)    stop("No serotypes selected.")

# ---- directory skeleton ---------------------------------------------------
for (s in names(serotypes)) {
  for (sub in c("fasta", "genbank_full", "metadata")) {
    ensure_dir(file.path(project_root, "results", s, sub))
  }
}
for (i in seq_len(nrow(countries_df))) {
  for (s in names(serotypes)) {
    for (sub in c("fasta", "genbank_full", "metadata")) {
      ensure_dir(file.path(cfg$country_view, countries_df$country[[i]], s, sub))
    }
  }
}

# ---- query builder --------------------------------------------------------
build_query <- function(serotype_terms, country_alias, max_year) {
  serotype_block <- paste(serotype_terms, collapse = " OR ")
  completeness_block <- paste(
    '"complete genome"[All Fields]',
    '"complete cds"[All Fields]',
    '"complete coding sequence"[All Fields]',
    sep = " OR "
  )
  exclusion_block <- paste(
    'partial[Title]', 'fragment[Title]', 'amplicon[Title]',
    'incomplete[Title]', 'synthetic[Title]', 'artificial[Title]',
    'construct[Title]', 'patent[All Fields]',
    sep = " OR "
  )
  sprintf('(%s) AND (%s) AND (%s) AND ("1900/01/01"[PDAT] : "%d/12/31"[PDAT]) NOT (%s)',
          serotype_block, completeness_block, country_alias,
          max_year, exclusion_block)
}

write_record_files <- function(record_text, record, country, serotype) {
  acc <- record$accession
  serotype_dir <- file.path(project_root, "results", serotype)
  country_dir  <- file.path(cfg$country_view, country, serotype)
  paths <- list(
    gb            = file.path(serotype_dir, "genbank_full", paste0(acc, ".gb")),
    fasta         = file.path(serotype_dir, "fasta",        paste0(acc, ".fasta")),
    country_gb    = file.path(country_dir, "genbank_full",  paste0(acc, ".gb")),
    country_fasta = file.path(country_dir, "fasta",         paste0(acc, ".fasta"))
  )
  for (p in paths) ensure_dir(dirname(p))
  
  writeLines(record_text, paths$gb, useBytes = TRUE)
  header <- sprintf("%s|serotype=%s|country=%s|year=%s",
                    clean_header_value(record$accession),
                    clean_header_value(serotype),
                    clean_header_value(country),
                    clean_header_value(if (is.na(record$year)) "NA" else record$year))
  write_fasta(record$sequence, header, paths$fasta)
  file.copy(paths$gb,    paths$country_gb,    overwrite = TRUE)
  file.copy(paths$fasta, paths$country_fasta, overwrite = TRUE)
  paths
}

# ---- record processing ----------------------------------------------------
seen_accessions <- new.env(parent = emptyenv())
seen_sequences  <- new.env(parent = emptyenv())

process_record <- function(parsed, record_text, serotype, target_country) {
  row <- genbank_to_metadata_row(parsed, serotype, target_country)
  
  if (!identical(row$country, target_country)) {
    return(list(status = "rejected",
                reason = sprintf("country_mismatch:%s", row$country)))
  }
  if (!is_dengue_record(parsed$definition, parsed$organism, parsed$source)) {
    return(list(status = "rejected", reason = "taxonomy_not_dengue"))
  }
  detected <- parse_serotype(parsed$definition, parsed$organism, parsed$source)
  if (!is.na(detected) && !identical(detected, serotype)) {
    return(list(status = "rejected",
                reason = sprintf("serotype_mismatch:%s", detected)))
  }
  hreason <- header_exclusion_reason(row, cfg$max_year)
  if (!is.na(hreason)) {
    return(list(status = "rejected", reason = hreason))
  }
  if (exists(parsed$accession, envir = seen_accessions, inherits = FALSE)) {
    return(list(status = "rejected", reason = "duplicate_accession"))
  }
  seq <- parsed$sequence
  if (!nzchar(seq)) {
    return(list(status = "rejected", reason = "empty_sequence"))
  }
  if (cfg$dedupe_seq && exists(seq, envir = seen_sequences, inherits = FALSE)) {
    first_acc <- get(seq, envir = seen_sequences)
    return(list(status = "rejected",
                reason = sprintf("duplicate_exact_sequence_of:%s", first_acc)))
  }
  
  assign(parsed$accession, TRUE, envir = seen_accessions)
  if (cfg$dedupe_seq) assign(seq, parsed$accession, envir = seen_sequences)
  
  files <- write_record_files(record_text, parsed, target_country, serotype)
  row$genbank_relative_path         <- rel_path(files$gb,            project_root)
  row$fasta_relative_path           <- rel_path(files$fasta,         project_root)
  row$country_genbank_relative_path <- rel_path(files$country_gb,    project_root)
  row$country_fasta_relative_path   <- rel_path(files$country_fasta, project_root)
  row$source_database <- "NCBI GenBank"
  row$notes <- "accepted"
  
  list(status = "accepted", reason = "accepted", metadata = row)
}

# ---- main loop ------------------------------------------------------------
client <- new_ncbi_client()
all_metadata   <- list()
all_filters    <- list()
all_stats      <- list()
query_failures <- 0L

for (serotype in names(serotypes)) {
  for (i in seq_len(nrow(countries_df))) {
    target_country <- countries_df$country[[i]]
    query <- build_query(serotypes[[serotype]], countries_df$query_alias[[i]],
                         cfg$max_year)
    log$info("Query | %s | %s", serotype, target_country)
    
    ids <- tryCatch(
      esearch_all_ids(client, query, retmax = cfg$retmax, max_ids = cfg$max_ids),
      error = function(e) {
        query_failures <<- query_failures + 1L
        log$error("ESearch failed | %s | %s | %s",
                  serotype, target_country, conditionMessage(e))
        character(0)
      }
    )
    log$info("Found %d IDs", length(ids))
    
    counts <- list(accepted = 0L, dup_acc = 0L, dup_seq = 0L,
                   complete_genome = 0L, complete_cds = 0L,
                   lengths = integer(0))
    
    if (length(ids) > 0L) {
      records <- tryCatch(
        efetch_genbank(client, ids, batch_size = cfg$batch_size),
        error = function(e) {
          log$error("EFetch failed | %s | %s | %s",
                    serotype, target_country, conditionMessage(e))
          list()
        }
      )
      log$info("Fetched %d records", length(records))
      
      for (record_text in records) {
        parsed <- tryCatch(parse_genbank_record(record_text),
                           error = function(e) NULL)
        if (is.null(parsed) || is.na(parsed$accession)) {
          all_filters[[length(all_filters) + 1L]] <- data.frame(
            serotype = serotype, query_country = target_country,
            accession = NA_character_, status = "rejected",
            reason = "parse_error", stringsAsFactors = FALSE)
          next
        }
        result <- process_record(parsed, record_text, serotype, target_country)
        all_filters[[length(all_filters) + 1L]] <- data.frame(
          serotype = serotype, query_country = target_country,
          accession = parsed$accession, status = result$status,
          reason = result$reason, stringsAsFactors = FALSE)
        
        if (identical(result$status, "accepted")) {
          all_metadata[[length(all_metadata) + 1L]] <- result$metadata
          counts$accepted <- counts$accepted + 1L
          if (identical(result$metadata$completeness, "complete_genome")) {
            counts$complete_genome <- counts$complete_genome + 1L
          } else if (identical(result$metadata$completeness, "complete_cds")) {
            counts$complete_cds <- counts$complete_cds + 1L
          }
          counts$lengths <- c(counts$lengths, parsed$sequence_length)
        } else if (identical(result$reason, "duplicate_accession")) {
          counts$dup_acc <- counts$dup_acc + 1L
        } else if (startsWith(result$reason, "duplicate_exact_sequence")) {
          counts$dup_seq <- counts$dup_seq + 1L
        }
      }
    }
    
    all_stats[[length(all_stats) + 1L]] <- data.frame(
      serotype = serotype, country = target_country,
      total_records_found        = length(ids),
      total_valid_sequences      = counts$accepted,
      total_complete_genomes     = counts$complete_genome,
      total_complete_cds         = counts$complete_cds,
      duplicate_accessions       = counts$dup_acc,
      duplicate_exact_sequences  = counts$dup_seq,
      min_length  = if (length(counts$lengths)) min(counts$lengths,  na.rm = TRUE) else NA_real_,
      max_length  = if (length(counts$lengths)) max(counts$lengths,  na.rm = TRUE) else NA_real_,
      mean_length = if (length(counts$lengths)) mean(counts$lengths, na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE)
    log$info("Query done | accepted=%d | dup_acc=%d | dup_seq=%d",
             counts$accepted, counts$dup_acc, counts$dup_seq)
  }
}

if (query_failures == nrow(countries_df) * length(serotypes)) {
  stop("All ESearch queries failed; refusing to overwrite final outputs.")
}

# ---- assemble outputs -----------------------------------------------------
empty_metadata <- function() {
  cols <- c("accession", "serotype", "query_country", "country", "region",
            "locality", "geo_loc_name", "country_field", "organism",
            "collection_date", "year", "host", "isolate", "strain",
            "sequence_length", "molecule_type", "completeness",
            "genome_or_cds", "definition", "source_line", "journal",
            "authors", "submission_date", "genbank_relative_path",
            "fasta_relative_path", "country_genbank_relative_path",
            "country_fasta_relative_path", "source_database", "notes")
  df <- as.data.frame(matrix(NA_character_, nrow = 0L, ncol = length(cols),
                             dimnames = list(NULL, cols)),
                      stringsAsFactors = FALSE)
  df
}

metadata_df <- if (length(all_metadata) > 0L) {
  do.call(rbind, all_metadata)
} else empty_metadata()

if (nrow(metadata_df) > 0L) {
  metadata_df <- metadata_df[order(metadata_df$country, metadata_df$serotype,
                                   metadata_df$year, metadata_df$accession,
                                   na.last = TRUE), , drop = FALSE]
  rownames(metadata_df) <- NULL
}

stats_df <- if (length(all_stats) > 0L) do.call(rbind, all_stats) else data.frame()
filters_df <- if (length(all_filters) > 0L) do.call(rbind, all_filters) else data.frame()

# Per-serotype globals
stats_by_serotype <- if (nrow(stats_df) > 0L) {
  do.call(rbind, lapply(split(stats_df, stats_df$serotype), function(df) {
    data.frame(
      serotype = df$serotype[[1L]],
      total_valid_sequences   = sum(df$total_valid_sequences, na.rm = TRUE),
      total_complete_genomes  = sum(df$total_complete_genomes, na.rm = TRUE),
      total_complete_cds      = sum(df$total_complete_cds, na.rm = TRUE),
      countries_with_records  = sum(df$total_valid_sequences > 0L, na.rm = TRUE),
      countries_without_records = sum(df$total_valid_sequences == 0L, na.rm = TRUE),
      duplicate_accessions    = sum(df$duplicate_accessions, na.rm = TRUE),
      duplicate_exact_sequences = sum(df$duplicate_exact_sequences, na.rm = TRUE),
      stringsAsFactors = FALSE)
  }))
} else data.frame()

# Per-country totals
country_totals <- if (nrow(metadata_df) > 0L) {
  out <- as.data.frame(table(metadata_df$country), stringsAsFactors = FALSE)
  names(out) <- c("country", "total_sequences")
  out[order(out$country), , drop = FALSE]
} else data.frame(country = countries_df$country, total_sequences = 0L,
                  stringsAsFactors = FALSE)

# Missing countries per serotype
missing_lines <- vapply(names(serotypes), function(s) {
  zero <- if (nrow(stats_df) > 0L) {
    stats_df$country[stats_df$serotype == s & stats_df$total_valid_sequences == 0L]
  } else countries_df$country
  if (length(zero) == 0L) sprintf("%s: none", s)
  else sprintf("%s: %s", s, paste(zero, collapse = ", "))
}, character(1L))

# ---- write summary files --------------------------------------------------
write.csv(metadata_df, file.path(cfg$resumen_dir,
                                 "metadata_DENV_south_america_up_to_2026.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(stats_df, file.path(cfg$resumen_dir,
                              "genome_statistics_DENV_south_america.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
if (nrow(filters_df) > 0L) {
  write.csv(filters_df, file.path(cfg$resumen_dir, "records_filter_log.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
}
writeLines(missing_lines, file.path(cfg$resumen_dir,
                                    "countries_without_records.txt"), useBytes = TRUE)

write_excel(file.path(cfg$resumen_dir, "metadata_DENV_south_america_up_to_2026.xlsx"),
            list(metadata = metadata_df))
write_excel(file.path(cfg$resumen_dir, "genome_statistics_DENV_south_america.xlsx"),
            list(by_serotype_country  = stats_df,
                 global_serotype    = stats_by_serotype,
                 global_country        = country_totals))

# Per-country metadata sidecars
for (i in seq_len(nrow(countries_df))) {
  ctry <- countries_df$country[[i]]
  for (s in names(serotypes)) {
    sub <- metadata_df[metadata_df$country == ctry & metadata_df$serotype == s, , drop = FALSE]
    out <- file.path(cfg$country_view, ctry, s, "metadata", sprintf("metadata_%s_%s.csv", ctry, s))
    ensure_dir(dirname(out))
    write.csv(sub, out, row.names = FALSE, fileEncoding = "UTF-8")
  }
}
# Per-serotype metadata sidecars
for (s in names(serotypes)) {
  sub <- metadata_df[metadata_df$serotype == s, , drop = FALSE]
  out <- file.path(project_root, "results", s, "metadata", sprintf("metadata_%s.csv", s))
  ensure_dir(dirname(out))
  write.csv(sub, out, row.names = FALSE, fileEncoding = "UTF-8")
}

# Consolidated FASTA per serotype
combined_dir <- file.path(cfg$resumen_dir, "consolidated_fasta_by_serotype")
for (s in names(serotypes)) {
  sub <- metadata_df[metadata_df$serotype == s, , drop = FALSE]
  out <- file.path(combined_dir,
                   sprintf("%s_south_america_up_to_%d.fasta", s, cfg$max_year))
  if (nrow(sub) == 0L) {
    writeLines(sprintf("; No valid sequences retained for %s", s), out)
    next
  }
  sub <- sub[order(sub$country, sub$year, sub$accession, na.last = TRUE), ,
             drop = FALSE]
  blocks <- vapply(seq_len(nrow(sub)), function(j) {
    fasta <- read_fasta(file.path(project_root, sub$fasta_relative_path[[j]]))
    header <- sprintf("%s|serotype=%s|country=%s|year=%s|host=%s|collection_date=%s",
                      clean_header_value(sub$accession[[j]]),
                      clean_header_value(sub$serotype[[j]]),
                      clean_header_value(sub$country[[j]]),
                      clean_header_value(sub$year[[j]]),
                      clean_header_value(sub$host[[j]]),
                      clean_header_value(sub$collection_date[[j]]))
    fasta_text(fasta$sequence, header)
  }, character(1L))
  writeLines(blocks, out, useBytes = TRUE)
}

# ---- verify ---------------------------------------------------------------
if (nrow(metadata_df) > 0L && nrow(stats_df) > 0L) {
  counts_meta <- as.data.frame(table(metadata_df$serotype, metadata_df$country),
                               stringsAsFactors = FALSE)
  names(counts_meta) <- c("serotype", "country", "metadata_n")
  merged <- merge(stats_df[, c("serotype", "country", "total_valid_sequences")],
                  counts_meta, by = c("serotype", "country"), all.x = TRUE)
  merged$metadata_n[is.na(merged$metadata_n)] <- 0L
  if (!isTRUE(all(merged$total_valid_sequences == merged$metadata_n))) {
    write.csv(merged, file.path(cfg$resumen_dir,
                                "verificacion_estadisticas_vs_metadata.csv"),
              row.names = FALSE)
    stop("Statistics verification failed: totals do not match metadata.")
  }
}

log$info("Done. Metadata rows: %d | Total accepted: %d",
         nrow(metadata_df), sum(stats_df$total_valid_sequences))

cat("\nTotal sequences per serotype:\n")
print(stats_by_serotype, row.names = FALSE)
cat("\nTotal sequences per country:\n")
print(country_totals, row.names = FALSE)
cat("\nConsolidated FASTA files:\n")
for (s in names(serotypes)) {
  cat(sprintf("  - %s\n", file.path(combined_dir,
                                    sprintf("%s_suramerica_hasta_%d.fasta", s, cfg$max_year))))
}