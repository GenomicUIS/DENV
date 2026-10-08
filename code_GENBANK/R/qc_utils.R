# Sequence QC, taxonomy checks, and country normalization.

sequence_qc <- function(sequence) {
  chars <- strsplit(toupper(sequence), "", fixed = TRUE)[[1L]]
  if (length(chars) == 1L && !nzchar(chars)) chars <- character(0)
  n_total <- length(chars)
  count_chr <- function(letters) sum(chars %in% letters)
  a <- count_chr("A"); c_ <- count_chr("C")
  g <- count_chr("G"); t_ <- count_chr("T")
  n <- count_chr("N"); gaps <- count_chr("-")
  iupac <- count_chr(c("R","Y","S","W","K","M","B","D","H","V"))
  other <- n_total - (a + c_ + g + t_ + n + gaps + iupac)
  non_acgt <- n + gaps + iupac + other
  
  list(
    coding_length = n_total,
    A = a, C = c_, G = g, T = t_,
    N = n, gaps = gaps,
    iupac_ambiguous = iupac,
    other_non_acgt = other,
    non_acgt_total = non_acgt,
    pct_N = if (n_total > 0L) 100 * n / n_total else NA_real_,
    pct_gaps = if (n_total > 0L) 100 * gaps / n_total else NA_real_,
    pct_non_acgt = if (n_total > 0L) 100 * non_acgt / n_total else NA_real_
  )
}

classify_completeness <- function(definition) {
  d <- tolower(definition %||% "")
  if (grepl("complete genome", d, fixed = TRUE)) {
    return(list(completeness = "complete_genome", genome_or_cds = "genome"))
  }
  if (grepl("complete cds", d, fixed = TRUE) ||
      grepl("complete coding sequence", d, fixed = TRUE)) {
    return(list(completeness = "complete_cds", genome_or_cds = "cds"))
  }
  list(completeness = "other", genome_or_cds = "other")
}

is_dengue_record <- function(definition, organism, source_line = "") {
  text <- tolower(paste(str_or(definition, ""),
                        str_or(organism, ""),
                        str_or(source_line, "")))
  grepl("dengue|denv|orthoflavivirus denguei", text, perl = TRUE)
}

parse_serotype <- function(definition, organism, source_line = "") {
  text <- tolower(paste(str_or(definition, ""),
                        str_or(organism, ""),
                        str_or(source_line, "")))
  matches <- character(0)
  if (grepl("type\\s*1|denv[- ]?1|dengue virus 1", text, perl = TRUE)) matches <- c(matches, "DENV1")
  if (grepl("type\\s*2|denv[- ]?2|dengue virus 2", text, perl = TRUE)) matches <- c(matches, "DENV2")
  if (grepl("type\\s*3|denv[- ]?3|dengue virus 3", text, perl = TRUE)) matches <- c(matches, "DENV3")
  if (grepl("type\\s*4|denv[- ]?4|dengue virus 4", text, perl = TRUE)) matches <- c(matches, "DENV4")
  matches <- unique(matches)
  if (length(matches) == 1L) matches[[1L]] else NA_character_
}

normalize_country <- function(x) {
  x <- trimws(tolower(str_or(x, "")))
  x <- sub(":.*$", "", x)
  x <- gsub("\\s+", " ", x)
  map <- c(
    "argentina" = "Argentina",
    "bolivia" = "Bolivia",
    "bolivia, plurinational state of" = "Bolivia",
    "brazil" = "Brazil",
    "brasil" = "Brazil",
    "chile" = "Chile",
    "colombia" = "Colombia",
    "ecuador" = "Ecuador",
    "guyana" = "Guyana",
    "paraguay" = "Paraguay",
    "peru" = "Peru",
    "suriname" = "Suriname",
    "surinam"  = "Suriname",
    "uruguay" = "Uruguay",
    "venezuela" = "Venezuela",
    "venezuela, bolivarian republic of" = "Venezuela"
  )
  out <- unname(map[x])
  ifelse(is.na(out), tools::toTitleCase(x), out)
}

#' Turn a parsed GenBank record into a flat metadata row.
genbank_to_metadata_row <- function(record, serotype, target_country) {
  classification <- classify_completeness(record$definition)
  geo_parts <- strsplit(record$geo_loc_name %||% "", ":", fixed = TRUE)[[1L]]
  country_raw <- if (length(geo_parts) >= 1L) geo_parts[[1L]] else ""
  rest <- if (length(geo_parts) >= 2L) paste(geo_parts[-1L], collapse = ":") else ""
  rest_parts <- trimws(strsplit(rest, ",", fixed = TRUE)[[1L]])
  rest_parts <- rest_parts[nzchar(rest_parts)]
  region <- if (length(rest_parts) >= 1L) rest_parts[[1L]] else ""
  locality <- if (length(rest_parts) >= 2L) paste(rest_parts[-1L], collapse = ", ") else ""
  
  data.frame(
    accession        = record$accession,
    serotype         = serotype,
    query_country    = target_country,
    country          = normalize_country(country_raw),
    region           = region,
    locality         = locality,
    geo_loc_name     = record$geo_loc_name  %||% NA_character_,
    country_field    = record$country_field %||% NA_character_,
    organism         = record$organism      %||% NA_character_,
    collection_date  = record$collection_date %||% NA_character_,
    year             = record$year,
    host             = record$host          %||% NA_character_,
    isolate          = record$isolate       %||% NA_character_,
    strain           = record$strain        %||% NA_character_,
    sequence_length  = record$sequence_length,
    molecule_type    = record$molecule_type %||% record$mol_type %||% NA_character_,
    completeness     = classification$completeness,
    genome_or_cds    = classification$genome_or_cds,
    definition       = record$definition    %||% NA_character_,
    source_line      = record$source        %||% NA_character_,
    journal          = record$journal       %||% NA_character_,
    authors          = record$authors       %||% NA_character_,
    submission_date  = record$submission_date %||% NA_character_,
    stringsAsFactors = FALSE
  )
}

#' Header-level exclusion rules applied before accepting a record.
#' Returns NA when the record passes, otherwise a reason string.
header_exclusion_reason <- function(row, max_year) {
  definition <- tolower(row$definition %||% "")
  header_text <- tolower(paste(row$definition %||% "",
                               row$source_line %||% "",
                               row$organism %||% ""))
  if (!is.na(row$year) && row$year > max_year) {
    return(sprintf("collection_year_after_%d", max_year))
  }
  if (row$completeness == "other") {
    return("not_complete_genome_or_cds")
  }
  if (grepl("\\bpartial\\b|\\bfragment\\b|\\bincomplete\\b|\\bamplicon\\b",
            definition, perl = TRUE)) {
    return("excluded_incomplete_definition")
  }
  if (grepl("\\bsynthetic\\b|\\bartificial\\b|\\bconstruct\\b",
            header_text, perl = TRUE)) {
    return("excluded_synthetic_or_construct")
  }
  if (grepl("\\bpatent\\b", tolower(row$journal %||% ""), perl = TRUE)) {
    return("excluded_patent")
  }
  NA_character_
}