# Line-based GenBank flat-file parser.

#' Parse a single GenBank record into a flat list.
parse_genbank_record <- function(text) {
  lines <- strsplit(text, "\n", fixed = TRUE)[[1L]]
  lines <- sub("\r$", "", lines)
  
  locus_line <- grep("^LOCUS", lines, value = TRUE)
  locus_line <- if (length(locus_line) > 0L) locus_line[[1L]] else ""
  tokens <- strsplit(trimws(locus_line), "\\s+")[[1L]]
  
  definition <- collect_multiline(lines, "^DEFINITION")
  accession  <- first_value(lines, "^ACCESSION")
  version    <- first_value(lines, "^VERSION")
  version    <- if (!is.na(version)) sub("\\s.*$", "", version) else NA_character_
  source_ln  <- first_value(lines, "^SOURCE")
  organism   <- first_value(lines, "^\\s+ORGANISM")
  
  # First REFERENCE block only: journal and authors.
  ref_idx <- grep("^REFERENCE", lines)
  journal <- NA_character_
  authors <- NA_character_
  if (length(ref_idx) > 0L) {
    end <- if (length(ref_idx) > 1L) ref_idx[[2L]] - 1L else length(lines)
    block <- lines[ref_idx[[1L]]:end]
    jl <- grep("^\\s+JOURNAL\\s+", block, value = TRUE)
    al <- grep("^\\s+AUTHORS\\s+", block, value = TRUE)
    if (length(jl) > 0L) journal <- trimws(gsub("^\\s+JOURNAL\\s+", "",
                                                paste(jl, collapse = " ")))
    if (length(al) > 0L) authors <- trimws(gsub("^\\s+AUTHORS\\s+", "",
                                                paste(al, collapse = " ")))
  }
  
  features <- parse_genbank_features(lines)
  source_feat <- NULL
  for (f in features) {
    if (identical(f$key, "source")) { source_feat <- f; break }
  }
  qual <- if (!is.null(source_feat)) source_feat$qualifiers else list()
  
  collection_date <- qual$collection_date %||% NA_character_
  year <- extract_year(collection_date)
  if (is.na(year)) year <- extract_year(definition)
  
  list(
    locus_name      = if (length(tokens) >= 2L) tokens[[2L]] else NA_character_,
    accession       = if (!is.na(version)) version else accession,
    accession_base  = accession,
    sequence_length = if (length(tokens) >= 3L)
      suppressWarnings(as.integer(tokens[[3L]]))
    else NA_integer_,
    molecule_type   = if (length(tokens) >= 5L) tokens[[5L]] else NA_character_,
    topology        = if (length(tokens) >= 6L) tokens[[6L]] else NA_character_,
    division        = if (length(tokens) >= 7L) tokens[[7L]] else NA_character_,
    submission_date = if (length(tokens) >= 8L) tokens[[8L]] else NA_character_,
    definition      = definition,
    source          = source_ln,
    organism        = organism,
    journal         = journal,
    authors         = authors,
    geo_loc_name    = qual$geo_loc_name    %||% NA_character_,
    country_field   = qual$country         %||% NA_character_,
    collection_date = collection_date,
    year            = year,
    host            = qual$host            %||% NA_character_,
    isolate         = qual$isolate         %||% NA_character_,
    strain          = qual$strain          %||% NA_character_,
    mol_type        = qual$mol_type        %||% NA_character_,
    features        = features,
    sequence        = extract_genbank_origin(lines)
  )
}

first_value <- function(lines, pattern) {
  idx <- grep(pattern, lines)
  if (length(idx) == 0L) return(NA_character_)
  trimws(sub(paste0(pattern, "\\s*"), "", lines[[idx[[1L]]]]))
}

collect_multiline <- function(lines, pattern) {
  start <- grep(pattern, lines)
  if (length(start) == 0L) return(NA_character_)
  start <- start[[1L]]
  value <- sub(paste0(pattern, "\\s*"), "", lines[[start]])
  i <- start + 1L
  while (i <= length(lines)) {
    line <- lines[[i]]
    if (grepl("^[A-Z]", line) && !grepl("^\\s", line)) break
    if (nzchar(trimws(line))) value <- paste(value, trimws(line))
    i <- i + 1L
  }
  trimws(value)
}

extract_year <- function(text) {
  if (is.null(text) || is.na(text)) return(NA_integer_)
  m <- regmatches(text, regexpr("(19|20)[0-9]{2}", text))
  if (length(m) == 0L || !nzchar(m)) return(NA_integer_)
  as.integer(m)
}

#' Parse the FEATURES section into a list of features.
#' Each feature: list(key, location, qualifiers = list(), last_qual).
parse_genbank_features <- function(lines) {
  feat_start <- grep("^FEATURES", lines)
  if (length(feat_start) == 0L) return(list())
  origin_start <- grep("^ORIGIN", lines)
  end <- if (length(origin_start) > 0L) origin_start[[1L]] - 1L else length(lines)
  if (end < feat_start[[1L]] + 1L) return(list())
  
  feat_lines <- lines[(feat_start[[1L]] + 1L):end]
  features <- list()
  current <- NULL
  
  for (raw in feat_lines) {
    line <- sub("\r$", "", raw)
    if (!nzchar(trimws(line))) next
    
    if (grepl("^     \\S", line)) {
      # New feature: "     CDS             join(97..10317)"
      if (!is.null(current)) features[[length(features) + 1L]] <- current
      key <- sub("^\\s+(\\S+).*$", "\\1", line)
      loc <- trimws(sub("^\\s+\\S+\\s+", "", line))
      current <- list(key = key, location = loc,
                      qualifiers = list(), last_qual = NA_character_)
    } else if (!is.null(current) && grepl("^\\s{6,}", line)) {
      qm <- regmatches(line, regexec(
        "^\\s+/([A-Za-z_][A-Za-z0-9_]*)(?:=(.*))?$", line))[[1L]]
      if (length(qm) >= 2L) {
        qname <- qm[[2L]]
        qval <- if (length(qm) >= 3L && !is.na(qm[[3L]])) strip_quotes(qm[[3L]]) else ""
        current$qualifiers[[qname]] <- qval
        current$last_qual <- qname
      } else if (is.na(current$last_qual)) {
        current$location <- paste0(current$location, trimws(line))
      } else {
        cont <- strip_quotes(trimws(line))
        current$qualifiers[[current$last_qual]] <-
          paste(current$qualifiers[[current$last_qual]], cont)
      }
    }
  }
  if (!is.null(current)) features[[length(features) + 1L]] <- current
  features
}

strip_quotes <- function(x) {
  x <- trimws(x)
  if (grepl('^".*"$', x)) return(substring(x, 2L, nchar(x) - 1L))
  if (startsWith(x, '"')) return(substring(x, 2L))
  if (endsWith(x, '"'))   return(substring(x, 1L, nchar(x) - 1L))
  x
}

extract_genbank_origin <- function(lines) {
  origin_idx <- grep("^ORIGIN", lines)
  if (length(origin_idx) == 0L) return("")
  end_idx <- grep("^//", lines)
  end_idx <- end_idx[end_idx > origin_idx[[1L]]]
  if (length(end_idx) == 0L) end_idx <- length(lines) + 1L
  origin_lines <- lines[(origin_idx[[1L]] + 1L):(end_idx[[1L]] - 1L)]
  seq <- paste(origin_lines, collapse = "")
  seq <- gsub("[^A-Za-z.-]", "", seq)
  seq <- gsub("\\.", "-", seq)
  toupper(seq)
}

#' Parse a GenBank location string into a list(ok, ranges, complement, fuzzy).
parse_genbank_location <- function(location) {
  loc <- gsub("\\s+", "", location %||% "")
  if (!nzchar(loc)) {
    return(list(ok = FALSE, ranges = NULL, complement = FALSE, fuzzy = FALSE))
  }
  complement <- FALSE
  fuzzy <- grepl("[<>]", loc)
  
  if (grepl("^complement\\(", loc)) {
    complement <- TRUE
    loc <- sub("^complement\\((.*)\\)$", "\\1", loc)
  }
  if (grepl("^(join|order)\\(", loc)) {
    loc <- sub("^(join|order)\\((.*)\\)$", "\\2", loc)
  }
  loc <- gsub("[<>]", "", loc)
  parts <- strsplit(loc, ",", fixed = TRUE)[[1L]]
  
  ranges <- lapply(parts, function(p) {
    if (grepl("^[0-9]+\\.\\.[0-9]+$", p)) {
      coords <- as.integer(strsplit(p, "..", fixed = TRUE)[[1L]])
      data.frame(start = coords[[1L]], end = coords[[2L]])
    } else if (grepl("^[0-9]+$", p)) {
      pos <- as.integer(p)
      data.frame(start = pos, end = pos)
    } else NULL
  })
  if (any(vapply(ranges, is.null, logical(1L)))) {
    return(list(ok = FALSE, ranges = NULL,
                complement = complement, fuzzy = fuzzy))
  }
  ranges <- do.call(rbind, ranges)
  list(
    ok = TRUE,
    ranges = ranges,
    complement = complement,
    fuzzy = fuzzy,
    start = min(ranges$start),
    end   = max(ranges$end)
  )
}

extract_from_location <- function(sequence, parsed) {
  n <- nchar(sequence)
  if (any(parsed$ranges$start < 1L) ||
      any(parsed$ranges$end > n) ||
      any(parsed$ranges$end < parsed$ranges$start)) {
    return(list(ok = FALSE, sequence = "", reason = "coords_out_of_bounds"))
  }
  pieces <- mapply(function(s, e) substr(sequence, s, e),
                   parsed$ranges$start, parsed$ranges$end,
                   SIMPLIFY = TRUE, USE.NAMES = FALSE)
  out <- paste(pieces, collapse = "")
  if (parsed$complement) out <- reverse_complement(out)
  list(ok = TRUE, sequence = out, reason = "ok")
}

reverse_complement <- function(sequence) {
  chars <- strsplit(toupper(sequence), "", fixed = TRUE)[[1L]]
  lookup <- c(A = "T", C = "G", G = "C", T = "A",
              R = "Y", Y = "R", S = "S", W = "W",
              K = "M", M = "K", B = "V", D = "H",
              H = "D", V = "B", N = "N", "-" = "-")
  mapped <- lookup[chars]
  mapped[is.na(mapped)] <- "N"
  paste(rev(mapped), collapse = "")
}

#' Pick a CDS from a parsed record and return its sequence.
#' Prefers polyprotein products, then the longest CDS.
extract_cds_sequence <- function(record, fallback_sequence = NULL) {
  cds_features <- Filter(function(f) identical(f$key, "CDS"), record$features)
  if (length(cds_features) == 0L) {
    return(list(status = "no_cds_found",
                sequence = fallback_sequence %||% "",
                location = NA_character_, product = NA_character_,
                note = NA_character_, protein_id = NA_character_,
                start = NA_integer_, end = NA_integer_, cds_length = NA_integer_,
                fuzzy = FALSE, source = "fallback"))
  }
  origin <- if (nzchar(record$sequence)) record$sequence else fallback_sequence
  if (is.null(origin) || !nzchar(origin)) {
    return(list(status = "no_sequence_available",
                sequence = "", location = NA_character_,
                product = NA_character_, note = NA_character_,
                protein_id = NA_character_, start = NA_integer_, end = NA_integer_,
                cds_length = NA_integer_, fuzzy = FALSE, source = "none"))
  }
  
  candidates <- lapply(cds_features, function(f) {
    parsed <- parse_genbank_location(f$location)
    if (!isTRUE(parsed$ok)) return(list(ok = FALSE))
    ext <- extract_from_location(origin, parsed)
    list(ok = ext$ok, feature = f, parsed = parsed,
         sequence = ext$sequence, reason = ext$reason)
  })
  valid <- Filter(function(c) isTRUE(c$ok), candidates)
  if (length(valid) == 0L) {
    return(list(status = "cds_parse_failed",
                sequence = fallback_sequence %||% "",
                location = NA_character_, product = NA_character_,
                note = NA_character_, protein_id = NA_character_,
                start = NA_integer_, end = NA_integer_, cds_length = NA_integer_,
                fuzzy = FALSE, source = "fallback"))
  }
  poly <- vapply(valid, function(c) {
    grepl("polyprotein",
          paste(c$feature$qualifiers$product %||% "",
                c$feature$qualifiers$note %||% ""),
          ignore.case = TRUE)
  }, logical(1L))
  lengths <- vapply(valid, function(c) nchar(c$sequence), integer(1L))
  ord <- order(-as.integer(poly), -lengths)
  best <- valid[[ord[[1L]]]]
  
  list(status = "cds_extracted",
       sequence = best$sequence,
       location = best$feature$location,
       product  = best$feature$qualifiers$product    %||% NA_character_,
       note     = best$feature$qualifiers$note       %||% NA_character_,
       protein_id = best$feature$qualifiers$protein_id %||% NA_character_,
       start = best$parsed$start,
       end   = best$parsed$end,
       cds_length = nchar(best$sequence),
       fuzzy = best$parsed$fuzzy,
       source = "genbank_cds")
}