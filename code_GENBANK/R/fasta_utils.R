# FASTA reading, writing, and record splitting.

#' Read a FASTA file into a list(header, sequence).
read_fasta <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  if (length(lines) == 0L) {
    return(list(header = NA_character_, sequence = ""))
  }
  header_idx <- which(startsWith(lines, ">"))
  if (length(header_idx) == 0L) {
    return(list(header = NA_character_, sequence = ""))
  }
  header <- sub("^>", "", lines[header_idx[[1L]]])
  sequence <- toupper(gsub("\\s+", "", paste(lines[-header_idx], collapse = "")))
  list(header = header, sequence = sequence)
}

#' Wrap a sequence to a fixed line width.
wrap_sequence <- function(sequence, width = 80L) {
  if (!nzchar(sequence)) return(character(0))
  n <- nchar(sequence)
  starts <- seq(1L, n, by = width)
  vapply(starts, function(s) substr(sequence, s, min(s + width - 1L, n)),
         character(1L))
}

#' Write a FASTA record to disk.
write_fasta <- function(sequence, header, path, width = 80L) {
  ensure_dir(dirname(path))
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con))
  writeLines(paste0(">", header), con)
  body <- wrap_sequence(sequence, width = width)
  if (length(body) > 0L) writeLines(body, con)
  invisible(path)
}

#' Build a FASTA text block (header + wrapped body) without touching disk.
fasta_text <- function(sequence, header, width = 80L) {
  body <- wrap_sequence(sequence, width = width)
  paste0(">", header, "\n", paste(body, collapse = "\n"), "\n")
}