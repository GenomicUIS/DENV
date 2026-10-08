# NCBI E-utilities client built on httr2.

new_ncbi_client <- function(api_key = trimws("XXXXXXXXXXXXXXXXXXXXXXXXXX"),
                            email   = trimws("XXXXXXXXXXXXXXXXXXXXXXXXXX"),
                            tool    = Sys.getenv("NCBI_TOOL",
                                                 "DENGUE_SURAMERICA"),
                            base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils",
                            timeout = 60L,
                            max_tries = 5L,
                            rate_limit = if (nzchar(api_key)) 10 else 3) {
  if (!requireNamespace("httr2", quietly = TRUE)) {
    stop("Package 'httr2' is required. Run scripts/00_setup.R first.")
  }
  structure(
    list(
      api_key = api_key,
      email   = email,
      tool    = tool,
      base_url = base_url,
      timeout = timeout,
      max_tries = max_tries,
      rate_limit = rate_limit
    ),
    class = "ncbi_client"
  )
}

print.ncbi_client <- function(x, ...) {
  cat("<ncbi_client>\n")
  cat("  base_url:   ", x$base_url, "\n", sep = "")
  cat("  tool:       ", x$tool, "\n", sep = "")
  cat("  email set:  ", nzchar(x$email), "\n", sep = "")
  cat("  api_key set:", nzchar(x$api_key), "\n", sep = "")
  cat("  rate_limit: ", x$rate_limit, " req/s\n", sep = "")
  invisible(x)
}

ncbi_build_request <- function(client, endpoint, params = list()) {
  req <- httr2::request(client$base_url)
  req <- httr2::req_url_path_append(req, endpoint)
  req <- do.call(httr2::req_url_query, c(list(req), params))
  req <- httr2::req_user_agent(req, client$tool)
  req <- httr2::req_timeout(req, client$timeout)
  req <- httr2::req_retry(
    req,
    max_tries = client$max_tries,
    is_transient = function(resp) {
      httr2::resp_status(resp) %in% c(429L, 500L, 502L, 503L, 504L)
    }
  )
  req <- httr2::req_throttle(req, rate = client$rate_limit, realm = "ncbi-eutils")
  if (nzchar(client$email))   req <- httr2::req_url_query(req, email   = client$email)
  if (nzchar(client$tool))    req <- httr2::req_url_query(req, tool    = client$tool)
  if (nzchar(client$api_key)) req <- httr2::req_url_query(req, api_key = client$api_key)
  req
}

ncbi_perform_text <- function(req, api_key = "", email = "") {
  resp <- tryCatch(httr2::req_perform(req),
                   error = function(e) {
                     msg <- redact_secrets(conditionMessage(e), api_key, email)
                     stop(msg, call. = FALSE)
                   })
  httr2::resp_check_status(resp)
  httr2::resp_body_string(resp)
}

#' Single ESearch page. Returns list(count, ids).
esearch_page <- function(client, query, retstart = 0L, retmax = 500L) {
  req <- ncbi_build_request(client, "esearch.fcgi", list(
    db = "nucleotide",
    term = query,
    retstart = retstart,
    retmax = retmax,
    retmode = "xml"
  ))
  body <- ncbi_perform_text(req, client$api_key, client$email)
  doc <- xml2::read_xml(body)
  count <- as.integer(xml2::xml_text(xml2::xml_find_first(doc, ".//Count")))
  ids <- xml2::xml_text(xml2::xml_find_all(doc, ".//IdList/Id"))
  ids <- ids[nzchar(ids)]
  list(count = if (is.na(count)) 0L else count, ids = ids)
}

#' Paginated ESearch; returns a de-duplicated character vector of IDs.
esearch_all_ids <- function(client, query, retmax = 500L, max_ids = Inf) {
  first <- esearch_page(client, query, retstart = 0L, retmax = retmax)
  if (first$count == 0L) return(character(0))
  total <- min(first$count, max_ids)
  ids <- first$ids
  if (length(ids) < total) {
    offsets <- seq(retmax, total - 1L, by = retmax)
    for (off in offsets) {
      page <- esearch_page(client, query, retstart = off, retmax = retmax)
      ids <- c(ids, page$ids)
      if (length(ids) >= total) break
    }
  }
  ids <- unique(ids)
  ids[seq_len(min(length(ids), total))]
}

#' Batched EFetch of GenBank records. Returns a list of record strings.
efetch_genbank <- function(client, ids, batch_size = 100L) {
  ids <- as.character(ids)
  if (length(ids) == 0L) return(list())
  batches <- split(ids, ceiling(seq_along(ids) / batch_size))
  out <- list()
  for (batch in batches) {
    req <- ncbi_build_request(client, "efetch.fcgi", list(
      db = "nucleotide",
      id = paste(batch, collapse = ","),
      rettype = "gbwithparts",
      retmode = "text"
    ))
    body <- ncbi_perform_text(req, client$api_key, client$email)
    out <- c(out, split_genbank_records(body))
  }
  out
}

#' Split a concatenated GenBank response into individual records.
split_genbank_records <- function(text) {
  if (!nzchar(text)) return(list())
  lines <- strsplit(text, "\n", fixed = TRUE)[[1L]]
  end_idx <- which(lines == "//")
  if (length(end_idx) == 0L) return(list(text))
  starts <- c(1L, end_idx[-length(end_idx)] + 1L)
  mapply(function(s, e) paste(lines[s:e], collapse = "\n"),
         starts, end_idx, SIMPLIFY = FALSE, USE.NAMES = FALSE)
}
