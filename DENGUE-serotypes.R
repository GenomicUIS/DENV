###############################################################################
# Download of complete genomic sequences of Dengue virus serotype 1.

# Required libraries
library(httr)
library(jsonlite)

# Constants and configuration
BASE_URL <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
API_KEY <- "0142866f7ab066271cca2802824901de9408"
BATCH_SIZE <- 100
RETRY_ATTEMPTS <- 5
RETRY_DELAY <- 3
OUTPUT_DIR <- "C:/Users/caden/Documents/dengue_genomes"

# Create output directory
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# Function to perform safe HTTP requests
safe_http_get <- function(url) {
  for (attempt in seq_len(RETRY_ATTEMPTS)) {
    response <- tryCatch({
      GET(url)
    }, error = function(e) NULL)
    if (!is.null(response) && status_code(response) == 200) {
      return(content(response, as = "text", encoding = "UTF-8"))
    }
    Sys.sleep(RETRY_DELAY * attempt)
  }
  return(NULL)
}

# Function to search for genome IDs
search_genomes <- function(retstart, query) {
  url <- sprintf(
    "%s/esearch.fcgi?db=nucleotide&term=%s&retstart=%d&api_key=%s&retmax=%d&retmode=json",
    BASE_URL, URLencode(query), retstart, API_KEY, BATCH_SIZE
  )
  response <- safe_http_get(url)
  if (is.null(response)) return(character(0))
  result <- fromJSON(response, simplifyVector = TRUE)
  return(result$esearchresult$idlist)
}

# Function to download a genome
fetch_genome <- function(id) {
  url <- sprintf(
    "%s/efetch.fcgi?db=nucleotide&id=%s&rettype=fasta&retmode=text&api_key=%s",
    BASE_URL, id, API_KEY
  )
  safe_http_get(url)
}

# Function to validate sequences
is_valid_sequence <- function(sequence) {
  grepl("Dengue virus type 1|Dengue virus 1|DENV-1|Orthoflavivirus denguei", sequence) &&
    grepl("complete genome|complete cds", sequence) &&
    !grepl("partial", sequence)
}

# Function to process a genome
process_genome <- function(id) {
  sequence <- fetch_genome(id)
  if (is.null(sequence)) return(list(included = 0, excluded = 0, errors = 1))
  if (is_valid_sequence(sequence)) {
    filepath <- file.path(OUTPUT_DIR, sprintf("dengue_genome_%s.fasta", id))
    writeLines(sequence, filepath)
    return(list(included = 1, excluded = 0, errors = 0))
  }
  return(list(included = 0, excluded = 1, errors = 0))
}

# Main function to download genomes
download_dengue_genomes <- function() {
  # Query with keywords at organism level and sequence characteristics
  query <- paste(
    "(Dengue virus type 1[Organism] OR Dengue virus 1[Organism] OR DENV-1[Organism] OR Orthoflavivirus denguei[Organism])",
    "AND (complete genome[Title] OR complete cds[Title])",
    sep = " "
  )
  
  retstart <- 0
  stats <- list(included = 0, excluded = 0, errors = 0)
  
  repeat {
    message(sprintf("Searching sequences... Retstart: %d", retstart))
    ids <- search_genomes(retstart, query)
    if (length(ids) == 0) break
    
    for (id in ids) {
      result <- process_genome(id)
      stats$included <- stats$included + result$included
      stats$excluded <- stats$excluded + result$excluded
      stats$errors <- stats$errors + result$errors
    }
    
    retstart <- retstart + BATCH_SIZE
    Sys.sleep(1)
  }
  
  message(sprintf(
    "Processing completed: Included %d, Excluded %d, Errors %d",
    stats$included, stats$excluded, stats$errors
  ))
  return(stats)
}

# Execute the main function
download_dengue_genomes()

###############################################################################
# Merge all downloaded sequences into a single FASTA file for serotype 1

# Required libraries
library(tidyverse)

# Directory where FASTA files are stored
input_dir <- "C:\\Users\\caden\\Documents\\dengue_genomes_DENV1"
output_file <- "C:\\Users\\caden\\Documents\\dengue_genomes_DENV1\\dengue_combined.fasta"

# Function to combine FASTA sequences
combine_fasta_files <- function(input_dir, output_file) {
  # Get all FASTA files from the directory
  fasta_files <- list.files(input_dir, pattern = "\\.fasta$", full.names = TRUE)
  
  if (length(fasta_files) == 0) {
    stop("No FASTA files found in the provided directory.")
  }
  
  # Read and merge the contents of FASTA files
  combined_fasta <- fasta_files %>%
    map_chr(~ readLines(.x) %>% paste(collapse = "\n")) %>%
    paste(collapse = "\n")
  
  # Write the combined file
  writeLines(combined_fasta, output_file)
  message(sprintf("Combined file saved at: %s", output_file))
}

# Execute the function
combine_fasta_files(input_dir, output_file)

################################################################################
# Remove gaps (-) and undetermined nucleotides (N) for serotype 1.

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("Biostrings")

# Required libraries
library(Biostrings)

# Input and output files
input_fasta <- "C:\\Users\\caden\\Documents\\dengue_genomes_DENV1\\dengue_combined.FASTA"
output_fasta <- "C:\\Users\\caden\\Documents\\dengue_genomes_DENV1\\dengue_purged.FASTA"

# Function to load and purge FASTA sequences
purge_fasta <- function(input_fasta, output_fasta) {
  # Read the FASTA file
  sequences <- readDNAStringSet(filepath = input_fasta, format = "fasta")
  
  # Initialize counters
  total_sequences <- length(sequences)
  valid_sequences <- 0
  removed_sequences <- 0
  
  # Filter valid sequences
  purged_sequences <- DNAStringSet()
  for (i in seq_along(sequences)) {
    seq <- as.character(sequences[[i]])
    if (!grepl("[-N]", seq)) { # Check for gaps (-) and undetermined nucleotides (N)
      purged_sequences <- append(purged_sequences, sequences[i])
      valid_sequences <- valid_sequences + 1
    } else {
      removed_sequences <- removed_sequences + 1
    }
  }
  
  # Save valid sequences to output file
  writeXStringSet(purged_sequences, filepath = output_fasta, format = "fasta")
  
  # Print statistics
  message(sprintf("Total sequences reviewed: %d", total_sequences))
  message(sprintf("Valid sequences (without gaps or undetermined nucleotides): %d", valid_sequences))
  message(sprintf("Removed sequences: %d", removed_sequences))
}

# Execute the function
purge_fasta(input_fasta, output_fasta)

###############################################################################
# Download of complete genomic sequences of Dengue virus serotype 2.
# Required libraries
library(httr)
library(jsonlite)

# Constants and configuration
BASE_URL <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
API_KEY <- "0142866f7ab066271cca2802824901de9408"
BATCH_SIZE <- 100
RETRY_ATTEMPTS <- 5
RETRY_DELAY <- 3
OUTPUT_DIR <- "C:/Users/caden/Documents/dengue_genomes_DENV2"

# Create output directory
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# Function to perform secure HTTP requests
safe_http_get <- function(url) {
  for (attempt in seq_len(RETRY_ATTEMPTS)) {
    response <- tryCatch({
      GET(url)
    }, error = function(e) NULL)
    if (!is.null(response) && status_code(response) == 200) {
      return(content(response, as = "text", encoding = "UTF-8"))
    }
    Sys.sleep(RETRY_DELAY * attempt)
  }
  return(NULL)
}

# Function to search for genome IDs
search_genomes <- function(retstart, query) {
  url <- sprintf(
    "%s/esearch.fcgi?db=nucleotide&term=%s&retstart=%d&api_key=%s&retmax=%d&retmode=json",
    BASE_URL, URLencode(query), retstart, API_KEY, BATCH_SIZE
  )
  response <- safe_http_get(url)
  if (is.null(response)) return(character(0))
  result <- fromJSON(response, simplifyVector = TRUE)
  return(result$esearchresult$idlist)
}

# Function to download a genome
fetch_genome <- function(id) {
  url <- sprintf(
    "%s/efetch.fcgi?db=nucleotide&id=%s&rettype=fasta&retmode=text&api_key=%s",
    BASE_URL, id, API_KEY
  )
  safe_http_get(url)
}

# Function to validate sequences
is_valid_sequence <- function(sequence) {
  grepl("Dengue virus type 2|Dengue virus 2|DENV-2|DENV2|Orthoflavivirus denguei", sequence) &&
    grepl("complete genome|complete cds", sequence) &&
    !grepl("partial", sequence)
}

# Function to process a genome
process_genome <- function(id) {
  sequence <- fetch_genome(id)
  if (is.null(sequence)) return(list(included = 0, excluded = 0, errors = 1))
  if (is_valid_sequence(sequence)) {
    filepath <- file.path(OUTPUT_DIR, sprintf("dengue_genome_%s.fasta", id))
    writeLines(sequence, filepath)
    return(list(included = 1, excluded = 0, errors = 0))
  }
  return(list(included = 0, excluded = 1, errors = 0))
}

# Main function to download genomes
download_dengue_genomes <- function() {
  # Query with organism-level keywords and sequence characteristics
  query <- paste(
    "(Dengue virus type 2[Organism] OR Dengue virus 2[Organism] OR DENV-2[Organism] OR DENV2[Organism] OR Orthoflavivirus denguei[Organism])",
    "AND (complete genome[Title] OR complete cds[Title])",
    sep = " "
  )
  
  retstart <- 0
  stats <- list(included = 0, excluded = 0, errors = 0)
  
  repeat {
    message(sprintf("Searching sequences... Retstart: %d", retstart))
    ids <- search_genomes(retstart, query)
    if (length(ids) == 0) break
    
    for (id in ids) {
      result <- process_genome(id)
      stats$included <- stats$included + result$included
      stats$excluded <- stats$excluded + result$excluded
      stats$errors <- stats$errors + result$errors
    }
    
    retstart <- retstart + BATCH_SIZE
    Sys.sleep(1)
  }
  
  message(sprintf(
    "Processing completed: Included %d, Excluded %d, Errors %d",
    stats$included, stats$excluded, stats$errors
  ))
  return(stats)
}

# Execute the main function
download_dengue_genomes()

###############################################################################
# Merge all downloaded sequences into a single FASTA file for serotype 2
# Required libraries
library(tidyverse)

# Directory where FASTA files are stored
input_dir <- "C:\\Users\\caden\\Documents\\dengue_genomes_DENV2"
output_file <- "C:\\Users\\caden\\Documents\\dengue_genomes_DENV2\\dengue_combined.fasta"

# Function to combine FASTA sequences
combine_fasta_files <- function(input_dir, output_file) {
  # Get all FASTA files from the directory
  fasta_files <- list.files(input_dir, pattern = "\\.fasta$", full.names = TRUE)
  
  if (length(fasta_files) == 0) {
    stop("No FASTA files found in the provided directory.")
  }
  
  # Read and combine the contents of FASTA files
  combined_fasta <- fasta_files %>%
    map_chr(~ readLines(.x) %>% paste(collapse = "\n")) %>%
    paste(collapse = "\n")
  
  # Write the combined file
  writeLines(combined_fasta, output_file)
  message(sprintf("Combined file saved at: %s", output_file))
}

# Execute the function
combine_fasta_files(input_dir, output_file)

#####################################################################################
# Remove gaps (-) and undetermined nucleotides (N) for serotype 2.
# Required libraries
library(Biostrings)

# Input and output file
input_fasta <- "C:\\Users\\caden\\Documents\\dengue_genomes_DENV2\\dengue_combined.FASTA"
output_fasta <- "C:\\Users\\caden\\Documents\\dengue_genomes_DENV2\\dengue_purged.FASTA"

# Function to load and purge FASTA sequences
purge_fasta <- function(input_fasta, output_fasta) {
  # Read the FASTA file
  sequences <- readDNAStringSet(filepath = input_fasta, format = "fasta")
  
  # Initialize counters
  total_sequences <- length(sequences)
  valid_sequences <- 0
  removed_sequences <- 0
  
  # Filter valid sequences
  purged_sequences <- DNAStringSet()
  for (i in seq_along(sequences)) {
    seq <- as.character(sequences[[i]])
    if (!grepl("[-N]", seq)) { # Check for gaps (-) and undetermined nucleotides (N)
      purged_sequences <- append(purged_sequences, sequences[i])
      valid_sequences <- valid_sequences + 1
    } else {
      removed_sequences <- removed_sequences + 1
    }
  }
  
  # Save valid sequences to the output file
  writeXStringSet(purged_sequences, filepath = output_fasta, format = "fasta")
  
  # Print statistics
  message(sprintf("Total sequences reviewed: %d", total_sequences))
  message(sprintf("Valid sequences (without gaps or undetermined nucleotides): %d", valid_sequences))
  message(sprintf("Removed sequences: %d", removed_sequences))
}

# Execute the function
purge_fasta(input_fasta, output_fasta)

#################################################################################
## Download of complete genomic sequences of Dengue virus serotype 3.
library(httr)
library(jsonlite)

# Constants and configuration
BASE_URL <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
API_KEY <- "0142866f7ab066271cca2802824901de9408"
BATCH_SIZE <- 100
RETRY_ATTEMPTS <- 5
RETRY_DELAY <- 3
OUTPUT_DIR <- "C:/Users/caden/Documents/dengue_genomes_DENV3"

# Create output directory
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# Function to perform secure HTTP requests
safe_http_get <- function(url) {
  for (attempt in seq_len(RETRY_ATTEMPTS)) {
    response <- tryCatch({
      GET(url)
    }, error = function(e) NULL)
    if (!is.null(response) && status_code(response) == 200) {
      return(content(response, as = "text", encoding = "UTF-8"))
    }
    Sys.sleep(RETRY_DELAY * attempt)
  }
  return(NULL)
}

# Function to search for genome IDs
search_genomes <- function(retstart, query) {
  url <- sprintf(
    "%s/esearch.fcgi?db=nucleotide&term=%s&retstart=%d&api_key=%s&retmax=%d&retmode=json",
    BASE_URL, URLencode(query), retstart, API_KEY, BATCH_SIZE
  )
  response <- safe_http_get(url)
  if (is.null(response)) return(character(0))
  result <- fromJSON(response, simplifyVector = TRUE)
  return(result$esearchresult$idlist)
}

# Function to download a genome
fetch_genome <- function(id) {
  url <- sprintf(
    "%s/efetch.fcgi?db=nucleotide&id=%s&rettype=fasta&retmode=text&api_key=%s",
    BASE_URL, id, API_KEY
  )
  safe_http_get(url)
}

# Function to validate sequences
is_valid_sequence <- function(sequence) {
  grepl("Dengue virus type 3|Dengue virus 3|DENV-3|DENV3|Orthoflavivirus denguei", sequence) &&
    grepl("complete genome|complete cds", sequence) &&
    !grepl("partial", sequence)
}

# Function to process a genome
process_genome <- function(id) {
  sequence <- fetch_genome(id)
  if (is.null(sequence)) return(list(included = 0, excluded = 0, errors = 1))
  if (is_valid_sequence(sequence)) {
    filepath <- file.path(OUTPUT_DIR, sprintf("dengue_genome_%s.fasta", id))
    writeLines(sequence, filepath)
    return(list(included = 1, excluded = 0, errors = 0))
  }
  return(list(included = 0, excluded = 1, errors = 0))
}

# Main function to download genomes
download_dengue_genomes <- function() {
  # Query with organism-level keywords and sequence characteristics
  query <- paste(
    "(Dengue virus type 3[Organism] OR Dengue virus 3[Organism] OR DENV-3[Organism] OR DENV3[Organism] OR Orthoflavivirus denguei[Organism])",
    "AND (complete genome[Title] OR complete cds[Title])",
    sep = " "
  )
  
  retstart <- 0
  stats <- list(included = 0, excluded = 0, errors = 0)
  
  repeat {
    message(sprintf("Searching sequences... Retstart: %d", retstart))
    ids <- search_genomes(retstart, query)
    if (length(ids) == 0) break
    
    for (id in ids) {
      result <- process_genome(id)
      stats$included <- stats$included + result$included
      stats$excluded <- stats$excluded + result$excluded
      stats$errors <- stats$errors + result$errors
    }
    
    retstart <- retstart + BATCH_SIZE
    Sys.sleep(1)
  }
  
  message(sprintf(
    "Processing completed: Included %d, Excluded %d, Errors %d",
    stats$included, stats$excluded, stats$errors
  ))
  return(stats)
}

# Execute main function
download_dengue_genomes()

###############################################################################
# Merge all downloaded sequences into a single FASTA file for serotype 3.
library(tidyverse)

# Directory where FASTA files are stored
input_dir <- "C:\\Users\\caden\\Documents\\dengue_genomes_DENV3"
output_file <- "C:\\Users\\caden\\Documents\\dengue_genomes_DENV3\\dengue_combined.fasta"

# Function to combine FASTA sequences
combine_fasta_files <- function(input_dir, output_file) {
  # Get all FASTA files in the directory
  fasta_files <- list.files(input_dir, pattern = "\\.fasta$", full.names = TRUE)
  
  if (length(fasta_files) == 0) {
    stop("No FASTA files found in the provided directory.")
  }
  
  # Read and merge FASTA file contents
  combined_fasta <- fasta_files %>%
    map_chr(~ readLines(.x) %>% paste(collapse = "\n")) %>%
    paste(collapse = "\n")
  
  # Write the merged file
  writeLines(combined_fasta, output_file)
  message(sprintf("Merged file saved at: %s", output_file))
}

# Execute the function
combine_fasta_files(input_dir, output_file)

#####################################################################################
# Remove gaps (-) and undetermined nucleotides (N) for serotype 3.
library(Biostrings)

# Input and output files
input_fasta <- "C:\\Users\\caden\\Documents\\dengue_genomes_DENV3\\dengue_combined.FASTA"
output_fasta <- "C:\\Users\\caden\\Documents\\dengue_genomes_DENV3\\dengue_purged.FASTA"

# Function to load and purge FASTA sequences
purge_fasta <- function(input_fasta, output_fasta) {
  # Read the FASTA file
  sequences <- readDNAStringSet(filepath = input_fasta, format = "fasta")
  
  # Initialize counters
  total_sequences <- length(sequences)
  valid_sequences <- 0
  removed_sequences <- 0
  
  # Filter valid sequences
  purged_sequences <- DNAStringSet()
  for (i in seq_along(sequences)) {
    seq <- as.character(sequences[[i]])
    if (!grepl("[-N]", seq)) { # Check for gaps (-) and undetermined nucleotides (N)
      purged_sequences <- append(purged_sequences, sequences[i])
      valid_sequences <- valid_sequences + 1
    } else {
      removed_sequences <- removed_sequences + 1
    }
  }
  
  # Save valid sequences to output file
  writeXStringSet(purged_sequences, filepath = output_fasta, format = "fasta")
  
  # Print statistics
  message(sprintf("Total sequences reviewed: %d", total_sequences))
  message(sprintf("Valid sequences (no gaps or undetermined nucleotides): %d", valid_sequences))
  message(sprintf("Removed sequences: %d", removed_sequences))
}

# Execute the function
purge_fasta(input_fasta, output_fasta)

#################################################################################
## Download complete genomic sequences of Dengue virus serotype 4.
library(httr)
library(jsonlite)

# Constants and configuration
BASE_URL <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
API_KEY <- "0142866f7ab066271cca2802824901de9408"
BATCH_SIZE <- 100
RETRY_ATTEMPTS <- 5
RETRY_DELAY <- 3
OUTPUT_DIR <- "C:/Users/caden/Documents/dengue_genomes_DENV4"

# Create output directory
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# Function for safe HTTP requests
safe_http_get <- function(url) {
  for (attempt in seq_len(RETRY_ATTEMPTS)) {
    response <- tryCatch({
      GET(url)
    }, error = function(e) NULL)
    if (!is.null(response) && status_code(response) == 200) {
      return(content(response, as = "text", encoding = "UTF-8"))
    }
    Sys.sleep(RETRY_DELAY * attempt)
  }
  return(NULL)
}

# Function to search for genome IDs
search_genomes <- function(retstart, query) {
  url <- sprintf(
    "%s/esearch.fcgi?db=nucleotide&term=%s&retstart=%d&api_key=%s&retmax=%d&retmode=json",
    BASE_URL, URLencode(query), retstart, API_KEY, BATCH_SIZE
  )
  response <- safe_http_get(url)
  if (is.null(response)) return(character(0))
  result <- fromJSON(response, simplifyVector = TRUE)
  return(result$esearchresult$idlist)
}

# Function to download a genome
fetch_genome <- function(id) {
  url <- sprintf(
    "%s/efetch.fcgi?db=nucleotide&id=%s&rettype=fasta&retmode=text&api_key=%s",
    BASE_URL, id, API_KEY
  )
  safe_http_get(url)
}

# Function to validate sequences
is_valid_sequence <- function(sequence) {
  grepl("Dengue virus type 4|Dengue virus 4|DENV-4|DENV4|Orthoflavivirus denguei", sequence) &&
    grepl("complete genome|complete cds", sequence) &&
    !grepl("partial", sequence)
}

# Function to process a genome
process_genome <- function(id) {
  sequence <- fetch_genome(id)
  if (is.null(sequence)) return(list(included = 0, excluded = 0, errors = 1))
  if (is_valid_sequence(sequence)) {
    filepath <- file.path(OUTPUT_DIR, sprintf("dengue_genome_%s.fasta", id))
    writeLines(sequence, filepath)
    return(list(included = 1, excluded = 0, errors = 0))
  }
  return(list(included = 0, excluded = 1, errors = 0))
}

# Main function to download genomes
download_dengue_genomes <- function() {
  # Query with keywords at the organism level and sequence features
  query <- paste(
    "(Dengue virus type 4[Organism] OR Dengue virus 4[Organism] OR DENV-4[Organism] OR DENV4[Organism] OR Orthoflavivirus denguei[Organism])",
    "AND (complete genome[Title] OR complete cds[Title])",
    sep = " "
  )
  
  retstart <- 0
  stats <- list(included = 0, excluded = 0, errors = 0)
  
  repeat {
    message(sprintf("Searching sequences... Retstart: %d", retstart))
    ids <- search_genomes(retstart, query)
    if (length(ids) == 0) break
    
    for (id in ids) {
      result <- process_genome(id)
      stats$included <- stats$included + result$included
      stats$excluded <- stats$excluded + result$excluded
      stats$errors <- stats$errors + result$errors
    }
    
    retstart <- retstart + BATCH_SIZE
    Sys.sleep(1)
  }
  
  message(sprintf(
    "Processing completed: Included %d, Excluded %d, Errors %d",
    stats$included, stats$excluded, stats$errors
  ))
  return(stats)
}

# Execute the main function
download_dengue_genomes()

###############################################################################
# Combine all downloaded sequences into a single FASTA file for serotype 4.
# Required libraries
library(tidyverse)

# Directory where FASTA files are stored
input_dir <- "C:\\Users\\caden\\Documents\\dengue_genomes_DENV4"
output_file <- "C:\\Users\\caden\\Documents\\dengue_genomes_DENV4\\dengue_combined.fasta"

# Function to combine FASTA files
combine_fasta_files <- function(input_dir, output_file) {
  # Get all FASTA files in the directory
  fasta_files <- list.files(input_dir, pattern = "\\.fasta$", full.names = TRUE)
  
  if (length(fasta_files) == 0) {
    stop("No FASTA files found in the provided directory.")
  }
  
  # Read and combine the contents of the FASTA files
  combined_fasta <- fasta_files %>%
    map_chr(~ readLines(.x) %>% paste(collapse = "\n")) %>%
    paste(collapse = "\n")
  
  # Write the combined file
  writeLines(combined_fasta, output_file)
  message(sprintf("Combined file saved at: %s", output_file))
}

# Execute the function
combine_fasta_files(input_dir, output_file)

#####################################################################################
# Remove gaps (-) and indeterminate nucleotides (N) for serotype 4.
# Required libraries
library(Biostrings)

# Input and output files
input_fasta <- "C:\\Users\\caden\\Documents\\dengue_genomes_DENV4\\dengue_combined.FASTA"
output_fasta <- "C:\\Users\\caden\\Documents\\dengue_genomes_DENV4\\dengue_purged.FASTA"

# Function to load and purge FASTA sequences
purge_fasta <- function(input_fasta, output_fasta) {
  # Read the FASTA file
  sequences <- readDNAStringSet(filepath = input_fasta, format = "fasta")
  
  # Initialize counters
  total_sequences <- length(sequences)
  valid_sequences <- 0
  removed_sequences <- 0
  
  # Filter valid sequences
  purged_sequences <- DNAStringSet()
  for (i in seq_along(sequences)) {
    seq <- as.character(sequences[[i]])
    if (!grepl("[-N]", seq)) { # Check for gaps (-) and indeterminate nucleotides (N)
      purged_sequences <- append(purged_sequences, sequences[i])
      valid_sequences <- valid_sequences + 1
    } else {
      removed_sequences <- removed_sequences + 1
    }
  }
  
  # Save valid sequences to the output file
  writeXStringSet(purged_sequences, filepath = output_fasta, format = "fasta")
  
  # Print statistics
  message(sprintf("Total sequences reviewed: %d", total_sequences))
  message(sprintf("Valid sequences (no gaps or indeterminate nucleotides): %d", valid_sequences))
  message(sprintf("Removed sequences: %d", removed_sequences))
}

# Execute the function
purge_fasta(input_fasta, output_fasta)
