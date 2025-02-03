###############################################################################
#Descarga de secuencias genomicas completas del virus del dengue serotipo 1.
# Librerías necesarias
library(httr)
library(jsonlite)

# Constantes y configuración
BASE_URL <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
API_KEY <- "0142866f7ab066271cca2802824901de9408"
BATCH_SIZE <- 100
RETRY_ATTEMPTS <- 5
RETRY_DELAY <- 3
OUTPUT_DIR <- "C:/Users/caden/Documents/dengue_genomes"

# Crear directorio de salida
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# Función para realizar solicitudes HTTP seguras
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

# Función para buscar IDs de genomas
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

# Función para descargar un genoma
fetch_genome <- function(id) {
  url <- sprintf(
    "%s/efetch.fcgi?db=nucleotide&id=%s&rettype=fasta&retmode=text&api_key=%s",
    BASE_URL, id, API_KEY
  )
  safe_http_get(url)
}

# Función para validar secuencias
is_valid_sequence <- function(sequence) {
  grepl("Dengue virus type 1|Dengue virus 1|DENV-1|Orthoflavivirus denguei", sequence) &&
    grepl("complete genome|complete cds", sequence) &&
    !grepl("partial", sequence)
}

# Función para procesar un genoma
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

# Función principal para descargar genomas
download_dengue_genomes <- function() {
  # Query con las palabras clave a nivel de organismo y características de la secuencia
  query <- paste(
    "(Dengue virus type 1[Organism] OR Dengue virus 1[Organism] OR DENV-1[Organism] OR Orthoflavivirus denguei[Organism])",
    "AND (complete genome[Title] OR complete cds[Title])",
    sep = " "
  )
  
  retstart <- 0
  stats <- list(included = 0, excluded = 0, errors = 0)
  
  repeat {
    message(sprintf("Buscando secuencias... Retstart: %d", retstart))
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
    "Procesamiento completado: Incluidos %d, Excluidos %d, Errores %d",
    stats$included, stats$excluded, stats$errors
  ))
  return(stats)
}

# Ejecutar la función principal
download_dengue_genomes()

###############################################################################
#Unir todas las secuencias descargadas en un unico archivo FASTA serotipo 1
# Librerías necesarias
library(tidyverse)

# Directorio donde están almacenados los archivos FASTA
input_dir <- "C:\\Users\\caden\\Documents\\dengue_genomes_DENV1"
output_file <- "C:\\Users\\caden\\Documents\\dengue_genomes_DENV1\\dengue_combined.fasta"

# Función para combinar secuencias FASTA
combine_fasta_files <- function(input_dir, output_file) {
  # Obtener todos los archivos FASTA del directorio
  fasta_files <- list.files(input_dir, pattern = "\\.fasta$", full.names = TRUE)
  
  if (length(fasta_files) == 0) {
    stop("No se encontraron archivos FASTA en el directorio proporcionado.")
  }
  
  # Leer y combinar los contenidos de los archivos FASTA
  combined_fasta <- fasta_files %>%
    map_chr(~ readLines(.x) %>% paste(collapse = "\n")) %>%
    paste(collapse = "\n")
  
  # Escribir el archivo combinado
  writeLines(combined_fasta, output_file)
  message(sprintf("Archivo combinado guardado en: %s", output_file))
}

# Ejecutar la función
combine_fasta_files(input_dir, output_file)

################################################################################
#Eliminar gaps (-) y nucleótidos indeterminados (N) para serotipo 1.

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("Biostrings")

# Librerías necesarias
library(Biostrings)

# Archivo de entrada y salida
input_fasta <- "C:\\Users\\caden\\Documents\\dengue_genomes_DENV1\\dengue_combined.FASTA"
output_fasta <- "C:\\Users\\caden\\Documents\\dengue_genomes_DENV1\\dengue_purged.FASTA"

# Función para cargar y purgar secuencias FASTA
purge_fasta <- function(input_fasta, output_fasta) {
  # Leer el archivo FASTA
  sequences <- readDNAStringSet(filepath = input_fasta, format = "fasta")
  
  # Inicializar contadores
  total_sequences <- length(sequences)
  valid_sequences <- 0
  removed_sequences <- 0
  
  # Filtrar secuencias válidas
  purged_sequences <- DNAStringSet()
  for (i in seq_along(sequences)) {
    seq <- as.character(sequences[[i]])
    if (!grepl("[-N]", seq)) { # Revisión de gaps (-) y nucleótidos indeterminados (N)
      purged_sequences <- append(purged_sequences, sequences[i])
      valid_sequences <- valid_sequences + 1
    } else {
      removed_sequences <- removed_sequences + 1
    }
  }
  
  # Guardar las secuencias válidas en el archivo de salida
  writeXStringSet(purged_sequences, filepath = output_fasta, format = "fasta")
  
  # Imprimir estadísticas
  message(sprintf("Total de secuencias revisadas: %d", total_sequences))
  message(sprintf("Secuencias válidas (sin gaps ni nucleótidos indeterminados): %d", valid_sequences))
  message(sprintf("Secuencias eliminadas: %d", removed_sequences))
}

# Ejecutar la función
purge_fasta(input_fasta, output_fasta)

###############################################################################
#Descarga de secuencias genomicas completas del virus del dengue serotipo 2.
# Librerías necesarias
library(httr)
library(jsonlite)

# Constantes y configuración
BASE_URL <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
API_KEY <- "0142866f7ab066271cca2802824901de9408"
BATCH_SIZE <- 100
RETRY_ATTEMPTS <- 5
RETRY_DELAY <- 3
OUTPUT_DIR <- "C:/Users/caden/Documents/dengue_genomes_DENV2"

# Crear directorio de salida
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# Función para realizar solicitudes HTTP seguras
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

# Función para buscar IDs de genomas
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

# Función para descargar un genoma
fetch_genome <- function(id) {
  url <- sprintf(
    "%s/efetch.fcgi?db=nucleotide&id=%s&rettype=fasta&retmode=text&api_key=%s",
    BASE_URL, id, API_KEY
  )
  safe_http_get(url)
}

# Función para validar secuencias
is_valid_sequence <- function(sequence) {
  grepl("Dengue virus type 2|Dengue virus 2|DENV-2|DENV2|Orthoflavivirus denguei", sequence) &&
    grepl("complete genome|complete cds", sequence) &&
    !grepl("partial", sequence)
}

# Función para procesar un genoma
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

# Función principal para descargar genomas
download_dengue_genomes <- function() {
  # Query con las palabras clave a nivel de organismo y características de la secuencia
  query <- paste(
    "(Dengue virus type 2[Organism] OR Dengue virus 2[Organism] OR DENV-2[Organism] OR DENV2[Organism] OR Orthoflavivirus denguei[Organism])",
    "AND (complete genome[Title] OR complete cds[Title])",
    sep = " "
  )
  
  retstart <- 0
  stats <- list(included = 0, excluded = 0, errors = 0)
  
  repeat {
    message(sprintf("Buscando secuencias... Retstart: %d", retstart))
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
    "Procesamiento completado: Incluidos %d, Excluidos %d, Errores %d",
    stats$included, stats$excluded, stats$errors
  ))
  return(stats)
}

# Ejecutar la función principal
download_dengue_genomes()

###############################################################################
#Unir todas las secuencias descargadas en un unico archivo FASTA para serotipo 2
# Librerías necesarias
library(tidyverse)

# Directorio donde están almacenados los archivos FASTA
input_dir <- "C:\\Users\\caden\\Documents\\dengue_genomes_DENV2"
output_file <- "C:\\Users\\caden\\Documents\\dengue_genomes_DENV2\\dengue_combined.fasta"

# Función para combinar secuencias FASTA
combine_fasta_files <- function(input_dir, output_file) {
  # Obtener todos los archivos FASTA del directorio
  fasta_files <- list.files(input_dir, pattern = "\\.fasta$", full.names = TRUE)
  
  if (length(fasta_files) == 0) {
    stop("No se encontraron archivos FASTA en el directorio proporcionado.")
  }
  
  # Leer y combinar los contenidos de los archivos FASTA
  combined_fasta <- fasta_files %>%
    map_chr(~ readLines(.x) %>% paste(collapse = "\n")) %>%
    paste(collapse = "\n")
  
  # Escribir el archivo combinado
  writeLines(combined_fasta, output_file)
  message(sprintf("Archivo combinado guardado en: %s", output_file))
}

# Ejecutar la función
combine_fasta_files(input_dir, output_file)

#####################################################################################
#Eliminar gaps (-) y nucleótidos indeterminados (N) para serotipo 2.
# Librerías necesarias
library(Biostrings)

# Archivo de entrada y salida
input_fasta <- "C:\\Users\\caden\\Documents\\dengue_genomes_DENV2\\dengue_combined.FASTA"
output_fasta <- "C:\\Users\\caden\\Documents\\dengue_genomes_DENV2\\dengue_purged.FASTA"

# Función para cargar y purgar secuencias FASTA
purge_fasta <- function(input_fasta, output_fasta) {
  # Leer el archivo FASTA
  sequences <- readDNAStringSet(filepath = input_fasta, format = "fasta")
  
  # Inicializar contadores
  total_sequences <- length(sequences)
  valid_sequences <- 0
  removed_sequences <- 0
  
  # Filtrar secuencias válidas
  purged_sequences <- DNAStringSet()
  for (i in seq_along(sequences)) {
    seq <- as.character(sequences[[i]])
    if (!grepl("[-N]", seq)) { # Revisión de gaps (-) y nucleótidos indeterminados (N)
      purged_sequences <- append(purged_sequences, sequences[i])
      valid_sequences <- valid_sequences + 1
    } else {
      removed_sequences <- removed_sequences + 1
    }
  }
  
  # Guardar las secuencias válidas en el archivo de salida
  writeXStringSet(purged_sequences, filepath = output_fasta, format = "fasta")
  
  # Imprimir estadísticas
  message(sprintf("Total de secuencias revisadas: %d", total_sequences))
  message(sprintf("Secuencias válidas (sin gaps ni nucleótidos indeterminados): %d", valid_sequences))
  message(sprintf("Secuencias eliminadas: %d", removed_sequences))
}

# Ejecutar la función
purge_fasta(input_fasta, output_fasta)

#################################################################################
##Descarga de secuencias genomicas completas del virus del dengue serotipo 3.
library(httr)
library(jsonlite)

# Constantes y configuración
BASE_URL <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
API_KEY <- "0142866f7ab066271cca2802824901de9408"
BATCH_SIZE <- 100
RETRY_ATTEMPTS <- 5
RETRY_DELAY <- 3
OUTPUT_DIR <- "C:/Users/caden/Documents/dengue_genomes_DENV3"

# Crear directorio de salida
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# Función para realizar solicitudes HTTP seguras
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

# Función para buscar IDs de genomas
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

# Función para descargar un genoma
fetch_genome <- function(id) {
  url <- sprintf(
    "%s/efetch.fcgi?db=nucleotide&id=%s&rettype=fasta&retmode=text&api_key=%s",
    BASE_URL, id, API_KEY
  )
  safe_http_get(url)
}

# Función para validar secuencias
is_valid_sequence <- function(sequence) {
  grepl("Dengue virus type 3|Dengue virus 3|DENV-3|DENV3|Orthoflavivirus denguei", sequence) &&
    grepl("complete genome|complete cds", sequence) &&
    !grepl("partial", sequence)
}

# Función para procesar un genoma
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

# Función principal para descargar genomas
download_dengue_genomes <- function() {
  # Query con las palabras clave a nivel de organismo y características de la secuencia
  query <- paste(
    "(Dengue virus type 3[Organism] OR Dengue virus 3[Organism] OR DENV-3[Organism] OR DENV3[Organism] OR Orthoflavivirus denguei[Organism])",
    "AND (complete genome[Title] OR complete cds[Title])",
    sep = " "
  )
  
  retstart <- 0
  stats <- list(included = 0, excluded = 0, errors = 0)
  
  repeat {
    message(sprintf("Buscando secuencias... Retstart: %d", retstart))
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
    "Procesamiento completado: Incluidos %d, Excluidos %d, Errores %d",
    stats$included, stats$excluded, stats$errors
  ))
  return(stats)
}

# Ejecutar la función principal
download_dengue_genomes()

###############################################################################
#Unir todas las secuencias descargadas en un unico archivo FASTA para serotipo 3.
# Librerías necesarias
library(tidyverse)

# Directorio donde están almacenados los archivos FASTA
input_dir <- "C:\\Users\\caden\\Documents\\dengue_genomes_DENV3"
output_file <- "C:\\Users\\caden\\Documents\\dengue_genomes_DENV3\\dengue_combined.fasta"

# Función para combinar secuencias FASTA
combine_fasta_files <- function(input_dir, output_file) {
  # Obtener todos los archivos FASTA del directorio
  fasta_files <- list.files(input_dir, pattern = "\\.fasta$", full.names = TRUE)
  
  if (length(fasta_files) == 0) {
    stop("No se encontraron archivos FASTA en el directorio proporcionado.")
  }
  
  # Leer y combinar los contenidos de los archivos FASTA
  combined_fasta <- fasta_files %>%
    map_chr(~ readLines(.x) %>% paste(collapse = "\n")) %>%
    paste(collapse = "\n")
  
  # Escribir el archivo combinado
  writeLines(combined_fasta, output_file)
  message(sprintf("Archivo combinado guardado en: %s", output_file))
}

# Ejecutar la función
combine_fasta_files(input_dir, output_file)

#####################################################################################
#Eliminar gaps (-) y nucleótidos indeterminados (N) para serotipo 3.
# Librerías necesarias
library(Biostrings)

# Archivo de entrada y salida
input_fasta <- "C:\\Users\\caden\\Documents\\dengue_genomes_DENV3\\dengue_combined.FASTA"
output_fasta <- "C:\\Users\\caden\\Documents\\dengue_genomes_DENV3\\dengue_purged.FASTA"

# Función para cargar y purgar secuencias FASTA
purge_fasta <- function(input_fasta, output_fasta) {
  # Leer el archivo FASTA
  sequences <- readDNAStringSet(filepath = input_fasta, format = "fasta")
  
  # Inicializar contadores
  total_sequences <- length(sequences)
  valid_sequences <- 0
  removed_sequences <- 0
  
  # Filtrar secuencias válidas
  purged_sequences <- DNAStringSet()
  for (i in seq_along(sequences)) {
    seq <- as.character(sequences[[i]])
    if (!grepl("[-N]", seq)) { # Revisión de gaps (-) y nucleótidos indeterminados (N)
      purged_sequences <- append(purged_sequences, sequences[i])
      valid_sequences <- valid_sequences + 1
    } else {
      removed_sequences <- removed_sequences + 1
    }
  }
  
  # Guardar las secuencias válidas en el archivo de salida
  writeXStringSet(purged_sequences, filepath = output_fasta, format = "fasta")
  
  # Imprimir estadísticas
  message(sprintf("Total de secuencias revisadas: %d", total_sequences))
  message(sprintf("Secuencias válidas (sin gaps ni nucleótidos indeterminados): %d", valid_sequences))
  message(sprintf("Secuencias eliminadas: %d", removed_sequences))
}

# Ejecutar la función
purge_fasta(input_fasta, output_fasta)

#################################################################################
##Descarga de secuencias genomicas completas del virus del dengue serotipo 4.
library(httr)
library(jsonlite)

# Constantes y configuración
BASE_URL <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
API_KEY <- "0142866f7ab066271cca2802824901de9408"
BATCH_SIZE <- 100
RETRY_ATTEMPTS <- 5
RETRY_DELAY <- 3
OUTPUT_DIR <- "C:/Users/caden/Documents/dengue_genomes_DENV4"

# Crear directorio de salida
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# Función para realizar solicitudes HTTP seguras
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

# Función para buscar IDs de genomas
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

# Función para descargar un genoma
fetch_genome <- function(id) {
  url <- sprintf(
    "%s/efetch.fcgi?db=nucleotide&id=%s&rettype=fasta&retmode=text&api_key=%s",
    BASE_URL, id, API_KEY
  )
  safe_http_get(url)
}

# Función para validar secuencias
is_valid_sequence <- function(sequence) {
  grepl("Dengue virus type 4|Dengue virus 4|DENV-4|DENV4|Orthoflavivirus denguei", sequence) &&
    grepl("complete genome|complete cds", sequence) &&
    !grepl("partial", sequence)
}

# Función para procesar un genoma
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

# Función principal para descargar genomas
download_dengue_genomes <- function() {
  # Query con las palabras clave a nivel de organismo y características de la secuencia
  query <- paste(
    "(Dengue virus type 4[Organism] OR Dengue virus 4[Organism] OR DENV-4[Organism] OR DENV4[Organism] OR Orthoflavivirus denguei[Organism])",
    "AND (complete genome[Title] OR complete cds[Title])",
    sep = " "
  )
  
  retstart <- 0
  stats <- list(included = 0, excluded = 0, errors = 0)
  
  repeat {
    message(sprintf("Buscando secuencias... Retstart: %d", retstart))
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
    "Procesamiento completado: Incluidos %d, Excluidos %d, Errores %d",
    stats$included, stats$excluded, stats$errors
  ))
  return(stats)
}

# Ejecutar la función principal
download_dengue_genomes()

###############################################################################
#Unir todas las secuencias descargadas en un unico archivo FASTA para serotipo 4.
# Librerías necesarias
library(tidyverse)

# Directorio donde están almacenados los archivos FASTA
input_dir <- "C:\\Users\\caden\\Documents\\dengue_genomes_DENV4"
output_file <- "C:\\Users\\caden\\Documents\\dengue_genomes_DENV4\\dengue_combined.fasta"

# Función para combinar secuencias FASTA
combine_fasta_files <- function(input_dir, output_file) {
  # Obtener todos los archivos FASTA del directorio
  fasta_files <- list.files(input_dir, pattern = "\\.fasta$", full.names = TRUE)
  
  if (length(fasta_files) == 0) {
    stop("No se encontraron archivos FASTA en el directorio proporcionado.")
  }
  
  # Leer y combinar los contenidos de los archivos FASTA
  combined_fasta <- fasta_files %>%
    map_chr(~ readLines(.x) %>% paste(collapse = "\n")) %>%
    paste(collapse = "\n")
  
  # Escribir el archivo combinado
  writeLines(combined_fasta, output_file)
  message(sprintf("Archivo combinado guardado en: %s", output_file))
}

# Ejecutar la función
combine_fasta_files(input_dir, output_file)

#####################################################################################
#Eliminar gaps (-) y nucleótidos indeterminados (N) para serotipo 4.
# Librerías necesarias
library(Biostrings)

# Archivo de entrada y salida
input_fasta <- "C:\\Users\\caden\\Documents\\dengue_genomes_DENV4\\dengue_combined.FASTA"
output_fasta <- "C:\\Users\\caden\\Documents\\dengue_genomes_DENV4\\dengue_purged.FASTA"

# Función para cargar y purgar secuencias FASTA
purge_fasta <- function(input_fasta, output_fasta) {
  # Leer el archivo FASTA
  sequences <- readDNAStringSet(filepath = input_fasta, format = "fasta")
  
  # Inicializar contadores
  total_sequences <- length(sequences)
  valid_sequences <- 0
  removed_sequences <- 0
  
  # Filtrar secuencias válidas
  purged_sequences <- DNAStringSet()
  for (i in seq_along(sequences)) {
    seq <- as.character(sequences[[i]])
    if (!grepl("[-N]", seq)) { # Revisión de gaps (-) y nucleótidos indeterminados (N)
      purged_sequences <- append(purged_sequences, sequences[i])
      valid_sequences <- valid_sequences + 1
    } else {
      removed_sequences <- removed_sequences + 1
    }
  }
  
  # Guardar las secuencias válidas en el archivo de salida
  writeXStringSet(purged_sequences, filepath = output_fasta, format = "fasta")
  
  # Imprimir estadísticas
  message(sprintf("Total de secuencias revisadas: %d", total_sequences))
  message(sprintf("Secuencias válidas (sin gaps ni nucleótidos indeterminados): %d", valid_sequences))
  message(sprintf("Secuencias eliminadas: %d", removed_sequences))
}

# Ejecutar la función
purge_fasta(input_fasta, output_fasta)