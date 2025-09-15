# This script contains the core logic for downloading reference files from IMGT
# and the IPD-IMGT/HLA database.

# --- Constants and Mappings for IMGT TCR/BCR Queries ---

.species_map <- list(
  human = "Homo+sapiens",
  mouse = "Mus",
  rat = "Rattus+norvegicus",
  rabbit = "Oryctolagus+cuniculus",
  pig = "Sus+scrofa",
  dog = "Canis+lupus+familiaris", 
  "rhesus_monkey" = "Macaca+mulatta",
  "cyno_monkey" = "Macaca+fascicularis"
)

.species_replace_map <- list(
  human = "Homo sapiens",
  mouse = "Mus musculus",
  rat = "Rattus norvegicus",
  rabbit = "Oryctolagus cuniculus",
  pig = "Sus+scrofa",
  dog = "Canis+lupus+familiaris",
  "rhesus_monkey" = "Macaca mulatta",
  "cyno_monkey" = "Macaca fascicularis"
)

# --- HLA Download Logic ---

.hla_base_url <- "ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/"

#' @title Fetch a Single FASTA File from an IMGT Query
#' @description Internal function to perform the HTTP request and parsing for a
#' single IMGT GENE-DB query.
#' @param url The full IMGT query URL.
#' @param dest_file The path where the downloaded FASTA file should be saved.
#' @param species_key The key for the species (e.g., "human") for cleaning headers.
#' @importFrom httr RETRY status_code content
#' @importFrom rvest read_html html_nodes html_text
#' @noRd
.fetch_imgt_query <- function(url, dest_file, species_key) {
  tryCatch({
    resp <- httr::RETRY("GET", url, times = 3, pause_base = 5, pause_cap = 60)
    if (httr::status_code(resp) >= 400) {
      warning("Failed to download ", url, " with status ", httr::status_code(resp))
      return(invisible(NULL))
    }

    fasta_text_nodes <- httr::content(resp, "text", encoding = "UTF-8") |>
      rvest::read_html() |>
      rvest::html_nodes("pre")

    if (length(fasta_text_nodes) < 2) {
      warning("Could not find FASTA content in ", url)
      return(invisible(NULL))
    }
    fasta_text <- fasta_text_nodes[[2]] |> rvest::html_text()

    if (nchar(fasta_text) > 0) {
      clean_species_name <- gsub(" ", "_", .species_replace_map[[species_key]])
      fasta_text <- gsub(.species_replace_map[[species_key]], clean_species_name, fasta_text)
      writeLines(fasta_text, dest_file)
    } else {
      warning("Empty FASTA content from ", url)
    }
  }, error = function(e) {
    warning("Error downloading or processing ", url, ": ", e$message)
  })
  return(invisible(NULL))
}


#' @title Fetch a Single HLA FASTA File from GitHub
#' @description Internal function to download a pre-compiled HLA FASTA file.
#' @param file_type The type of HLA sequence to download. One of "nuc", "prot".
#' @param dest_file The path where the downloaded FASTA file should be saved.
#' @importFrom httr write_disk
#' @noRd
.fetch_hla_file <- function(file_type = "nuc", dest_file) {
  if (!file_type %in% c("nuc", "prot")) {
    stop("file_type must be one of 'nuc' or 'prot'.")
  }

  file_name <- paste0("hla_", file_type, ".fasta")
  url <- paste0(.hla_base_url, file_name)

  tryCatch({
    resp <- httr::RETRY("GET", url, times = 3, pause_base = 5, pause_cap = 60,
                        httr::write_disk(dest_file, overwrite = TRUE))
    if (httr::status_code(resp) >= 400) {
      warning("Failed to download ", url, " with status ", httr::status_code(resp))
    }
  }, error = function(e) {
    warning("Error downloading ", url, ": ", e$message)
  })
  return(invisible(NULL))
}
