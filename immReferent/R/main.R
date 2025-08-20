# This file contains the main user-facing functions for the immReferent package.

#' @import Biostrings
#' @importFrom methods new
#' @importFrom utils read.csv

# A helper to show the license message once per session
.show_license_message <- function(suppress = FALSE) {
  if (!suppress && !isTRUE(getOption("immReferent.license.shown"))) {
    message("Data from IMGT (https://www.imgt.org) is for non-commercial use only, under a CC-BY-NC-ND 4.0 license. Attribution is required.")
    message("HLA data from IPD-IMGT/HLA (https://www.ebi.ac.uk/ipd/imgt/hla/) is provided under the CC-BY-ND 4.0 license.")
    options(immReferent.license.shown = TRUE)
  }
}

#' @title Download and Load Immune Receptor and HLA Sequences
#' @description This is the main function to download and load reference sequences from IMGT
#' and the IPD-IMGT/HLA database. It handles caching of downloaded files.
#'
#' @param species The species for which to download data. Required for TCR/BCR queries.
#'        Currently supported: "human", "mouse", "rat", "rabbit", "rhesus_monkey". Defaults to "human" for HLA.
#' @param gene The gene or locus to download. For TCR/BCR, this can be a specific
#'        chain (e.g., "IGHV", "TRBJ") or a group (e.g., "IGH", "TCR"). For HLA, use "HLA".
#' @param type The type of sequence to retrieve. Either "NUC" for nucleotide or
#'        "PROT" for protein sequences. This primarily distinguishes between VDJ nucleotide
#'        and V-region amino acid sequences for TCR/BCR genes.
#' @param refresh Logical. If `TRUE`, forces a re-download of the data even if it
#'        exists in the cache.
#' @param suppressMessages Logical. If `TRUE`, suppresses the license and other
#'        informational messages.
#'
#' @return A `DNAStringSet` or `AAStringSet` object containing the requested sequences.
#' @export
#' @examples
#' \dontrun{
#'   # Download human IGHV nucleotide sequences
#'   ighv_nuc <- getIMGT(species = "human", gene = "IGHV", type = "NUC")
#'
#'   # Download all HLA protein sequences
#'   hla_prot <- getIMGT(gene = "HLA", type = "PROT")
#'
#'   # Download all mouse TRB genes
#'   trb_mouse <- getIMGT(species = "mouse", gene = "TRB", type = "NUC")
#' }
getIMGT <- function(species = "human", gene, type = c("NUC", "PROT"),
                    refresh = FALSE, suppressMessages = FALSE) {

  .show_license_message(suppress = suppressMessages)
  type <- match.arg(type)
  gene_orig <- gene
  gene <- toupper(gene)

  cache_dir <- .get_cache_dir()
  .ensure_cache_dir()

  # --- Handle HLA data ---
  if (startsWith(gene, "HLA")) {
    if (!suppressMessages) message("-> Handling HLA request.")
    hla_file_type <- ifelse(type == "NUC", "nuc", "prot")
    hla_dir <- file.path(cache_dir, "human", "hla")
    if (!dir.exists(hla_dir)) dir.create(hla_dir, recursive = TRUE)
    dest_file <- file.path(hla_dir, paste0("hla_", hla_file_type, ".fasta"))
    if (refresh || !file.exists(dest_file)) {
      if (!is_imgt_available()) stop("Internet connection unavailable.")
      if (!suppressMessages) message("-> Downloading HLA ", type, " sequences...")
      .fetch_hla_file(file_type = hla_file_type, dest_file = dest_file)
      .update_cache_metadata(species = "human", gene = gene_orig)
    } else {
      if (!suppressMessages) message("-> Loading HLA ", type, " sequences from cache.")
    }
    if (!file.exists(dest_file)) stop("Failed to find or download HLA data.")
    if (type == "NUC") return(Biostrings::readDNAStringSet(dest_file))
    else return(Biostrings::readAAStringSet(dest_file))
  }

  # --- Handle TCR/BCR data ---
  else {
    if (!species %in% names(.species_map)) {
      stop("Species '", species, "' is not supported for TCR/BCR queries.")
    }

    target_genes <- if (gene %in% names(.gene_groups)) .gene_groups[[gene]] else gene
    query_type <- ifelse(type == "NUC", "NUC", "PROT")
    files_to_load <- c()

    for (g in target_genes) {
      params <- .imgt_db_map[.imgt_db_map$gene == g & .imgt_db_map$type == query_type, ]
      if (nrow(params) == 0) {
        params <- .imgt_db_map[.imgt_db_map$gene == g & .imgt_db_map$type == "NUC-C", ]
        if (nrow(params) == 0) {
            if (!suppressMessages) warning("No rule to download '", g, "' with type '", type, "'. Skipping.")
            next
        }
      }

      target_dir <- file.path(cache_dir, species, params$cache_subdir)
      if (!dir.exists(target_dir)) dir.create(target_dir, recursive = TRUE)
      file_name <- paste0(params$filename_prefix, species, "_", g, ".fasta")
      dest_file <- file.path(target_dir, file_name)
      files_to_load <- c(files_to_load, dest_file)

      if (refresh || !file.exists(dest_file)) {
        if (!is_imgt_available()) stop("Internet connection unavailable.")
        if (!suppressMessages) message("-> Downloading ", g, " (", species, ", ", type, ")")

        species_query <- .species_map[[species]]
        query <- paste0(params$query_prefix, "+", g)
        url <- paste0("https://www.imgt.org/genedb/GENElect?query=", query, "&species=", species_query)
        if (params$query_label != "") url <- paste0(url, "&IMGTlabel=", params$query_label)

        .fetch_imgt_query(url = url, dest_file = dest_file, species_key = species)
        .update_cache_metadata(species = species, gene = g)
      } else {
        if (!suppressMessages) message("-> Found ", g, " in cache.")
      }
    }

    if (length(files_to_load) == 0) {
      warning("No data could be found or downloaded for the request.")
      return(NULL)
    }

    if (!suppressMessages) message("-> Loading ", length(files_to_load), " file(s) into Biostrings object.")

    all_sets <- lapply(files_to_load, function(f) {
      if (!file.exists(f) || file.info(f)$size == 0) return(NULL)
      if (type == "NUC") Biostrings::readDNAStringSet(f)
      else Biostrings::readAAStringSet(f)
    })

    combined_set <- do.call(c, all_sets[!sapply(all_sets, is.null)])
    return(combined_set)
  }
}


#' @title Load Cached IMGT/HLA Sequences
#' @description Loads sequences from the local cache without attempting to download.
#' This function relies on `getIMGT(refresh = FALSE)`. If the data is not found
#' in the cache, it will be downloaded unless an internet connection is unavailable.
#'
#' @inheritParams getIMGT
#' @return A `DNAStringSet` or `AAStringSet` object.
#' @export
loadIMGT <- function(species = "human", gene, type = c("NUC", "PROT"), suppressMessages = FALSE) {
  if (!suppressMessages) {
    message("-> `loadIMGT` is a wrapper for `getIMGT(refresh = FALSE)`. It will load from cache if available.")
  }
  getIMGT(species = species, gene = gene, type = type, refresh = FALSE, suppressMessages = suppressMessages)
}


#' @title Force Re-download of IMGT/HLA Sequences
#' @description A convenience wrapper for `getIMGT(..., refresh = TRUE)` to ensure that
#' the local cache is updated with the latest versions of the requested sequences.
#'
#' @inheritParams getIMGT
#' @return A `DNAStringSet` or `AAStringSet` object.
#' @export
refreshIMGT <- function(species = "human", gene, type = c("NUC", "PROT"), suppressMessages = FALSE) {
  getIMGT(species = species, gene = gene, type = type, refresh = TRUE, suppressMessages = suppressMessages)
}


#' @title List Datasets in the Local Cache
#' @description Scans the cache directory and returns a list of available datasets.
#'
#' @return A character vector of file paths for the cached datasets.
#' @export
listIMGT <- function() {
  cache_dir <- .get_cache_dir()
  if (!dir.exists(cache_dir)) {
    message("Cache directory does not exist. No datasets found.")
    return(character(0))
  }
  list.files(cache_dir, recursive = TRUE, full.names = TRUE)
}
