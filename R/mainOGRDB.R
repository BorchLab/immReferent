# This file contains the main OGRDB user-facing functions for the immReferent package.


# Accept flexible inputs and normalize to OGRDB expectations.
.ogrdb_species_map <- list(
  "human" = list(api="Human", taxid="9606", binomial="Homo sapiens"),
  "homo sapiens" = list(api="Human", taxid="9606", binomial="Homo sapiens"),
  "mouse" = list(api="Mouse", taxid="10090", binomial="Mus musculus"),
  "mus musculus" = list(api="Mouse", taxid="10090", binomial="Mus musculus")
  # TODO expand as more species added
)

# Convenience: map locus -> default set name (Human IG reference sets)
.ogrdb_default_set <- list(
  IGH = "IGH_VDJ",
  IGK = "IGKappa_VJ",
  IGL = "IGLambda_VJ"
  #TODO When TCR sets arrive, add e.g. TRA="TRA_VJ", TRB="TRB_VDJ"
)

# ---- URL builder for OGRDB REST ---------------------------------------------
.ogrdb_build_url <- function(species_api, 
                             set_name, 
                             version = c("published","latest"),
                             format = c("gapped","ungapped","airr"),
                             species_subgroup = NULL,
                             human_reference = TRUE,
                             use_api = TRUE) {
  version <- match.arg(version)
  format  <- match.arg(format)
  
  # For Human Reference Set, append _ex to the FORMAT (not the set name)
  fmt_final <- if (identical(species_api, "Human") && human_reference) {
    paste0(format, "_ex")
  } else format
  
  parts <- c(if (use_api) "https://ogrdb.airr-community.org/api/germline/set"
             else          "https://ogrdb.airr-community.org/download_germline_set",
             utils::URLencode(species_api, reserved = TRUE))
  
  if (!is.null(species_subgroup) && nzchar(species_subgroup)) {
    subgroup_enc <- gsub("/", "%252f", species_subgroup, fixed = TRUE)
    parts <- c(parts, utils::URLencode(subgroup_enc, reserved = TRUE))
  }
  
  parts <- c(parts,
             utils::URLencode(set_name, reserved = TRUE),
             version,
             fmt_final)
  
  paste(parts, collapse = "/")
}


# ---- Internal fetcher --------------------------------------------------------
.ogrdb_fetch_to_file <- function(url, dest_file) {
  try_one <- function(u) {
    tmp <- paste0(dest_file, ".tmp")
    status <- try(utils::download.file(u, tmp, quiet = TRUE, mode = "wb"), silent = TRUE)
    if (!inherits(status, "try-error") && file.exists(tmp) && file.info(tmp)$size > 0) {
      dir.create(dirname(dest_file), recursive = TRUE, showWarnings = FALSE)
      file.rename(tmp, dest_file)
      return(TRUE)
    }
    if (file.exists(tmp)) unlink(tmp)
    FALSE
  }
  
  # 1) Try API form
  ok <- try_one(url)
  if (!ok) {
    # 2) Try download form (swap prefix)
    fallback <- sub("/api/germline/set/", "/download_germline_set/", url, fixed = TRUE)
    ok <- try_one(fallback)
  }
  if (!ok) stop("Failed to download from OGRDB: ", url)
  dest_file
}

#' @title Download and Load Immune Receptor Germline Sequences from OGRDB (AIRR)
#' 
#' @description Downloads AIRR-compliant germline sets (or FASTA) from OGRDB and 
#' returns sequences as `DNAStringSet` (NUC) or attempts AA translation for V 
#' genes.
#' 
#' @param species Species string. Accepts "human"/"Homo sapiens"/"mouse"/"Mus 
#'    musculus".
#' @param locus Either a locus short code ("IGH","IGK","IGL", etc.) OR NULL if 
#'    you pass a `set_name` explicitly.
#' @param set_name Optional explicit OGRDB set name (e.g., "IGH_VDJ"). If 
#'    provided, overrides `locus`.
#' @param type "NUC" (default) or "PROT". PROT will translate V-gene CDS; only 
#'    supported for FASTA or AIRR records that include a valid CDS.
#' @param format "FASTA_GAPPED", "FASTA_UNGAPPED", or "AIRR". Default 
#'    "FASTA_GAPPED".
#' @param version "published" (default), "latest", or a specific revision 
#'    number as character/number.
#' @param species_subgroup Optional subgroup (e.g., a mouse strain like 
#'    "C57BL/6"). If it contains '/', OGRDB requires it encoded as "\%252f".
#' @param refresh If TRUE, redownload even if cached.
#' @param suppressMessages If TRUE, be quiet.
#' 
#' @examples
#' if (is_ogrdb_available()) {
#'   # Download human IGH nucleotide sequences (gapped FASTA)
#'   igh_nuc <- getOGRDB(species = "human",
#'                       locus   = "IGH",
#'                       type    = "NUC",
#'                       format  = "FASTA_GAPPED")
#'
#'   # Download human IGK sequences in AIRR JSON format
#'   igk_airr <- getOGRDB(species = "human",
#'                        locus   = "IGK",
#'                        type    = "NUC",
#'                        format  = "AIRR")
#'
#'   # Download human IGL sequences and translate to AA
#'   igl_prot <- getOGRDB(species = "human",
#'                        locus   = "IGL",
#'                        type    = "PROT",
#'                        format  = "FASTA_UNGAPPED")
#'
#'   # Example using an explicit set name (instead of locus)
#'   igh_explicit <- getOGRDB(species  = "human",
#'                            set_name = "IGH_VDJ",
#'                            type     = "NUC",
#'                            format   = "FASTA_GAPPED")
#' }
#' 
#' @return `DNAStringSet` for NUC; if `type="PROT"`, returns `AAStringSet` 
#'    where possible.
#' @importFrom jsonlite fromJSON
#' @export
getOGRDB <- function(species = "human",
                     locus = c("IGH","IGK","IGL"),
                     set_name = NULL,
                     type = c("NUC", "PROT"),
                     format = c("FASTA_GAPPED","FASTA_UNGAPPED","AIRR"),
                     version = c("published","latest"),
                     species_subgroup = NULL,
                     refresh = FALSE,
                     suppressMessages = FALSE) {
  
  .ensure_cache_dir()
  type   <- match.arg(type)
  format <- toupper(match.arg(format))
  version <- match.arg(version)
  
  # normalize species
  sp_key <- tolower(species)
  if (!sp_key %in% names(.ogrdb_species_map)) {
    stop("Species '", species, "' not recognized. Supported keys: ", 
         paste(names(.ogrdb_species_map), collapse = ", "))
  }
  species_api <- .ogrdb_species_map[[sp_key]]$api
  
  # resolve set name from locus if needed
  if (is.null(set_name)) {
    locus <- toupper(if (length(locus)) locus[1] else "")
    if (!nzchar(locus)) stop("Provide a locus (e.g., 'IGH') or an explicit set_name.")
    if (!locus %in% names(.ogrdb_default_set)) {
      stop("No default set mapping for locus '", locus,
           "'. Provide `set_name` explicitly (see OGRDB germline sets).")
    }
    set_name <- .ogrdb_default_set[[locus]]
  }
  
  # format mapping
  fmt_slug <- switch(toupper(format),
                     "FASTA_GAPPED"   = "gapped",
                     "FASTA_UNGAPPED" = "ungapped",
                     "AIRR"           = "airr"
  )
  
  # build cache path
  cache_dir <- file.path(.get_cache_dir(), species_api, "ogrdb")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  fname_core <- paste(c(species_api,
                        if (!is.null(species_subgroup) && nzchar(species_subgroup)) gsub("[/ ]", "_", species_subgroup) else NULL,
                        set_name, version, fmt_slug), collapse = "_")
  dest_file <- file.path(cache_dir, paste0(fname_core, if (fmt_slug == "airr") ".json" else ".fasta"))
  
  # download if needed
  if (refresh || !file.exists(dest_file) || file.info(dest_file)$size == 0) {
    url <- .ogrdb_build_url(species_api, set_name,
                            version = version,
                            format  = fmt_slug,
                            species_subgroup = species_subgroup,
                            human_reference = TRUE,   # keeps _ex on format for Human
                            use_api = TRUE)
    if (!suppressMessages) message("-> Downloading OGRDB: ", url)
    .ogrdb_fetch_to_file(url, dest_file)
    .update_cache_metadata(species = species_api, gene = set_name)
  } else if (!suppressMessages) {
    message("-> Using cached OGRDB file: ", dest_file)
  }
  
  # load & return sequences
  if (fmt_slug %in% c("gapped","ungapped")) {
    dna <- Biostrings::readDNAStringSet(dest_file)
    if (type == "NUC") return(dna)
    # Best-effort translation for V-gene CDS (warn that this is naive)
    if (!suppressMessages)
      message("-> Translating nucleotide sequences to AA (PROT). This assumes in-frame CDS; verify for your use case.")
    aa <- Biostrings::translate(dna, if.fuzzy.codon = "X")
    return(aa)
  } else { # AIRR JSON
    txt <- readLines(dest_file, warn = FALSE)
    if (length(txt) == 0L) stop("AIRR JSON is empty: ", dest_file)
    obj <- jsonlite::fromJSON(paste(txt, collapse = "\n"), simplifyVector = FALSE)
    
    gs <- NULL
    if (!is.null(obj$germline_set))            gs <- obj$germline_set
    else if (!is.null(obj$GermlineSet))       gs <- if (is.list(obj$GermlineSet) && !is.null(obj$GermlineSet[[1]])) obj$GermlineSet[[1]] else obj$GermlineSet
    else if (!is.null(obj[[1]]$germline_set)) gs <- obj[[1]]$germline_set
    else if (!is.null(obj[[1]]$GermlineSet))  gs <- obj[[1]]$GermlineSet
    
    # pull allele descriptions from either the set or top level
    alleles <- NULL
    if (!is.null(gs$allele_descriptions)) alleles <- gs$allele_descriptions
    else if (!is.null(obj$allele_descriptions)) alleles <- obj$allele_descriptions
    
    if (is.null(alleles) || length(alleles) == 0) {
      stop("AIRR JSON schema not recognized: missing allele_descriptions. ",
           "Top-level fields: ", paste(names(obj), collapse = ", "))
    }
    
    # extract sequences and labels with multiple fallbacks
    get_seq <- function(a) {
      a$sequence %||% a$sequence_ungapped %||% a$sequence_gapped %||% NA_character_
    }
    get_name <- function(a) {
      a$allele_description_name %||% a$label %||% a$name %||%
        a$allele_id %||% a$sequence_id %||% a$allele_description_id %||% ""
    }
    
    seqs  <- vapply(alleles, get_seq,  character(1), USE.NAMES = FALSE)
    names <- vapply(alleles, get_name, character(1), USE.NAMES = FALSE)
    
    keep <- !is.na(seqs) & nzchar(seqs)
    dna  <- Biostrings::DNAStringSet(seqs[keep])
    names(dna) <- names[keep]
    
    if (type == "NUC") return(dna)
    if (!suppressMessages) message("-> Translating AIRR JSON nucleotide sequences to AA (PROT).")
    aa <- Biostrings::translate(dna, if.fuzzy.codon = "X")
    return(aa)
    
  }
}

`%||%` <- function(x, y) if (is.null(x)) y else x

#' @title Load Cached OGRDB Sequences
#' @inheritParams getOGRDB
#' @description Loads sequences from the local cache without attempting to download.
#' This function relies on `getOGRDB(refresh = FALSE)`. If the data is not found
#' in the cache, it will be downloaded unless an internet connection is unavailable.
#' 
#' @examples
#' 
#' if (is_ogrdb_available()) {
#'   # First, ensure the file is cached
#'   getOGRDB(species = "human", locus = "IGH",
#'            type = "NUC", format = "FASTA_GAPPED",
#'            suppressMessages = TRUE)
#'
#'   # Now load from cache only
#'   igh_cached <- loadOGRDB(species = "human",
#'                           locus   = "IGH",
#'                           type    = "NUC",
#'                           format  = "FASTA_GAPPED")
#' }
#' @export
loadOGRDB <- function(species = "human", locus = c("IGH","IGK","IGL"),
                      set_name = NULL, type = c("NUC","PROT"),
                      format = c("FASTA_GAPPED","FASTA_UNGAPPED","AIRR"),
                      version = c("published","latest"),
                      species_subgroup = NULL, suppressMessages = FALSE) {
  getOGRDB(species, locus, set_name, type, format, version,
           species_subgroup, refresh = FALSE, suppressMessages = suppressMessages)
}

#' @title Force Re-download of OGRDB Sequences
#' @description A convenience wrapper for `getOGRDB(..., refresh = TRUE)` to 
#' ensure that the local cache is updated with the latest versions of the 
#' requested sequences.
#
#' @inheritParams getOGRDB
#' @examples
#' if (is_ogrdb_available()) {
#'   # Force a re-download of the human IGK sequences
#'   igk_fresh <- refreshOGRDB(species = "human",
#'                             locus   = "IGK",
#'                             type    = "NUC",
#'                             format  = "FASTA_GAPPED")
#' }
#' @export
refreshOGRDB <- function(species = "human", 
                         locus = c("IGH","IGK","IGL"),
                         set_name = NULL, 
                         type = c("NUC","PROT"),
                         format = c("FASTA_GAPPED","FASTA_UNGAPPED","AIRR"),
                         version = c("published","latest"),
                         species_subgroup = NULL, 
                         suppressMessages = FALSE) {
  getOGRDB(species, locus, set_name, type, format, version,
           species_subgroup, refresh = TRUE, suppressMessages = suppressMessages)
}

#' @title List OGRDB Datasets in Local Cache
#' @description Scans the cache directory and returns a list of available datasets.
#' 
#' @examples
#' if (is_ogrdb_available()) {
#'   # List cached OGRDB files
#'   cached_files <- listOGRDB()
#'   head(cached_files)
#' }
#' @export
listOGRDB <- function() {
  cache_dir <- .get_cache_dir()
  if (!dir.exists(cache_dir)) return(character(0))
  list.files(file.path(cache_dir, "*", "ogrdb"), recursive = TRUE, full.names = TRUE)
}
