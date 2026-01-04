# This file contains export functions for creating reference databases
# compatible with various immune repertoire analysis tools.

#' @title Export Reference Sequences to MiXCR Format
#' @description Exports a DNAStringSet or AAStringSet to FASTA files formatted
#' for use with MiXCR's `buildLibrary` command. The function creates separate
#' FASTA files for V, D, J, and C gene segments.
#'
#' @param sequences A `DNAStringSet` or `AAStringSet` object containing immune
#'        receptor sequences. Sequence names should follow IMGT nomenclature
#'        (e.g., "IGHV1-2*01", "TRBJ2-1*01").
#' @param output_dir The directory where output files will be written.
#' @param chain The chain type for the output files. One of "IGH", "IGK", "IGL",
#'        "TRA", "TRB", "TRD", or "TRG".
#'
#' @details
#' MiXCR expects FASTA files with simple headers containing only the gene name.
#' The function filters sequences by gene type (V, D, J, C) based on the gene

#' name pattern and writes separate files for each segment type.
#'
#' Output files follow the naming convention: `v-genes.<chain>.fasta`,
#' `d-genes.<chain>.fasta`, `j-genes.<chain>.fasta`, `c-genes.<chain>.fasta`.
#'
#' @return A named list containing the paths to the created files, invisibly.
#' @export
#' @seealso \url{https://mixcr.com/mixcr/guides/create-custom-library/}
#' @examples
#' # Create a small example DNAStringSet
#' seqs <- Biostrings::DNAStringSet(c(
#'   "ATGCGATCGATCGATCG",
#'   "ATGCGATCGATCG",
#'   "ATGCGATC",
#'   "ATGCGATCGATCGATCGATCG"
#' ))
#' names(seqs) <- c("IGHV1-2*01", "IGHD1-1*01", "IGHJ1*01", "IGHC*01")
#'
#' # Export to temporary directory
#' output_dir <- tempdir()
#' files <- exportMiXCR(seqs, output_dir, chain = "IGH")
#' print(files)
#'
#' # Clean up
#' unlink(unlist(files))
exportMiXCR <- function(sequences, output_dir, chain = c("IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG")) {
  chain <- match.arg(chain)

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Get sequence names
  seq_names <- names(sequences)
  if (is.null(seq_names)) {
    stop("Sequences must have names following IMGT nomenclature.", call. = FALSE)
  }

  # Categorize sequences by segment type
  # Pattern matches: IGHV, TRBV, etc. for V; IGHD, TRBD for D; IGHJ, TRBJ for J; IGHC, TRBC for C
  v_idx <- grep(paste0("^", chain, "V"), seq_names)
  d_idx <- grep(paste0("^", chain, "D"), seq_names)
  j_idx <- grep(paste0("^", chain, "J"), seq_names)
  c_idx <- grep(paste0("^", chain, "C"), seq_names)

  output_files <- list()

  # Helper function to write FASTA with MiXCR-compatible headers
  .write_mixcr_fasta <- function(seqs, path) {
    if (length(seqs) == 0) return(NULL)
    # MiXCR expects simple gene names as headers
    # Extract just the gene name (everything before any extra annotation)
    simple_names <- sub("\\|.*$", "", names(seqs))
    names(seqs) <- simple_names
    Biostrings::writeXStringSet(seqs, path)
    return(path)
  }

  # Write V genes
  if (length(v_idx) > 0) {
    v_path <- file.path(output_dir, paste0("v-genes.", chain, ".fasta"))
    output_files$v_genes <- .write_mixcr_fasta(sequences[v_idx], v_path)
  }

  # Write D genes
  if (length(d_idx) > 0) {
    d_path <- file.path(output_dir, paste0("d-genes.", chain, ".fasta"))
    output_files$d_genes <- .write_mixcr_fasta(sequences[d_idx], d_path)
  }

  # Write J genes
  if (length(j_idx) > 0) {
    j_path <- file.path(output_dir, paste0("j-genes.", chain, ".fasta"))
    output_files$j_genes <- .write_mixcr_fasta(sequences[j_idx], j_path)
  }

  # Write C genes
  if (length(c_idx) > 0) {
    c_path <- file.path(output_dir, paste0("c-genes.", chain, ".fasta"))
    output_files$c_genes <- .write_mixcr_fasta(sequences[c_idx], c_path)
  }

  if (length(output_files) == 0) {
    warning("No sequences matching chain '", chain, "' were found.", call. = FALSE)
  }

  invisible(output_files)
}


#' @title Export Reference Sequences to TRUST4 Format
#' @description Exports a DNAStringSet to a FASTA file formatted for use with
#' TRUST4. The output follows the format produced by TRUST4's `BuildImgtAnnot.pl`
#' script.
#'
#' @param sequences A `DNAStringSet` object containing immune receptor sequences.
#'        Sequence names should follow IMGT nomenclature (e.g., "IGHV1-2*01").
#' @param output_file The path to the output FASTA file.
#' @param include_constant Logical. If `TRUE`, include constant region sequences.
#'        TRUST4's IMGT+C.fa file includes constant regions. Default is `TRUE`.
#'
#' @details
#' TRUST4 expects FASTA files with headers containing only the allele name
#' (e.g., ">IGHV1-2*01"). The function reformats sequence headers to match
#' the output of TRUST4's `BuildImgtAnnot.pl` script.
#'
#' TRUST4 uses this reference for the `--ref` parameter in its analysis pipeline.
#'
#' @return The path to the created file, invisibly.
#' @export
#' @seealso \url{https://github.com/liulab-dfci/TRUST4}
#' @examples
#' # Create a small example DNAStringSet
#' seqs <- Biostrings::DNAStringSet(c(
#'   "ATGCGATCGATCGATCG",
#'   "ATGCGATCGATCG",
#'   "ATGCGATC"
#' ))
#' names(seqs) <- c("IGHV1-2*01", "IGHJ1*01", "IGHC*01")
#'
#' # Export to temporary file
#' output_file <- tempfile(fileext = ".fa")
#' exportTRUST4(seqs, output_file)
#'
#' # View the result
#' cat(readLines(output_file), sep = "\n")
#'
#' # Clean up
#' unlink(output_file)
exportTRUST4 <- function(sequences, output_file, include_constant = TRUE) {
  if (!inherits(sequences, "DNAStringSet")) {
    stop("sequences must be a DNAStringSet object.", call. = FALSE)
  }

  seq_names <- names(sequences)
  if (is.null(seq_names)) {
    stop("Sequences must have names following IMGT nomenclature.", call. = FALSE)
  }

  # Filter out constant regions if requested
  if (!include_constant) {
    # Keep only sequences with V, D, or J in the name pattern
    keep_idx <- grep("(IG[HKL]|TR[ABDG])[VDJ]", seq_names)
    sequences <- sequences[keep_idx]
    seq_names <- names(sequences)
  }

  if (length(sequences) == 0) {
    stop("No valid sequences to export.", call. = FALSE)
  }

  # TRUST4's BuildImgtAnnot.pl extracts only the allele name from IMGT headers
  # Format: gene*allele (e.g., IGHV1-2*01)
  # Remove any additional annotation after the allele designation
  simple_names <- sub("\\|.*$", "", seq_names)
  simple_names <- sub(" .*$", "", simple_names)
  names(sequences) <- simple_names

  # Ensure output directory exists
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir) && output_dir != ".") {
    dir.create(output_dir, recursive = TRUE)
  }

  Biostrings::writeXStringSet(sequences, output_file)

  invisible(output_file)
}


#' @title Export Reference Sequences to Cell Ranger VDJ Format
#' @description Exports a DNAStringSet to FASTA format suitable for creating
#' a custom Cell Ranger VDJ reference. The function generates a FASTA file
#' with properly formatted headers for use with `cellranger mkvdjref`.
#'
#' @param sequences A `DNAStringSet` object containing immune receptor sequences.
#'        Sequence names should follow IMGT nomenclature (e.g., "IGHV1-2*01").
#' @param output_file The path to the output FASTA file.
#' @param gene_type The type of gene region. One of "V", "D", "J", or "C".
#'        If NULL (default), the function will attempt to infer the type from
#'        sequence names.
#'
#' @details
#' Cell Ranger's `mkvdjref` command expects FASTA files with specific header
#' formats. This function creates a FASTA file that can be used as input
#' to build a custom VDJ reference.
#'
#' Note: For a complete Cell Ranger VDJ reference, you also need a GTF file
#' with gene annotations. This function only creates the FASTA component.
#'
#' @return The path to the created file, invisibly.
#' @export
#' @seealso \url{https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-5p-references}
#' @examples
#' # Create a small example DNAStringSet
#' seqs <- Biostrings::DNAStringSet(c(
#'   "ATGCGATCGATCGATCG",
#'   "ATGCGATCGATCG"
#' ))
#' names(seqs) <- c("IGHV1-2*01", "IGHV1-3*01")
#'
#' # Export to temporary file
#' output_file <- tempfile(fileext = ".fa")
#' exportCellRanger(seqs, output_file)
#'
#' # View the result
#' cat(readLines(output_file), sep = "\n")
#'
#' # Clean up
#' unlink(output_file)
exportCellRanger <- function(sequences, output_file, gene_type = NULL) {
  if (!inherits(sequences, "DNAStringSet")) {
    stop("sequences must be a DNAStringSet object.", call. = FALSE)
  }

  seq_names <- names(sequences)
  if (is.null(seq_names)) {
    stop("Sequences must have names following IMGT nomenclature.", call. = FALSE)
  }

  if (length(sequences) == 0) {
    stop("No sequences to export.", call. = FALSE)
  }

  # Cell Ranger expects clean gene names
  # Format headers as: >gene_name
  simple_names <- sub("\\|.*$", "", seq_names)
  simple_names <- sub(" .*$", "", simple_names)
  names(sequences) <- simple_names

  # Ensure output directory exists
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir) && output_dir != ".") {
    dir.create(output_dir, recursive = TRUE)
  }

  Biostrings::writeXStringSet(sequences, output_file)

  invisible(output_file)
}


#' @title Export Reference Sequences to IgBLAST Format
#' @description Exports a DNAStringSet to FASTA files formatted for use with
#' IgBLAST. The function creates separate FASTA files for V, D, and J gene
#' segments with simplified headers compatible with IgBLAST's requirements.
#'
#' @param sequences A `DNAStringSet` object containing immune receptor sequences.
#'        Sequence names should follow IMGT nomenclature (e.g., "IGHV1-2*01").
#' @param output_dir The directory where output files will be written.
#' @param organism The organism name for the output files. Used in file naming.
#'        Default is "custom".
#' @param receptor_type The receptor type. One of "ig" for immunoglobulin or
#'        "tcr" for T-cell receptor. Default is "ig".
#'
#' @details
#' IgBLAST requires FASTA files with simplified headers containing only the
#' gene/allele name. This function mimics the output of IgBLAST's
#' `edit_imgt_file.pl` script, which truncates IMGT headers to keep only
#' the allele designation.
#'
#' Output files follow the naming convention used by IgBLAST:
#' `<organism>_<receptor_type>_v.fasta`, `<organism>_<receptor_type>_d.fasta`,
#' `<organism>_<receptor_type>_j.fasta`.
#'
#' After exporting, use `makeblastdb` with the `-parse_seqids` flag to create
#' the BLAST database:
#' ```
#' makeblastdb -parse_seqids -dbtype nucl -in <fasta_file> -out <db_name>
#' ```
#'
#' @return A named list containing the paths to the created files, invisibly.
#' @export
#' @seealso \url{https://ncbi.github.io/igblast/}
#' @examples
#' # Create a small example DNAStringSet
#' seqs <- Biostrings::DNAStringSet(c(
#'   "ATGCGATCGATCGATCG",
#'   "ATGCGATCGATCG",
#'   "ATGCGATC"
#' ))
#' names(seqs) <- c("IGHV1-2*01", "IGHD1-1*01", "IGHJ1*01")
#'
#' # Export to temporary directory
#' output_dir <- tempdir()
#' files <- exportIgBLAST(seqs, output_dir, organism = "human", receptor_type = "ig")
#' print(files)
#'
#' # Clean up
#' unlink(unlist(files))
exportIgBLAST <- function(sequences, output_dir, organism = "custom", receptor_type = c("ig", "tcr")) {
  receptor_type <- match.arg(receptor_type)

  if (!inherits(sequences, "DNAStringSet")) {
    stop("sequences must be a DNAStringSet object.", call. = FALSE)
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  seq_names <- names(sequences)
  if (is.null(seq_names)) {
    stop("Sequences must have names following IMGT nomenclature.", call. = FALSE)
  }

  # IgBLAST's edit_imgt_file.pl simplifies headers to just the allele name
  # Categorize by segment type
  if (receptor_type == "ig") {
    v_idx <- grep("^IG[HKL]V", seq_names)
    d_idx <- grep("^IGHD", seq_names)
    j_idx <- grep("^IG[HKL]J", seq_names)
  } else {
    v_idx <- grep("^TR[ABDG]V", seq_names)
    d_idx <- grep("^TR[BD]D", seq_names)
    j_idx <- grep("^TR[ABDG]J", seq_names)
  }

  output_files <- list()

  # Helper function to write FASTA with IgBLAST-compatible headers
  .write_igblast_fasta <- function(seqs, path) {
    if (length(seqs) == 0) return(NULL)
    # Simplify names to just allele designation
    simple_names <- sub("\\|.*$", "", names(seqs))
    simple_names <- sub(" .*$", "", simple_names)
    names(seqs) <- simple_names
    Biostrings::writeXStringSet(seqs, path)
    return(path)
  }

  # Write V genes
  if (length(v_idx) > 0) {
    v_path <- file.path(output_dir, paste0(organism, "_", receptor_type, "_v.fasta"))
    output_files$v_genes <- .write_igblast_fasta(sequences[v_idx], v_path)
  }

  # Write D genes
  if (length(d_idx) > 0) {
    d_path <- file.path(output_dir, paste0(organism, "_", receptor_type, "_d.fasta"))
    output_files$d_genes <- .write_igblast_fasta(sequences[d_idx], d_path)
  }

  # Write J genes
  if (length(j_idx) > 0) {
    j_path <- file.path(output_dir, paste0(organism, "_", receptor_type, "_j.fasta"))
    output_files$j_genes <- .write_igblast_fasta(sequences[j_idx], j_path)
  }

  if (length(output_files) == 0) {
    warning("No sequences matching receptor type '", receptor_type, "' were found.", call. = FALSE)
  }

  invisible(output_files)
}
