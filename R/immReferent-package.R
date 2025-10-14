#' immReferent: An Interface for Immune Receptor and HLA Gene Reference Data
#'
#' @description
#' **immReferent** provides a stable, reproducible, and lightweight interface to
#' reference sequences for immune receptors (TCR/BCR) and HLA genes sourced from
#' IMGT, IPD-IMGT/HLA, and the AIRR-C's OGRDB. It centralizes downloading,
#' caching, and querying of curated nucleotide and protein sequences, and plays a
#' foundational role in computational immunology workflows.
#'
#' @details
#' The package is designed as a common reference layer across immunoinformatics tools,
#' ensuring consistent provenance and offline reproducibility via caching.
#'
#' **Core functionality**
#' \itemize{
#'   \item Download and parse receptor and HLA sequences from IMGT and OGRDB.
#'   \item Local caching to support offline, reproducible analysis.
#'   \item Query by gene, allele, species, locus, and sequence type/format.
#'   \item Interoperability with Bioconductor classes such as
#'         \code{Biostrings::DNAStringSet} and \code{Biostrings::AAStringSet}.
#' }
#'
#' **Supported data sources**
#' \itemize{
#'   \item IMGT — The international ImMunoGeneTics information system:
#'         \url{https://www.imgt.org/}
#'   \item IPD-IMGT/HLA — The HLA Database:
#'         \url{https://www.ebi.ac.uk/ipd/imgt/hla/}
#'   \item OGRDB — Open Germline Receptor Database (AIRR-C):
#'         \url{https://ogrdb.airr-community.org/}
#' }
#'
##' \strong{Getting started}
#' \preformatted{
#' browseVignettes("immReferent")
#' }
#'
#' @section Attribution and Licensing:
#' Data obtained from IMGT and OGRDB must be cited according to their terms.
#' IMGT data are distributed under a
#' \href{https://creativecommons.org/licenses/by-nc-nd/4.0/}{CC BY-NC-ND 4.0}
#' license. Proper attribution is required, and derivative or commercial use is
#' restricted per IMGT policy. Always review the current licensing and citation
#' requirements of each resource prior to use.
#'
#'#' @seealso
#' \url{https://github.com/BorchLab/immReferent} \cr
#' \url{https://github.com/BorchLab/Ibex/immReferent}
#' 
#' @keywords package
#' @md
#' @name immReferent-package
#' @aliases immReferent immReferent-package
#' @docType package
"_PACKAGE"
