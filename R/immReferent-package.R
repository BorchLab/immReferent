#' immReferent: An Interface for Immune Receptor and HLA Gene Reference Data
#'
#' @description
#' \strong{immReferent} provides a stable, reproducible, and lightweight
#' interface to reference sequences for immune receptors (TCR/BCR) and HLA genes
#' sourced from IMGT, IPD-IMGT/HLA, and the AIRR-C's OGRDB. It centralizes
#' downloading, caching, and querying of curated nucleotide and protein
#' sequences, and plays a foundational role in computational immunology
#' workflows.
#'
#' @details
#' The package is designed as a common reference layer across immunoinformatics
#' tools, ensuring consistent provenance and offline reproducibility via
#' caching.
#'
#' \strong{Core functionality}
#' \itemize{
#'   \item Download and parse receptor and HLA sequences from IMGT and OGRDB
#'   \item Local caching to support offline, reproducible analysis
#'   \item Query by gene, allele, species, locus, and sequence type/format
#'   \item Export to popular analysis tools (MiXCR, TRUST4, Cell Ranger, IgBLAST)
#'   \item Interoperability with Bioconductor classes such as
#'         \code{\link[Biostrings]{DNAStringSet}} and
#'         \code{\link[Biostrings]{AAStringSet}}
#' }
#'
#' \strong{Data retrieval functions}
#' \itemize{
#'   \item \code{\link{getIMGT}}: Download sequences from IMGT
#'   \item \code{\link{getOGRDB}}: Download sequences from OGRDB
#'   \item \code{\link{loadIMGT}}, \code{\link{loadOGRDB}}: Load cached sequences
#'   \item \code{\link{refreshIMGT}}, \code{\link{refreshOGRDB}}: Force re-download
#' }
#'
#' \strong{Export functions}
#' \itemize{
#'   \item \code{\link{exportMiXCR}}: Export for MiXCR analysis
#'   \item \code{\link{exportTRUST4}}: Export for TRUST4 analysis
#'   \item \code{\link{exportCellRanger}}: Export for 10x Cell Ranger VDJ
#'   \item \code{\link{exportIgBLAST}}: Export for IgBLAST analysis
#' }
#'
#' \strong{Supported data sources}
#' \itemize{
#'   \item IMGT: The international ImMunoGeneTics information system
#'         (\url{https://www.imgt.org/})
#'   \item IPD-IMGT/HLA: The HLA Database
#'         (\url{https://www.ebi.ac.uk/ipd/imgt/hla/})
#'   \item OGRDB: Open Germline Receptor Database (AIRR-C)
#'         (\url{https://ogrdb.airr-community.org/})
#' }
#'
#' \strong{Getting started}
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
#' @seealso
#' \url{https://github.com/BorchLab/immReferent}
#'
#' @keywords package
#' @name immReferent-package
#' @aliases immReferent immReferent-package
#' @docType package
"_PACKAGE"
