#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

#' immReferent: A package for downloading and managing immune receptor and HLA gene reference data.
#'
#' The immReferent package provides a consistent and reliable interface to download,
#' cache, and load reference sequences from the IMGT and IPD-IMGT/HLA databases.
#' It is designed to be a core dependency in immunogenomics analysis pipelines.
#'
#' @section Main functions:
#' \itemize{
#'   \item \code{\link{getIMGT}}: The main function to download and load data with caching.
#'   \item \code{\link{loadIMGT}}: Loads data strictly from the cache.
#'   \item \code{\link{listIMGT}}: Lists datasets available in the local cache.
#'   \item \code{\link{refreshIMGT}}: Forces a re-download to update cached data.
#' }
#'
#' @docType package
#' @name immReferent
NULL
