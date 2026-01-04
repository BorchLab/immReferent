#' @title Check if IMGT Website is Available
#'
#' @description Sends a lightweight HEAD request to the main IMGT page to check
#' if the service is online and accessible. This function is used to
#' conditionally run examples and tests that require an internet connection.
#'
#' @return A logical value: \code{TRUE} if the IMGT website is accessible,
#'   \code{FALSE} otherwise.
#'
#' @export
#' @importFrom httr HEAD status_code
#' @seealso
#' \code{\link{is_ogrdb_available}} for checking OGRDB availability
#'
#' \code{\link{getIMGT}} which uses this function
#'
#' @examples
#' is_imgt_available()
is_imgt_available <- function() {
  tryCatch({
    # Use httr::HEAD for a lightweight request and set a timeout
    response <- httr::HEAD("https://www.imgt.org", timeout = 2)
    httr::status_code(response) >= 200 && httr::status_code(response) < 300
  }, error = function(e) {
    FALSE
  })
}

#' @title Check if OGRDB Website is Available
#'
#' @description Sends a lightweight HEAD request to the OGRDB API to check if
#' the service is online and accessible. This function is used to conditionally
#' run examples and tests that require an internet connection.
#'
#' @return A logical value: \code{TRUE} if the OGRDB website is accessible,
#'   \code{FALSE} otherwise.
#'
#' @export
#' @importFrom httr HEAD status_code
#' @seealso
#' \code{\link{is_imgt_available}} for checking IMGT availability
#'
#' \code{\link{getOGRDB}} which uses this function
#'
#' @examples
#' is_ogrdb_available()
is_ogrdb_available <- function() {
  tryCatch({
    response <- httr::HEAD("https://ogrdb.airr-community.org/api", httr::timeout(2))
    httr::status_code(response) >= 200 && httr::status_code(response) < 300
  }, error = function(e) {
    FALSE
  })
}
