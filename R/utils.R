#' @title Check if IMGT Website is Available
#' @description This function sends a lightweight HEAD request to the main IMGT page to
#' check if the service is online and accessible.
#' @return A logical value: `TRUE` if the IMGT website is accessible, `FALSE` otherwise.
#' @export
#' @importFrom httr HEAD status_code
#' @examples
#'   is_imgt_available()
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
#' @description This function sends a lightweight HEAD request to the main OGRDB page to
#' check if the service is online and accessible.
#' @return A logical value: `TRUE` if the IMGT website is accessible, `FALSE` otherwise.
#' @export
#' @importFrom httr HEAD status_code
#' @examples
#'   is_imgt_available()
is_ogrdb_available <- function() {
  tryCatch({
    response <- httr::HEAD("https://ogrdb.airr-community.org/api", httr::timeout(2))
    httr::status_code(response) >= 200 && httr::status_code(response) < 300
  }, error = function(e) {
    FALSE
  })
}
