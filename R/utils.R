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
    response <- httr::HEAD("https://www.imgt.org", timeout = 5)
    # Check for a successful status code (2xx)
    httr::status_code(response) >= 200 && httr::status_code(response) < 300
  }, error = function(e) {
    # If any error occurs (e.g., timeout, DNS issue), assume it's unavailable
    FALSE
  })
}
