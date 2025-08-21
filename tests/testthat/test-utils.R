test_that("is_imgt_available() returns TRUE on 2xx status", {
  with_mocked_bindings(
    {
      expect_true(is_imgt_available())
    },
    # Patch httr namespace (use unqualified names + .package = "httr")
    HEAD        = function(url, ...) "resp",
    status_code = function(resp) 200L,
    .package    = "httr"
  )
  
  with_mocked_bindings(
    {
      expect_true(is_imgt_available())  # 204 is 2xx
    },
    HEAD        = function(url, ...) "resp",
    status_code = function(resp) 204L,
    .package    = "httr"
  )
})

# tests/testthat/test-utils.R

test_that("is_imgt_available() returns FALSE on non-2xx status and on error", {
  # 500 -> FALSE
  with_mocked_bindings(
    {
      expect_false(is_imgt_available())
    },
    HEAD        = function(url, ...) "resp",
    status_code = function(resp) 500L,
    .package    = "httr"
  )
  
  # Redirect 301 is not 2xx -> FALSE
  with_mocked_bindings(
    {
      expect_false(is_imgt_available())
    },
    HEAD        = function(url, ...) "resp",
    status_code = function(resp) 301L,
    .package    = "httr"
  )
  
  # Error path -> FALSE
  with_mocked_bindings(
    {
      expect_false(is_imgt_available())
    },
    HEAD        = function(url, ...) stop("timeout"),
    # status_code won't be called
    .package    = "httr"
  )
})