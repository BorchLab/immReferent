# tests/testthat/test-download.R

test_that(".fetch_imgt_query writes cleaned FASTA when HTML has two <pre> blocks", {
  tmp <- withr::local_tempfile(fileext = ".fasta")
  
  html <- paste0(
    "<html><body>",
    "<pre>ignored</pre>",
    "<pre>>seq1 | Homo sapiens\nATGC\n</pre>",
    "</body></html>"
  )
  
  with_mocked_bindings(
    {
      .fetch_imgt_query("http://fake", tmp, "human")
    },
    # Patch the httr namespace (CRITICAL): use unqualified names + .package="httr"
    RETRY       = function(...) "resp",
    status_code = function(resp) 200L,
    content     = function(resp, as, encoding) html,
    .package    = "httr"
  )
  
  expect_true(file.exists(tmp))
  txt <- readLines(tmp, warn = FALSE)
  expect_true(any(grepl("Homo_sapiens", txt), na.rm = TRUE))
  expect_false(any(grepl("Homo sapiens",  txt), na.rm = TRUE))
})

test_that(".fetch_imgt_query warns and does not write on HTTP >= 400", {
  tmp <- withr::local_tempfile(fileext = ".fasta")
  
  expect_warning(
    with_mocked_bindings(
      {
        .fetch_imgt_query("http://fake404", tmp, "human")
      },
      RETRY       = function(...) "resp",
      status_code = function(resp) 404L,
      content     = function(...) "",
      .package    = "httr"
    ),
    regexp = "Failed to download"
  )
  
  expect_false(file.exists(tmp))
})

test_that(".fetch_imgt_query warns when <pre> content is missing (length < 2)", {
  tmp <- withr::local_tempfile(fileext = ".fasta")
  html <- "<html><body><pre>only one</pre></body></html>"
  
  expect_warning(
    with_mocked_bindings(
      {
        .fetch_imgt_query("http://fake", tmp, "human")
      },
      RETRY       = function(...) "resp",
      status_code = function(resp) 200L,
      content     = function(resp, as, encoding) html,
      .package    = "httr"
    ),
    regexp = "Could not find FASTA content"
  )
  
  expect_false(file.exists(tmp))
})

test_that(".fetch_imgt_query warns when second <pre> is empty", {
  tmp <- withr::local_tempfile(fileext = ".fasta")
  html <- "<html><body><pre>ignored</pre><pre></pre></body></html>"
  
  expect_warning(
    with_mocked_bindings(
      {
        .fetch_imgt_query("http://fake", tmp, "human")
      },
      RETRY       = function(...) "resp",
      status_code = function(resp) 200L,
      content     = function(resp, as, encoding) html,
      .package    = "httr"
    ),
    regexp = "Empty FASTA content"
  )
  
  expect_false(file.exists(tmp))
})

test_that(".fetch_hla_file validates file_type", {
  tmp <- withr::local_tempfile(fileext = ".fasta")
  
  expect_error(
    .fetch_hla_file("dna", tmp),
    "file_type must be one of 'nuc' or 'prot'"
  )
  
  with_mocked_bindings(
    {
      expect_silent(.fetch_hla_file("nuc",  tmp))
      expect_silent(.fetch_hla_file("prot", tmp))
    },
    RETRY       = function(..., times, pause_base, pause_cap) "resp",
    status_code = function(resp) 200L,
    .package    = "httr"
  )
})

test_that(".fetch_hla_file writes file (simulated) and warns on HTTP >= 400", {
  tmp_ok  <- withr::local_tempfile(fileext = ".fasta")
  tmp_bad <- withr::local_tempfile(fileext = ".fasta")
  
  # Success path: simulate write_disk by writing in RETRY mock
  with_mocked_bindings(
    {
      .fetch_hla_file("nuc", tmp_ok)
    },
    RETRY = function(method, url, times, pause_base, pause_cap, ...) {
      writeLines(">HLA-A*01:01\nATGC", tmp_ok)
      "resp"
    },
    status_code = function(resp) 200L,
    .package = "httr"
  )
  expect_true(file.exists(tmp_ok))
  expect_match(readLines(tmp_ok, warn = FALSE)[1], "^>HLA")
  
  # Failure path: HTTP 500 -> warn; no file content written
  expect_warning(
    with_mocked_bindings(
      {
        .fetch_hla_file("nuc", tmp_bad)
      },
      RETRY       = function(...) "resp",
      status_code = function(resp) 500L,
      .package    = "httr"
    ),
    regexp = "Failed to download"
  )
  expect_false(file.exists(tmp_bad))
})
