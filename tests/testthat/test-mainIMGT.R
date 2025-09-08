# tests/testthat/test-mainIMGT.R

testthat::test_that(".show_license_message shows once and respects suppression", {
  # Reset option
  old <- getOption("immReferent.license.shown", NULL)
  withr::defer(options(immReferent.license.shown = old))
  
  options(immReferent.license.shown = NULL)
  
  # First call: should message and set the option
  expect_snapshot(
    .show_license_message(suppress = FALSE),
    transform = function(x) gsub("https?://\\S+", "<url>", x) # keep snapshot stable
  )
  expect_true(isTRUE(getOption("immReferent.license.shown")))
  
  # Second call: no messages (already shown)
  expect_silent(.show_license_message(suppress = FALSE))
  
  # Suppressed: no messages, does not alter option
  options(immReferent.license.shown = NULL)
  expect_silent(.show_license_message(suppress = TRUE))
  expect_null(getOption("immReferent.license.shown", NULL))
})

# Utility: small DNA / AA FASTA writers for simulated downloads
.write_dna_fasta <- function(path) {
  writeLines(c(">seq1", "ATGC", ">seq2", "ATGCGG"), path)
}
.write_aa_fasta <- function(path) {
  writeLines(c(">seqA", "MSTK", ">seqB", "VVVVAA"), path)
}

testthat::test_that("getIMGT() — HLA branch: downloads when missing, loads when cached", {
  cache <- withr::local_tempdir()
  # Prepare a place for HLA cache to be written
  hla_dir <- file.path(cache, "human", "hla")
  # We’ll mock the *package* internals, so set .package to your package name.
  # If your package is immReferent (as shown), keep "immReferent" below.
  pkg <- "immReferent"
  
  # 1) Fresh download (no cache yet) — NUC
  with_mocked_bindings(
    {
      res <- getIMGT(gene = "HLA", type = "NUC", suppressMessages = TRUE)
      expect_s4_class(res, "DNAStringSet")
      expect_identical(length(res) > 0, TRUE)
    },
    # point cache to temp dir and skip real ensure
    .get_cache_dir = function() cache,
    .ensure_cache_dir = function() NULL,
    # pretend internet is available
    is_imgt_available = function() TRUE,
    # simulate download: write DNA FASTA
    .fetch_hla_file = function(file_type, dest_file) {
      expect_true(file_type %in% c("nuc", "prot"))
      dir.create(dirname(dest_file), recursive = TRUE, showWarnings = FALSE)
      if (file_type == "nuc") .write_dna_fasta(dest_file) else .write_aa_fasta(dest_file)
      invisible(NULL)
    },
    .update_cache_metadata = function(species, gene) NULL,
    .package = pkg
  )
  
  # 2) Uses cache when present — PROT (do not call .fetch_hla_file)
  # Precreate the cached file
  dir.create(hla_dir, recursive = TRUE, showWarnings = FALSE)
  prot_path <- file.path(hla_dir, "hla_prot.fasta")
  .write_aa_fasta(prot_path)
  
  with_mocked_bindings(
    {
      expect_silent({
        res <- getIMGT(gene = "HLA", type = "PROT", suppressMessages = TRUE)
      })
      expect_s4_class(res, "AAStringSet")
      expect_identical(length(res) > 0, TRUE)
    },
    .get_cache_dir = function() cache,
    .ensure_cache_dir = function() NULL,
    # If code tried to download, we'd error to catch it
    is_imgt_available = function() stop("should not be called"),
    .fetch_hla_file = function(...) stop("should not be called"),
    .update_cache_metadata = function(...) stop("should not be called"),
    .package = pkg
  )
})

testthat::test_that("getIMGT() — HLA branch: errors when internet unavailable and cache missing", {
  cache <- withr::local_tempdir()
  pkg <- "immReferent"
  
  expect_error(
    with_mocked_bindings(
      {
        getIMGT(gene = "HLA", type = "NUC", suppressMessages = TRUE)
      },
      .get_cache_dir = function() cache,
      .ensure_cache_dir = function() NULL,
      is_imgt_available = function() FALSE,
      .fetch_hla_file = function(...) stop("no download expected"),
      .update_cache_metadata = function(...) NULL,
      .package = pkg
    ),
    "Internet connection unavailable\\."
  )
})

testthat::test_that("getIMGT() — TCR/BCR branch: downloads per mapping and loads combined Biostrings", {
  cache <- withr::local_tempdir()
  pkg <- "immReferent"
  
  # Minimal rules for a single gene 'IGHV' NUC mapping
  .fake_map <- data.frame(
    gene = c("IGHV"),
    type = c("NUC"),
    cache_subdir = c(file.path("bcr", "ighv")),
    filename_prefix = c("imgt_"),
    query_prefix = c("V-REGION+GENE+and+allele="),
    query_label = c(""),
    stringsAsFactors = FALSE
  )
  
  with_mocked_bindings(
    {
      res <- getIMGT(species = "human", gene = "IGHV", type = "NUC", suppressMessages = TRUE)
      expect_s4_class(res, "DNAStringSet")
      # Our mock writes two sequences
      expect_gte(length(res), 2L)
    },
    .get_cache_dir = function() cache,
    .ensure_cache_dir = function() NULL,
    is_imgt_available = function() TRUE,
    # supply only the rows above
    .imgt_db_map = .fake_map,
    .gene_groups = list(),               # no grouping
    # simulate TCR/BCR download writer
    .fetch_imgt_query = function(url, dest_file, species_key) {
      dir.create(dirname(dest_file), recursive = TRUE, showWarnings = FALSE)
      .write_dna_fasta(dest_file)
      invisible(NULL)
    },
    .update_cache_metadata = function(species, gene) NULL,
    .package = pkg
  )
})

testthat::test_that("getIMGT() — TCR/BCR branch: falls back to NUC-C rule when exact type missing", {
  cache <- withr::local_tempdir()
  pkg <- "immReferent"
  
  # Only NUC-C available; request PROT -> should fall back and still fetch
  .fake_map <- data.frame(
    gene = c("TRB"),
    type = c("NUC-C"),  # consolidation rule
    cache_subdir = c(file.path("tcr", "trb")),
    filename_prefix = c("imgt_"),
    query_prefix = c("V-REGION+GENE+and+allele="),
    query_label = c(""),
    stringsAsFactors = FALSE
  )
  
  with_mocked_bindings(
    {
      res <- getIMGT(species = "human", gene = "TRB", type = "PROT", suppressMessages = TRUE)
      # Even though we requested PROT, the file content we write can be AA
      expect_s4_class(res, "AAStringSet")
      expect_gte(length(res), 2L)
    },
    .get_cache_dir = function() cache,
    .ensure_cache_dir = function() NULL,
    is_imgt_available = function() TRUE,
    .imgt_db_map = .fake_map,
    .gene_groups = list(),
    .fetch_imgt_query = function(url, dest_file, species_key) {
      dir.create(dirname(dest_file), recursive = TRUE, showWarnings = FALSE)
      .write_aa_fasta(dest_file)
      invisible(NULL)
    },
    .update_cache_metadata = function(species, gene) NULL,
    .package = pkg
  )
})

testthat::test_that("getIMGT() — TCR/BCR branch: warns and returns NULL when no rules found", {
  cache <- withr::local_tempdir()
  pkg <- "immReferent"
  
  expect_warning(
    res <- with_mocked_bindings(
      {
        getIMGT(species = "human", gene = "FAKEGENE", type = "NUC", suppressMessages = TRUE)
      },
      .get_cache_dir = function() cache,
      .ensure_cache_dir = function() NULL,
      is_imgt_available = function() TRUE,
      .imgt_db_map = data.frame(
        gene = character(), type = character(), cache_subdir = character(),
        filename_prefix = character(), query_prefix = character(), query_label = character(),
        stringsAsFactors = FALSE
      ),
      .gene_groups = list(),
      .fetch_imgt_query = function(...) stop("should not be called"),
      .update_cache_metadata = function(...) NULL,
      .package = pkg
    ),
    "No data could be found or downloaded for the request\\."
  )
  expect_null(res)
})

testthat::test_that("getIMGT() — TCR/BCR branch: errors on unsupported species", {
  cache <- withr::local_tempdir()
  pkg <- "immReferent"
  
  expect_error(
    with_mocked_bindings(
      {
        getIMGT(species = "platypus", gene = "IGHV", type = "NUC", suppressMessages = TRUE)
      },
      .get_cache_dir = function() cache,
      .ensure_cache_dir = function() NULL,
      # Use the package's real .species_map for coverage if available; otherwise stub a small one:
      .species_map = list(human = "Homo+sapiens"),
      .imgt_db_map = data.frame(),
      .gene_groups = list(),
      .package = pkg
    ),
    "Species 'platypus' is not supported for TCR/BCR queries\\."
  )
})

testthat::test_that("loadIMGT() calls getIMGT(refresh = FALSE) and returns sequences", {
  cache <- withr::local_tempdir()
  pkg <- "immReferent"
  
  # Put a cached IGHV file into place
  target_dir <- file.path(cache, "human", "bcr", "ighv")
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  fpath <- file.path(target_dir, "imgt_human_IGHV.fasta")
  .write_dna_fasta(fpath)
  
  # Use a minimal map so loadIMGT -> getIMGT can find the cached file (no download)
  .fake_map <- data.frame(
    gene = c("IGHV"),
    type = c("NUC"),
    cache_subdir = c(file.path("bcr", "ighv")),
    filename_prefix = c("imgt_"),
    query_prefix = c("V-REGION+GENE+and+allele="),
    query_label = c(""),
    stringsAsFactors = FALSE
  )
  
  with_mocked_bindings(
    {
      res <- loadIMGT(species = "human", gene = "IGHV", type = "NUC", suppressMessages = TRUE)
      expect_s4_class(res, "DNAStringSet")
      expect_gte(length(res), 2L)
    },
    .get_cache_dir = function() cache,
    .ensure_cache_dir = function() NULL,
    # If code tried to download, we'd error—ensures wrapper loads cache
    is_imgt_available = function() stop("should not be called"),
    .imgt_db_map = .fake_map,
    .gene_groups = list(),
    .fetch_imgt_query = function(...) stop("should not be called"),
    .update_cache_metadata = function(...) NULL,
    .package = pkg
  )
})

testthat::test_that("refreshIMGT() calls getIMGT(refresh = TRUE)", {
  cache <- withr::local_tempdir()
  pkg <- "immReferent"
  
  # Minimal mapping for IGHV, but force a "download" by mocking refresh path
  .fake_map <- data.frame(
    gene = c("IGHV"),
    type = c("NUC"),
    cache_subdir = c(file.path("bcr", "ighv")),
    filename_prefix = c("imgt_"),
    query_prefix = c("V-REGION+GENE+and+allele="),
    query_label = c(""),
    stringsAsFactors = FALSE
  )
  
  with_mocked_bindings(
    {
      res <- refreshIMGT(species = "human", gene = "IGHV", type = "NUC", suppressMessages = TRUE)
      expect_s4_class(res, "DNAStringSet")
      expect_gte(length(res), 2L)
    },
    .get_cache_dir = function() cache,
    .ensure_cache_dir = function() NULL,
    is_imgt_available = function() TRUE,  # required because refresh forces download
    .imgt_db_map = .fake_map,
    .gene_groups = list(),
    .fetch_imgt_query = function(url, dest_file, species_key) {
      dir.create(dirname(dest_file), recursive = TRUE, showWarnings = FALSE)
      .write_dna_fasta(dest_file)
      invisible(NULL)
    },
    .update_cache_metadata = function(...) NULL,
    .package = pkg
  )
})

testthat::test_that("listIMGT() returns files or empty vector with message", {
  cache <- withr::local_tempdir()
  pkg <- "immReferent"
  
  # 1) No cache dir
  tmp_noexist <- file.path(cache, "nope", "cache")
  expect_snapshot(
    with_mocked_bindings(
      {
        out <- listIMGT()
        expect_identical(out, character(0))
      },
      .get_cache_dir = function() tmp_noexist,
      .package = pkg
    )
  )
  
  # 2) Cache exists with files
  target_dir <- file.path(cache, "human", "bcr", "ighv")
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  f1 <- file.path(target_dir, "imgt_human_IGHV.fasta")
  .write_dna_fasta(f1)
  
  out2 <- with_mocked_bindings(
    {
      listIMGT()
    },
    .get_cache_dir = function() cache,
    .package = pkg
  )
  expect_true(length(out2) >= 1L)
  expect_true(any(grepl("IGHV\\.fasta$", out2)))
})
