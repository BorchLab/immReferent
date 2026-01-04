# tests/testthat/test-mainOGRDB.R

# --- Small writers for simulated downloads (FASTA / AIRR JSON) ----------------
.write_dna_fasta <- function(path) {
  writeLines(c(">IGHV1-1*01", "ATGGCTGCT", ">IGHV1-2*01", "ATGGCAGCT"), path)
}
.write_airr_json <- function(path) {
  airr <- list(
    germline_set = list(
      allele_descriptions = list(
        list(label = "IGHV1-1*01", sequence = "ATGGCTGCT"),
        list(label = "IGHV1-2*01", sequence = "ATGGCAGCT")
      )
    )
  )
  jsonlite::write_json(airr, path, auto_unbox = TRUE, pretty = FALSE)
}

# ------------------------------------------------------------------------------
testthat::test_that("getOGRDB() — FASTA gapped NUC: downloads when missing, loads DNAStringSet", {
  cache <- withr::local_tempdir()
  pkg <- "immReferent"
  
  with_mocked_bindings(
    {
      res <- getOGRDB(
        species = "human", locus = "IGH", type = "NUC",
        format = "FASTA_GAPPED", suppressMessages = TRUE
      )
      expect_s4_class(res, "DNAStringSet")
      expect_gte(length(res), 2L)
    },
    # mock internals
    .get_cache_dir = function() cache,
    .ensure_cache_dir = function() NULL,
    is_ogrdb_available = function() TRUE,
    .ogrdb_build_url = function(species_api, set_name, version, format, species_subgroup, human_reference, use_api) {
      # Return a fake URL (unused content-wise)
      paste("mock://", species_api, set_name, version, format, sep = "/")
    },
    .ogrdb_fetch_to_file = function(url, dest_file) {
      dir.create(dirname(dest_file), recursive = TRUE, showWarnings = FALSE)
      .write_dna_fasta(dest_file)
      dest_file
    },
    .update_cache_metadata = function(...) NULL,
    .package = pkg
  )
})

testthat::test_that("getOGRDB() — PROT translation from FASTA works and returns AAStringSet", {
  cache <- withr::local_tempdir()
  pkg <- "immReferent"
  
  with_mocked_bindings(
    {
      res <- getOGRDB(
        species = "human", locus = "IGH",
        type = "PROT", format = "FASTA_UNGAPPED",
        suppressMessages = TRUE
      )
      expect_s4_class(res, "AAStringSet")
      expect_gte(length(res), 2L)
    },
    .get_cache_dir = function() cache,
    .ensure_cache_dir = function() NULL,
    is_ogrdb_available = function() TRUE,
    .ogrdb_build_url = function(...) "mock://airr/fasta/ungapped",
    .ogrdb_fetch_to_file = function(url, dest_file) {
      dir.create(dirname(dest_file), recursive = TRUE, showWarnings = FALSE)
      .write_dna_fasta(dest_file)  # DNA in, code path translates to AA
      dest_file
    },
    .update_cache_metadata = function(...) NULL,
    .package = pkg
  )
})

testthat::test_that("getOGRDB() — AIRR JSON parsed to DNAStringSet/AAStringSet", {
  cache <- withr::local_tempdir()
  pkg <- "immReferent"
  
  # NUC from AIRR JSON
  with_mocked_bindings(
    {
      res <- getOGRDB(
        species = "human", locus = "IGK",
        type = "NUC", format = "AIRR",
        suppressMessages = TRUE
      )
      expect_s4_class(res, "DNAStringSet")
      expect_gte(length(res), 2L)
    },
    .get_cache_dir = function() cache,
    .ensure_cache_dir = function() NULL,
    is_ogrdb_available = function() TRUE,
    .ogrdb_build_url = function(...) "mock://airr/json",
    .ogrdb_fetch_to_file = function(url, dest_file) {
      dir.create(dirname(dest_file), recursive = TRUE, showWarnings = FALSE)
      .write_airr_json(dest_file)
      dest_file
    },
    .update_cache_metadata = function(...) NULL,
    .package = pkg
  )
  
  # PROT from AIRR JSON (translation path)
  with_mocked_bindings(
    {
      res <- getOGRDB(
        species = "human", locus = "IGL",
        type = "PROT", format = "AIRR",
        suppressMessages = TRUE
      )
      expect_s4_class(res, "AAStringSet")
      expect_gte(length(res), 2L)
    },
    .get_cache_dir = function() cache,
    .ensure_cache_dir = function() NULL,
    is_ogrdb_available = function() TRUE,
    .ogrdb_build_url = function(...) "mock://airr/json",
    .ogrdb_fetch_to_file = function(url, dest_file) {
      dir.create(dirname(dest_file), recursive = TRUE, showWarnings = FALSE)
      .write_airr_json(dest_file)
      dest_file
    },
    .update_cache_metadata = function(...) NULL,
    .package = pkg
  )
})

testthat::test_that("getOGRDB() — uses cache when present (no download attempted)", {
  cache <- withr::local_tempdir()
  pkg <- "immReferent"
  
  # Place a cached gapped FASTA file
  cache_dir <- file.path(cache, "Human", "ogrdb")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cached <- file.path(cache_dir, "Human_IGH_VDJ_published_gapped.fasta")
  .write_dna_fasta(cached)
  
  with_mocked_bindings(
    {
      expect_silent({
        res <- getOGRDB(
          species = "human", locus = "IGH",
          type = "NUC", format = "FASTA_GAPPED",
          suppressMessages = TRUE
        )
      })
      expect_s4_class(res, "DNAStringSet")
      expect_gte(length(res), 2L)
    },
    .get_cache_dir = function() cache,
    .ensure_cache_dir = function() NULL,
    # If code tries to download, error out to catch it
    is_ogrdb_available = function() stop("should not be called"),
    .ogrdb_build_url = function(...) stop("should not be called"),
    .ogrdb_fetch_to_file = function(...) stop("should not be called"),
    .update_cache_metadata = function(...) stop("should not be called"),
    .package = pkg
  )
})

testthat::test_that("getOGRDB() — species_subgroup is encoded correctly in URL", {
  # This calls the real internal URL builder; adjust if not exported
  skip_if_not(exists(".ogrdb_build_url"), "URL builder not found")
  
  url <- .ogrdb_build_url(
    species_api = "Mouse",
    set_name = "IGKV",
    version = "published",
    format = "gapped",
    species_subgroup = "C57BL/6J",   # has a '/'
    human_reference = FALSE,
    use_api = TRUE
  )
  # OGRDB requires '/' encoded as '%252f' in subgroup
  expect_match(url, "C57BL%252f6J")
  expect_match(url, "/gapped$")  # no _ex for mouse
})

testthat::test_that("getOGRDB() — errors on unsupported species", {
  cache <- withr::local_tempdir()
  pkg <- "immReferent"
  
  expect_error(
    with_mocked_bindings(
      {
        getOGRDB(
          species = "platypus", locus = "IGH",
          type = "NUC", format = "FASTA_GAPPED",
          suppressMessages = TRUE
        )
      },
      .get_cache_dir = function() cache,
      .ensure_cache_dir = function() NULL,
      # minimal map to only include human/mouse
      .ogrdb_species_map = list(
        "human" = list(api = "Human"),
        "mouse" = list(api = "Mouse")
      ),
      .package = pkg
    )
  )
})

testthat::test_that("loadOGRDB() calls getOGRDB(refresh = FALSE) and returns sequences from cache", {
  cache <- withr::local_tempdir()
  pkg <- "immReferent"
  
  # Put a cached file in place
  target_dir <- file.path(cache, "Human", "ogrdb")
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  fpath <- file.path(target_dir, "Human_IGH_VDJ_published_gapped.fasta")
  .write_dna_fasta(fpath)
  
  with_mocked_bindings(
    {
      res <- loadOGRDB(
        species = "human", locus = "IGH",
        type = "NUC", format = "FASTA_GAPPED",
        suppressMessages = TRUE
      )
      expect_s4_class(res, "DNAStringSet")
      expect_gte(length(res), 2L)
    },
    .get_cache_dir = function() cache,
    .ensure_cache_dir = function() NULL,
    is_ogrdb_available = function() stop("should not be called"),
    .ogrdb_build_url = function(...) stop("should not be called"),
    .ogrdb_fetch_to_file = function(...) stop("should not be called"),
    .update_cache_metadata = function(...) stop("should not be called"),
    .package = pkg
  )
})

testthat::test_that("refreshOGRDB() forces re-download", {
  cache <- withr::local_tempdir()
  pkg <- "immReferent"
  
  with_mocked_bindings(
    {
      res <- refreshOGRDB(
        species = "human", locus = "IGK",
        type = "NUC", format = "FASTA_GAPPED",
        suppressMessages = TRUE
      )
      expect_s4_class(res, "DNAStringSet")
      expect_gte(length(res), 2L)
    },
    .get_cache_dir = function() cache,
    .ensure_cache_dir = function() NULL,
    is_ogrdb_available = function() TRUE,
    .ogrdb_build_url = function(...) "mock://force/refresh",
    .ogrdb_fetch_to_file = function(url, dest_file) {
      dir.create(dirname(dest_file), recursive = TRUE, showWarnings = FALSE)
      .write_dna_fasta(dest_file)
      dest_file
    },
    .update_cache_metadata = function(...) NULL,
    .package = pkg
  )
})

testthat::test_that("is_ogrdb_available() returns a single logical (no strict assumption on network)", {
  # Don't make a hard TRUE/FALSE claim (CRAN/network variability); just shape.
  res <- try(is_ogrdb_available(), silent = TRUE)
  # If the function errored (rare), skip
  if (inherits(res, "try-error")) testthat::skip("HEAD failed in this environment")
  expect_type(res, "logical")
  expect_length(res, 1L)
})
