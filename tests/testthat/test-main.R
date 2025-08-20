# This file contains tests for the main user-facing functions.

# Helper to get path to test data
get_test_data_path <- function() {
  system.file("extdata", "test_seqs.fasta", package = "immReferent")
}

test_that("getIMGT works for HLA data (offline)", {
  # Create a temporary cache directory for this test
  temp_cache <- withr::local_tempdir()
  withr::local_options(list(immReferent.cache = temp_cache))

  # Mock the download function to just copy our test file
  mock_fetch_hla <- function(file_type, dest_file) {
    file.copy(get_test_data_path(), dest_file)
  }

  with_mocked_bindings(
    .fetch_hla_file = mock_fetch_hla,
    {
      # Call the function, which should now use the mock
      hla_seqs <- getIMGT(gene = "HLA", type = "NUC", suppressMessages = TRUE)

      # Assertions
      expect_s4_class(hla_seqs, "DNAStringSet")
      expect_equal(length(hla_seqs), 2) # Both seqs from the file
      expect_true(any(grepl("HLA:HLA00001", names(hla_seqs))))
    }
  )
})

test_that("getIMGT works for TCR/BCR data (offline)", {
  temp_cache <- withr::local_tempdir()
  withr::local_options(list(immReferent.cache = temp_cache))

  # Mock the IMGT query function
  mock_fetch_query <- function(url, dest_file, species_key) {
    file.copy(get_test_data_path(), dest_file)
  }

  with_mocked_bindings(
    .fetch_imgt_query = mock_fetch_query,
    {
      ighv_seqs <- getIMGT(species = "human", gene = "IGHV", type = "NUC", suppressMessages = TRUE)

      # Assertions
      expect_s4_class(ighv_seqs, "DNAStringSet")
      expect_equal(length(ighv_seqs), 2)
      expect_true(any(grepl("IGHV1-2\\*01", names(ighv_seqs))))
    }
  )
})

test_that("listIMGT correctly lists cached files", {
  temp_cache <- withr::local_tempdir()
  withr::local_options(list(immReferent.cache = temp_cache))

  # Create some dummy files and directories
  dir.create(file.path(temp_cache, "human", "vdj"), recursive = TRUE)
  file.create(file.path(temp_cache, "human", "vdj", "test1.fasta"))
  file.create(file.path(temp_cache, "human", "vdj", "test2.fasta"))

  cached_files <- listIMGT()

  expect_equal(length(cached_files), 2)
  expect_true(any(grepl("test1.fasta", cached_files)))
})

test_that("is_imgt_available works when online", {
  # This test will only run if there is an internet connection
  skip_if_offline(host = "www.imgt.org")

  # We expect it to return a single logical value
  expect_is(is_imgt_available(), "logical")
})

test_that("refreshIMGT calls getIMGT with refresh=TRUE", {
  # We can test this by mocking getIMGT itself
  mock_getIMGT <- mockery::mock(NULL)

  with_mocked_bindings(
    getIMGT = mock_getIMGT,
    {
      refreshIMGT(species = "human", gene = "IGHV", type = "NUC")
    }
  )

  # Check that getIMGT was called with refresh = TRUE
  args <- mockery::mock_args(mock_getIMGT)
  expect_true(args[[1]]$refresh)
})
