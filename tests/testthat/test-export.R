# tests/testthat/test-export.R

# Helper function to create test DNAStringSet
.create_test_dna_seqs <- function() {
  seqs <- Biostrings::DNAStringSet(c(
    "ATGCGATCGATCGATCGATCGATCGATCG",
    "ATGCGATCGATCGATCGATCG",
    "ATGCGATCGATCGATCG",
    "ATGCGATCGATC",
    "ATGCGATC",
    "ATGC"
  ))
  names(seqs) <- c(
    "IGHV1-2*01",
    "IGHV3-11*02",
    "IGHD1-1*01",
    "IGHJ1*01",
    "IGHJ4*02",
    "IGHC*01"
  )
  seqs
}

# Helper function to create TCR test sequences
.create_test_tcr_seqs <- function() {
  seqs <- Biostrings::DNAStringSet(c(
    "ATGCGATCGATCGATCGATCGATCGATCG",
    "ATGCGATCGATCGATCGATCG",
    "ATGCGATCGATCGATCG",
    "ATGCGATC"
  ))
  names(seqs) <- c(
    "TRBV1-1*01",
    "TRBD1*01",
    "TRBJ1-1*01",
    "TRBC1*01"
  )
  seqs
}

# Helper function to create AAStringSet
.create_test_aa_seqs <- function() {
  seqs <- Biostrings::AAStringSet(c(
    "MSTKVLRQFG",
    "MSTKVLRQ",
    "MSTKV"
  ))
  names(seqs) <- c("IGHV1-2*01", "IGHJ1*01", "IGHC*01")
  seqs
}

# ==============================================================================
# Tests for exportMiXCR()
# ==============================================================================

testthat::test_that("exportMiXCR() creates separate files for V, D, J, C segments", {
  seqs <- .create_test_dna_seqs()
  output_dir <- withr::local_tempdir()

  result <- exportMiXCR(seqs, output_dir, chain = "IGH")

  # Check that all expected files were created
  expect_true(!is.null(result$v_genes))
  expect_true(!is.null(result$d_genes))
  expect_true(!is.null(result$j_genes))
  expect_true(!is.null(result$c_genes))

  # Check files exist

  expect_true(file.exists(result$v_genes))
  expect_true(file.exists(result$d_genes))
  expect_true(file.exists(result$j_genes))
  expect_true(file.exists(result$c_genes))

  # Check file naming convention
  expect_match(basename(result$v_genes), "^v-genes\\.IGH\\.fasta$")
  expect_match(basename(result$d_genes), "^d-genes\\.IGH\\.fasta$")
  expect_match(basename(result$j_genes), "^j-genes\\.IGH\\.fasta$")
  expect_match(basename(result$c_genes), "^c-genes\\.IGH\\.fasta$")
})

testthat::test_that("exportMiXCR() writes correct FASTA content", {
  seqs <- .create_test_dna_seqs()
  output_dir <- withr::local_tempdir()

  result <- exportMiXCR(seqs, output_dir, chain = "IGH")

  # Read V gene file and check content
  v_content <- readLines(result$v_genes)
  expect_true(any(grepl("^>IGHV1-2\\*01$", v_content)))
  expect_true(any(grepl("^>IGHV3-11\\*02$", v_content)))
  expect_true(any(grepl("^ATGCGATCGATCGATCGATCGATCGATCG$", v_content)))

  # Read D gene file
  d_content <- readLines(result$d_genes)
  expect_true(any(grepl("^>IGHD1-1\\*01$", d_content)))

  # Read J gene file
  j_content <- readLines(result$j_genes)
  expect_true(any(grepl("^>IGHJ1\\*01$", j_content)))
  expect_true(any(grepl("^>IGHJ4\\*02$", j_content)))
})

testthat::test_that("exportMiXCR() handles AAStringSet", {
  seqs <- .create_test_aa_seqs()
  output_dir <- withr::local_tempdir()

  result <- exportMiXCR(seqs, output_dir, chain = "IGH")

  expect_true(!is.null(result$v_genes))
  expect_true(!is.null(result$j_genes))
  expect_true(!is.null(result$c_genes))

  # Check amino acid content
  v_content <- readLines(result$v_genes)
  expect_true(any(grepl("^MSTKVLRQFG$", v_content)))
})

testthat::test_that("exportMiXCR() creates output directory if needed", {
  seqs <- .create_test_dna_seqs()
  temp_base <- withr::local_tempdir()
  output_dir <- file.path(temp_base, "new", "nested", "dir")

  expect_false(dir.exists(output_dir))

  result <- exportMiXCR(seqs, output_dir, chain = "IGH")

  expect_true(dir.exists(output_dir))
  expect_true(file.exists(result$v_genes))
})

testthat::test_that("exportMiXCR() errors on sequences without names", {
  seqs <- Biostrings::DNAStringSet(c("ATGC", "GCTA"))

  output_dir <- withr::local_tempdir()

  expect_error(
    exportMiXCR(seqs, output_dir, chain = "IGH"),
    "Sequences must have names"
  )
})

testthat::test_that("exportMiXCR() warns when no matching sequences found", {
  seqs <- .create_test_dna_seqs()
  output_dir <- withr::local_tempdir()

  expect_warning(
    result <- exportMiXCR(seqs, output_dir, chain = "TRB"),
    "No sequences matching chain 'TRB'"
  )

  expect_equal(length(result), 0)
})

testthat::test_that("exportMiXCR() handles all supported chains", {
  chains <- c("IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG")

  for (chain in chains) {
    # Create sequences for this chain
    seqs <- Biostrings::DNAStringSet(c("ATGCGATCGATC", "ATGCGATC"))
    names(seqs) <- c(paste0(chain, "V1*01"), paste0(chain, "J1*01"))

    output_dir <- withr::local_tempdir()
    result <- exportMiXCR(seqs, output_dir, chain = chain)

    expect_true(!is.null(result$v_genes), info = paste("Failed for chain:", chain))
    expect_true(!is.null(result$j_genes), info = paste("Failed for chain:", chain))
  }
})

testthat::test_that("exportMiXCR() removes extra annotation from headers", {
  seqs <- Biostrings::DNAStringSet(c("ATGCGATCGATC"))
  names(seqs) <- c("IGHV1-2*01|Homo sapiens|F|V-REGION|1..296")

  output_dir <- withr::local_tempdir()
  result <- exportMiXCR(seqs, output_dir, chain = "IGH")

  content <- readLines(result$v_genes)
  expect_true(any(grepl("^>IGHV1-2\\*01$", content)))
  expect_false(any(grepl("Homo sapiens", content)))
})

# ==============================================================================
# Tests for exportTRUST4()
# ==============================================================================

testthat::test_that("exportTRUST4() creates FASTA file with correct format", {
  seqs <- .create_test_dna_seqs()
  output_file <- tempfile(fileext = ".fa")
  withr::defer(unlink(output_file))

  result <- exportTRUST4(seqs, output_file)

  expect_equal(result, output_file)
  expect_true(file.exists(output_file))

  content <- readLines(output_file)
  # Check headers are simplified
  expect_true(any(grepl("^>IGHV1-2\\*01$", content)))
  expect_true(any(grepl("^>IGHD1-1\\*01$", content)))
  expect_true(any(grepl("^>IGHJ1\\*01$", content)))
  expect_true(any(grepl("^>IGHC\\*01$", content)))
})

testthat::test_that("exportTRUST4() excludes constant regions when requested", {
  seqs <- .create_test_dna_seqs()
  output_file <- tempfile(fileext = ".fa")
  withr::defer(unlink(output_file))

  result <- exportTRUST4(seqs, output_file, include_constant = FALSE)

  content <- readLines(output_file)
  # V, D, J should be present
  expect_true(any(grepl("^>IGHV", content)))
  expect_true(any(grepl("^>IGHD", content)))
  expect_true(any(grepl("^>IGHJ", content)))
  # C should be absent
  expect_false(any(grepl("^>IGHC", content)))
})

testthat::test_that("exportTRUST4() errors on AAStringSet input", {
  seqs <- .create_test_aa_seqs()
  output_file <- tempfile(fileext = ".fa")

  expect_error(
    exportTRUST4(seqs, output_file),
    "sequences must be a DNAStringSet"
  )
})

testthat::test_that("exportTRUST4() errors on sequences without names", {
  seqs <- Biostrings::DNAStringSet(c("ATGC", "GCTA"))
  output_file <- tempfile(fileext = ".fa")

  expect_error(
    exportTRUST4(seqs, output_file),
    "Sequences must have names"
  )
})

testthat::test_that("exportTRUST4() creates output directory if needed", {
  seqs <- .create_test_dna_seqs()
  temp_base <- withr::local_tempdir()
  output_file <- file.path(temp_base, "nested", "dir", "output.fa")

  expect_false(dir.exists(dirname(output_file)))

  result <- exportTRUST4(seqs, output_file)

  expect_true(dir.exists(dirname(output_file)))
  expect_true(file.exists(output_file))
})

testthat::test_that("exportTRUST4() removes extra annotation from headers", {
  seqs <- Biostrings::DNAStringSet(c("ATGCGATCGATC"))
  names(seqs) <- c("IGHV1-2*01|Homo sapiens|F|V-REGION|1..296")

  output_file <- tempfile(fileext = ".fa")
  withr::defer(unlink(output_file))

  exportTRUST4(seqs, output_file)

  content <- readLines(output_file)
  expect_true(any(grepl("^>IGHV1-2\\*01$", content)))
  expect_false(any(grepl("Homo sapiens", content)))
})

testthat::test_that("exportTRUST4() errors when no valid sequences after filtering", {
  # Create only constant region sequences
  seqs <- Biostrings::DNAStringSet(c("ATGCGATC"))
  names(seqs) <- c("IGHC*01")

  output_file <- tempfile(fileext = ".fa")

  expect_error(
    exportTRUST4(seqs, output_file, include_constant = FALSE),
    "No valid sequences to export"
  )
})

# ==============================================================================
# Tests for exportCellRanger()
# ==============================================================================

testthat::test_that("exportCellRanger() creates FASTA file with correct format", {
  seqs <- .create_test_dna_seqs()
  output_file <- tempfile(fileext = ".fa")
  withr::defer(unlink(output_file))

  result <- exportCellRanger(seqs, output_file)

  expect_equal(result, output_file)
  expect_true(file.exists(output_file))

  content <- readLines(output_file)
  # Check headers are simplified
  expect_true(any(grepl("^>IGHV1-2\\*01$", content)))
})

testthat::test_that("exportCellRanger() errors on AAStringSet input", {
  seqs <- .create_test_aa_seqs()
  output_file <- tempfile(fileext = ".fa")

  expect_error(
    exportCellRanger(seqs, output_file),
    "sequences must be a DNAStringSet"
  )
})

testthat::test_that("exportCellRanger() errors on empty sequences", {
  seqs <- Biostrings::DNAStringSet()
  output_file <- tempfile(fileext = ".fa")

  expect_error(
    exportCellRanger(seqs, output_file),
    "No sequences to export"
  )
})

testthat::test_that("exportCellRanger() errors on sequences without names", {
  seqs <- Biostrings::DNAStringSet(c("ATGC", "GCTA"))
  output_file <- tempfile(fileext = ".fa")

  expect_error(
    exportCellRanger(seqs, output_file),
    "Sequences must have names"
  )
})

testthat::test_that("exportCellRanger() creates output directory if needed", {
  seqs <- .create_test_dna_seqs()
  temp_base <- withr::local_tempdir()
  output_file <- file.path(temp_base, "nested", "dir", "output.fa")

  expect_false(dir.exists(dirname(output_file)))

  result <- exportCellRanger(seqs, output_file)

  expect_true(dir.exists(dirname(output_file)))
  expect_true(file.exists(output_file))
})

testthat::test_that("exportCellRanger() removes extra annotation from headers", {
  seqs <- Biostrings::DNAStringSet(c("ATGCGATCGATC"))
  names(seqs) <- c("IGHV1-2*01|Homo sapiens|F|V-REGION extra info")

  output_file <- tempfile(fileext = ".fa")
  withr::defer(unlink(output_file))

  exportCellRanger(seqs, output_file)

  content <- readLines(output_file)
  expect_true(any(grepl("^>IGHV1-2\\*01$", content)))
  expect_false(any(grepl("Homo sapiens", content)))
  expect_false(any(grepl("extra info", content)))
})

# ==============================================================================
# Tests for exportIgBLAST()
# ==============================================================================

testthat::test_that("exportIgBLAST() creates separate files for V, D, J segments", {
  seqs <- .create_test_dna_seqs()
  output_dir <- withr::local_tempdir()

  result <- exportIgBLAST(seqs, output_dir, organism = "human", receptor_type = "ig")

  # Check that expected files were created
  expect_true(!is.null(result$v_genes))
  expect_true(!is.null(result$d_genes))
  expect_true(!is.null(result$j_genes))

  # Check files exist
  expect_true(file.exists(result$v_genes))
  expect_true(file.exists(result$d_genes))
  expect_true(file.exists(result$j_genes))

  # Check file naming convention
  expect_match(basename(result$v_genes), "^human_ig_v\\.fasta$")
  expect_match(basename(result$d_genes), "^human_ig_d\\.fasta$")
  expect_match(basename(result$j_genes), "^human_ig_j\\.fasta$")
})

testthat::test_that("exportIgBLAST() writes correct FASTA content", {
  seqs <- .create_test_dna_seqs()
  output_dir <- withr::local_tempdir()

  result <- exportIgBLAST(seqs, output_dir, organism = "human", receptor_type = "ig")

  # Read V gene file and check content
  v_content <- readLines(result$v_genes)
  expect_true(any(grepl("^>IGHV1-2\\*01$", v_content)))
  expect_true(any(grepl("^>IGHV3-11\\*02$", v_content)))
  expect_true(any(grepl("^ATGCGATCGATCGATCGATCGATCGATCG$", v_content)))
})

testthat::test_that("exportIgBLAST() handles TCR sequences", {
  seqs <- .create_test_tcr_seqs()
  output_dir <- withr::local_tempdir()

  result <- exportIgBLAST(seqs, output_dir, organism = "human", receptor_type = "tcr")

  expect_true(!is.null(result$v_genes))
  expect_true(!is.null(result$d_genes))
  expect_true(!is.null(result$j_genes))

  # Check file naming for TCR
  expect_match(basename(result$v_genes), "^human_tcr_v\\.fasta$")
  expect_match(basename(result$d_genes), "^human_tcr_d\\.fasta$")
  expect_match(basename(result$j_genes), "^human_tcr_j\\.fasta$")

  # Check content
  v_content <- readLines(result$v_genes)
  expect_true(any(grepl("^>TRBV1-1\\*01$", v_content)))
})

testthat::test_that("exportIgBLAST() creates output directory if needed", {
  seqs <- .create_test_dna_seqs()
  temp_base <- withr::local_tempdir()
  output_dir <- file.path(temp_base, "new", "nested", "dir")

  expect_false(dir.exists(output_dir))

  result <- exportIgBLAST(seqs, output_dir, organism = "human")

  expect_true(dir.exists(output_dir))
  expect_true(file.exists(result$v_genes))
})

testthat::test_that("exportIgBLAST() errors on AAStringSet input", {
  seqs <- .create_test_aa_seqs()
  output_dir <- withr::local_tempdir()

  expect_error(
    exportIgBLAST(seqs, output_dir),
    "sequences must be a DNAStringSet"
  )
})

testthat::test_that("exportIgBLAST() errors on sequences without names", {
  seqs <- Biostrings::DNAStringSet(c("ATGC", "GCTA"))
  output_dir <- withr::local_tempdir()

  expect_error(
    exportIgBLAST(seqs, output_dir),
    "Sequences must have names"
  )
})

testthat::test_that("exportIgBLAST() warns when no matching sequences found", {
  seqs <- .create_test_dna_seqs()  # Contains IG sequences
  output_dir <- withr::local_tempdir()

  expect_warning(
    result <- exportIgBLAST(seqs, output_dir, receptor_type = "tcr"),
    "No sequences matching receptor type 'tcr'"
  )

  expect_equal(length(result), 0)
})

testthat::test_that("exportIgBLAST() removes extra annotation from headers", {
  seqs <- Biostrings::DNAStringSet(c("ATGCGATCGATC", "ATGC"))
  names(seqs) <- c(
    "IGHV1-2*01|Homo sapiens|F|V-REGION|1..296",
    "IGHJ1*01 some extra info"
  )

  output_dir <- withr::local_tempdir()
  result <- exportIgBLAST(seqs, output_dir, organism = "human", receptor_type = "ig")

  v_content <- readLines(result$v_genes)
  expect_true(any(grepl("^>IGHV1-2\\*01$", v_content)))
  expect_false(any(grepl("Homo sapiens", v_content)))

  j_content <- readLines(result$j_genes)
  expect_true(any(grepl("^>IGHJ1\\*01$", j_content)))
  expect_false(any(grepl("extra info", j_content)))
})

testthat::test_that("exportIgBLAST() uses custom organism name", {
  seqs <- .create_test_dna_seqs()
  output_dir <- withr::local_tempdir()

  result <- exportIgBLAST(seqs, output_dir, organism = "custom_species", receptor_type = "ig")

  expect_match(basename(result$v_genes), "^custom_species_ig_v\\.fasta$")
  expect_match(basename(result$d_genes), "^custom_species_ig_d\\.fasta$")
  expect_match(basename(result$j_genes), "^custom_species_ig_j\\.fasta$")
})

testthat::test_that("exportIgBLAST() does not create C gene file (not used by IgBLAST)", {
  seqs <- .create_test_dna_seqs()
  output_dir <- withr::local_tempdir()

  result <- exportIgBLAST(seqs, output_dir, organism = "human", receptor_type = "ig")

  # IgBLAST doesn't use C genes, so no c_genes file should be created
  expect_null(result$c_genes)
})
