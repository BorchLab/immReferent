# tests/testthat/test-cache.R

test_that(".get_cache_dir() honors option override and default path format", {
  # Option override
  tmp_cache <- withr::local_tempdir()
  withr::local_options(immReferent.cache = tmp_cache)
  expect_identical(.get_cache_dir(), tmp_cache)
  
  # Default path (option cleared)
  withr::local_options(immReferent.cache = NULL)
  expected <- file.path(path.expand("~"), ".immReferent")
  expect_identical(.get_cache_dir(), expected)
})

# -----------------------------
# .ensure_cache_dir()
# -----------------------------

test_that(".ensure_cache_dir() creates directory when missing and is idempotent", {
  cache <- withr::local_tempdir()
  target <- file.path(cache, "custom-cache-root")
  pkg <- "immReferent"
  
  # 1) Creates when missing + emits message
  expect_message(
    out <- with_mocked_bindings(
      {
        .ensure_cache_dir()
      },
      .get_cache_dir = function() target,
      .package = pkg
    ),
    regexp = "Creating immReferent cache directory"
  )
  expect_true(dir.exists(target))
  expect_identical(out, target)
  
  # 2) Second call: directory exists -> no message
  expect_silent(
    with_mocked_bindings(
      { .ensure_cache_dir() },
      .get_cache_dir = function() target,
      .package = pkg
    )
  )
})

# -----------------------------
# .update_cache_metadata()
# -----------------------------

test_that(".update_cache_metadata() updates an existing YAML file", {
  skip_if_not_installed("yaml")
  cache <- withr::local_tempdir()
  pkg <- "immReferent"
  
  # Pre-seed a minimal YAML so the branch that reads existing file is used
  yaml_file <- file.path(cache, "immReferent_log.yaml")
  yaml::write_yaml(
    list(
      source = "seed",
      package_version = "0.0.0",
      downloaded_items = list()
    ),
    yaml_file
  )
  
  # Run updater (should not need to mock utils::packageVersion)
  with_mocked_bindings(
    {
      .update_cache_metadata(species = "human", gene = "IGHV")
    },
    .get_cache_dir = function() cache,
    .package = pkg
  )
  
  # Validate YAML contents
  meta <- yaml::read_yaml(yaml_file)
  expect_true("last_updated" %in% names(meta))
  expect_true("downloaded_items" %in% names(meta))
  expect_true("human_IGHV" %in% names(meta$downloaded_items))
  expect_match(meta$downloaded_items[["human_IGHV"]], "^\\d{4}-\\d{2}-\\d{2}$")
})

