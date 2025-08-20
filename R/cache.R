#' @title Get the Root Cache Directory
#' @description Retrieves the path to the cache directory. The default location is
#' `~/.immReferent`, but it can be overridden by setting the
#' `immReferent.cache` option (e.g., `options(immReferent.cache = "your/path")`).
#' @return The path to the cache directory.
#' @noRd
.get_cache_dir <- function() {
  getOption("immReferent.cache", default = file.path(path.expand("~"), ".immReferent"))
}

#' @title Ensure the Cache Directory Exists
#' @description Checks if the cache directory exists at the location specified by
#' `.get_cache_dir()`. If it does not exist, it creates the directory.
#' @return Invisibly returns the path to the cache directory.
#' @noRd
.ensure_cache_dir <- function() {
  cache_dir <- .get_cache_dir()
  if (!dir.exists(cache_dir)) {
    message("Creating immReferent cache directory at: ", cache_dir)
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(cache_dir)
}

#' @title Update Cache Metadata YAML File
#' @description Internal function to write or update a YAML file in the cache root
#' that tracks download dates and package version.
#' @param species The species that was downloaded.
#' @param gene The gene/locus that was downloaded.
#' @importFrom yaml read_yaml write_yaml
#' @noRd
.update_cache_metadata <- function(species, gene) {
  cache_dir <- .get_cache_dir()
  yaml_file <- file.path(cache_dir, "immReferent_log.yaml")

  metadata <- if (file.exists(yaml_file)) {
    yaml::read_yaml(yaml_file)
  } else {
    list(
      source = "IMGT/GENE-DB (https://www.imgt.org) and IPD-IMGT/HLA (https://www.ebi.ac.uk/ipd/imgt/hla/)",
      package_version = as.character(utils::packageVersion("immReferent")),
      downloaded_items = list()
    )
  }

  metadata$last_updated <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

  item_key <- paste(species, gene, sep = "_")
  metadata$downloaded_items[[item_key]] <- as.character(Sys.Date())

  yaml::write_yaml(metadata, yaml_file)
}
