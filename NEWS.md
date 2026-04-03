# immReferent 1.0.2

## New Features

* Added OGRDB (AIRR Community) support with full suite of functions:
    - `getOGRDB()`: Download germline sets from OGRDB
    - `loadOGRDB()`: Load OGRDB sequences into environment
    - `refreshOGRDB()`: Update locally cached OGRDB sequences
    - `listOGRDB()`: List archived OGRDB germline set versions
    - `is_ogrdb_available()`: Check OGRDB service availability
* Added export functions for popular immune repertoire analysis tools:
    - `exportMiXCR()`: Export sequences for MiXCR custom library building
    - `exportTRUST4()`: Export sequences for TRUST4 analysis
    - `exportCellRanger()`: Export sequences for 10x Cell Ranger VDJ reference
    - `exportIgBLAST()`: Export sequences for IgBLAST database creation
* Export functions work with sequences from both IMGT and OGRDB sources
* Redirected HLA-based queries to IPD-IMGT/HLA FTP site for improved reliability
* Added support for additional species: cynomolgus monkey, pig, and dog

## Documentation Improvements

* Added package-level man page (`immReferent-package`)
* Renamed getting-started vignette to `immReferent` and expanded content
* Expanded caching vignette with additional guidance
* Standardized roxygen2 documentation across all functions with proper
  formatting using `\code{}`, `\strong{}`, `\itemize{}`, and `\url{}`
* Added comprehensive `@seealso` sections linking related functions
* Improved parameter documentation with explicit type specifications
* Added return value documentation to OGRDB functions
* Added Bioconductor installation instructions to README and vignette

## Infrastructure

* Set up GitHub Actions CI (R-CMD-check on macOS, Windows, Ubuntu)
* Added test coverage reporting via codecov
* Added PR command workflows for `/document` and `/style`
* CI now runs on `main`, `master`, and `devel` branches
* Updated R version dependency from 4.2.0 to 4.5.0
* Added spelling wordlist

## Testing

* Reorganized test suite: split monolithic test file into focused modules
    - `test-mainIMGT.R`, `test-mainOGRDB.R`, `test-export.R`
    - `test-cache.R`, `test-dowload.R`, `test-utils.R`
* Added comprehensive unit tests for export, OGRDB, cache, download, and
  utility functions

## Bug Fixes

* Fixed typo in `is_ogrdb_available()` documentation (was incorrectly
  referencing IMGT)
* Fixed BiocCheck findings from initial submission

# immReferent 0.99.0

* Initial release with IMGT sequence downloading, caching, and loading
