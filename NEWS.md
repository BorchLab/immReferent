# immReferent VERSION 0.99.6

## New Features

* Added export functions for popular immune repertoire analysis tools:
    - `exportMiXCR()`: Export sequences for MiXCR custom library building
    - `exportTRUST4()`: Export sequences for TRUST4 analysis
    - `exportCellRanger()`: Export sequences for 10x Cell Ranger VDJ reference
    - `exportIgBLAST()`: Export sequences for IgBLAST database creation
* Export functions work with sequences from both IMGT and OGRDB sources

## Documentation Improvements

* Standardized roxygen2 documentation across all functions with proper
  formatting using `\code{}`, `\strong{}`, `\itemize{}`, and `\url{}`
* Added comprehensive `@seealso` sections linking related functions
* Improved parameter documentation with explicit type specifications
* Enhanced package-level documentation with function overview
* Fixed typo in `is_ogrdb_available()` documentation (was incorrectly
  referencing IMGT)

# immReferent VERSION 0.99.5

* Added package level man page
* Updated R version dependency from 4.2.0 to 4.5.0
* Added installation instruction via Bioconductor to README and vignette
* Renamed getting-started vignette to immReferent

# immReferent VERSION 0.99.4

* Adding return documentation to OGRDB functions
* Requiring R 4.2 due to use of pipe

# immReferent VERSION 0.99.3

* Adding .Rhistory to gitignore
* Added immReferent to profile on support.bioconductor.org

# immReferent VERSION 0.99.3

* Adding .Rhistory to gitignore
* Added immReferent to profile on support.bioconductor.org

# immReferent VERSION 0.99.2

Add Support for OGRDB

* internal .fetch_airr_files()
* getAIRR main downloading function
* loadAIRR for loading sequences into environment
* refreshAIRR to update sequences
* listAIRR to list archived versions

# immReferent VERSION 0.99.1

* Redirected HLA-based query to IPD-HLA/IMGT FTP site


# immReferent VERSION 0.99.0

Initial release

