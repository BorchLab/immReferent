# immReferent

An R Package for Accessing and Managing IMGT Sequences

<img align="right" src="https://github.com/ncborcherding/immReferent/blob/main/www/immreferent_hex.png" width="305" height="352">

## Introduction 

immReferent is an R package designed to provide a stable, reproducible, and lightweight interface to IMGT immune receptor (TCR/BCR) and HLA sequences.
It serves as the backbone for computational immunology workflows by ensuring a consistent source of high-quality nucleotide and protein sequences.

### This package enables:

* Downloading IMGT data (receptor and HLA sequences).
* Caching for offline reproducibility.
* Querying metadata (allele, gene, species).
* Interoperability with Bioconductor objects (DNAStringSet, AAStringSet).

### Why immReferent?

As immune repertoire analysis expands, a centralized sequence reference layer is critical for:

* Integrating across tools like [scRepertoire](https://github.com/BorchLab/scRepertoire) and [immApex](https://github.com/BorchLab/immApex)
* Reducing redundancy by avoiding repeated hard-coded IMGT downloads.
* Guaranteeing reproducibility with cached versions.

# Installation

```
# Install development version from GitHub
devtools::install_github("BorchLab/immReferent")
```

# Getting Started

```
library(immReferent)

# Check if IMGT is online
is_imgt_available()

# Download human IGHV sequences
ighv <- getIMGT("IGHV")

# Load cached data
cached <- loadIMGT("IGHV")

# List available datasets
listIMGT()

# Refresh cache
refreshIMGT("IGHV")
```

# IMGT usage

IMGT is used as a reference for gene names and sequence information can be accessed via getIMGT(). Data from IMGT is under a CC BY-NC-ND 4.0 license. Please be aware that attribution is required for usage and it is the intent of IMGT to not allow derivative or commercial usage.

# Bug Reports/New Features

#### If you run into any issues or bugs please submit a [GitHub issue](https://github.com/BorchLab/immReferent/issues) with details of the issue.

If possible please include a [reproducible example](https://reprex.tidyverse.org/). 

#### Any requests for new features or enhancements can also be submitted as [GitHub issues](https://github.com/BorchLab/immReferent/issues).

#### [Pull Requests](https://github.com/BorchLab/immReferent/pulls) are welcome for bug fixes, new features, or enhancements.

## Citation
In Progress