# This file contains the mapping tables to translate user requests
# into IMGT query parameters.

# A data frame that maps gene, type, and data source to query parameters
#' @importFrom tibble tribble
.imgt_db_map <- tibble::tribble(
  ~gene,   ~type,  ~query_prefix, ~query_label,     ~cache_subdir,     ~filename_prefix,
  # VDJ Nucleotides
  "IGHV",  "NUC",  "7.1",         "",               "vdj",             "imgt_",
  "IGHD",  "NUC",  "7.1",         "",               "vdj",             "imgt_",
  "IGHJ",  "NUC",  "7.1",         "",               "vdj",             "imgt_",
  "IGKV",  "NUC",  "7.1",         "",               "vdj",             "imgt_",
  "IGKJ",  "NUC",  "7.1",         "",               "vdj",             "imgt_",
  "IGLV",  "NUC",  "7.1",         "",               "vdj",             "imgt_",
  "IGLJ",  "NUC",  "7.1",         "",               "vdj",             "imgt_",
  "TRAV",  "NUC",  "7.1",         "",               "vdj",             "imgt_",
  "TRAJ",  "NUC",  "7.1",         "",               "vdj",             "imgt_",
  "TRBV",  "NUC",  "7.1",         "",               "vdj",             "imgt_",
  "TRBD",  "NUC",  "7.1",         "",               "vdj",             "imgt_",
  "TRBJ",  "NUC",  "7.1",         "",               "vdj",             "imgt_",
  "TRDV",  "NUC",  "7.1",         "",               "vdj",             "imgt_",
  "TRDD",  "NUC",  "7.1",         "",               "vdj",             "imgt_",
  "TRDJ",  "NUC",  "7.1",         "",               "vdj",             "imgt_",
  "TRGV",  "NUC",  "7.1",         "",               "vdj",             "imgt_",
  "TRGJ",  "NUC",  "7.1",         "",               "vdj",             "imgt_",

  # V-region Amino Acids
  "IGHV",  "PROT", "7.3",         "",               "vdj_aa",          "imgt_aa_",
  "IGKV",  "PROT", "7.3",         "",               "vdj_aa",          "imgt_aa_",
  "IGLV",  "PROT", "7.3",         "",               "vdj_aa",          "imgt_aa_",
  "TRAV",  "PROT", "7.3",         "",               "vdj_aa",          "imgt_aa_",
  "TRBV",  "PROT", "7.3",         "",               "vdj_aa",          "imgt_aa_",
  "TRDV",  "PROT", "7.3",         "",               "vdj_aa",          "imgt_aa_",
  "TRGV",  "PROT", "7.3",         "",               "vdj_aa",          "imgt_aa_",

  # Spliced Leader + V-exon
  "IGHV",  "L-V-EXON", "8.1",    "L-PART1+V-EXON",  "leader_vexon",    "imgt_lv_",
  "IGKV",  "L-V-EXON", "8.1",    "L-PART1+V-EXON",  "leader_vexon",    "imgt_lv_",
  "IGLV",  "L-V-EXON", "8.1",    "L-PART1+V-EXON",  "leader_vexon",    "imgt_lv_",
  "TRAV",  "L-V-EXON", "8.1",    "L-PART1+V-EXON",  "leader_vexon",    "imgt_lv_",
  "TRBV",  "L-V-EXON", "8.1",    "L-PART1+V-EXON",  "leader_vexon",    "imgt_lv_",
  "TRDV",  "L-V-EXON", "8.1",    "L-PART1+V-EXON",  "leader_vexon",    "imgt_lv_",
  "TRGV",  "L-V-EXON", "8.1",    "L-PART1+V-EXON",  "leader_vexon",    "imgt_lv_",

  # Spliced Leader
  "IGHV", "LEADER", "8.1", "L-PART1+L-PART2", "leader", "imgt_IGHL", # Note: original was IGH, but needs V
  "IGKV", "LEADER", "8.1", "L-PART1+L-PART2", "leader", "imgt_IGKL",
  "IGLV", "LEADER", "8.1", "L-PART1+L-PART2", "leader", "imgt_IGLL",
  "TRAV", "LEADER", "8.1", "L-PART1+L-PART2", "leader", "imgt_TRAL",
  "TRBV", "LEADER", "8.1", "L-PART1+L-PART2", "leader", "imgt_TRBL",
  "TRDV", "LEADER", "8.1", "L-PART1+L-PART2", "leader", "imgt_TRDL",
  "TRGV", "LEADER", "8.1", "L-PART1+L-PART2", "leader", "imgt_TRGL",

  # Spliced Constant
  "IGHC",  "NUC-C",  "14.1",        "",               "constant",        "imgt_",
  "IGKC",  "NUC-C",  "7.5",         "",               "constant",        "imgt_",
  "IGLC",  "NUC-C",  "7.5",         "",               "constant",        "imgt_",
  "TRAC",  "NUC-C",  "14.1",        "",               "constant",        "imgt_",
  "TRBC",  "NUC-C",  "14.1",        "",               "constant",        "imgt_",
  "TRGC",  "NUC-C",  "14.1",        "",               "constant",        "imgt_",
  "TRDC",  "NUC-C",  "14.1",        "",               "constant",        "imgt_"
)

# A list to define gene groups for convenience
.gene_groups <- list(
    "IG" = c("IGHV", "IGHD", "IGHJ", "IGKV", "IGKJ", "IGLV", "IGLJ", "IGHC", "IGKC", "IGLC"),
    "TCR" = c("TRAV", "TRAJ", "TRBV", "TRBD", "TRBJ", "TRDV", "TRDD", "TRDJ", "TRGV", "TRGJ", "TRAC", "TRBC", "TRDC", "TRGC"),
    "IGH" = c("IGHV", "IGHD", "IGHJ", "IGHC"),
    "IGK" = c("IGKV", "IGKJ", "IGKC"),
    "IGL" = c("IGLV", "IGLJ", "IGLC"),
    "TRA" = c("TRAV", "TRAJ", "TRAC"),
    "TRB" = c("TRBV", "TRBD", "TRBJ", "TRBC"),
    "TRD" = c("TRDV", "TRDD", "TRDJ", "TRDC"),
    "TRG" = c("TRGV", "TRGJ", "TRGC")
)
