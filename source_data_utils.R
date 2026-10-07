# Helpers to export the source data of each figure.
#
# Sourced from the SectionN/ Rmds as source("../source_data_utils.R"), so the
# output directory is relative to the section folder (SectionN/tables/).
#
# Each figure gets Source_Data_<figure>.xlsx (one sheet per panel) plus one CSV
# per sheet, named Source_Data_<figure>_<sheet>.csv.

# ---- Display labels ----
tool_labels    <- c(bambu = "Bambu", flair = "FLAIR", isoquant = "IsoQuant",
                    kallisto = "Kallisto", oarfish = "Oarfish")
norm_labels    <- c(CPM = "CPM", TPM = "TPM", cqn = "CQN", EDA = "EDASeq",
                    ratio_counts = "SMR", ratio_correction = "SMR")
dataset_labels <- c(NIH = "Mouse", Neuron = "iNeuron")

relabel_with <- function(x, labels) {
  x <- as.character(x)
  out <- unname(labels[x])
  ifelse(is.na(out), x, out)
}
relabel_tool    <- function(x) relabel_with(x, tool_labels)
relabel_norm    <- function(x) relabel_with(x, norm_labels)
relabel_dataset <- function(x) relabel_with(x, dataset_labels)

# ---- Panels ----
# One sheet per panel from a long table, e.g.
#   split_panels(df, c("Fig 3a" = "Mouse", "Fig 3b" = "LongBench"))
split_panels <- function(df, panels, col = "Dataset") {
  lapply(panels, function(p) df[as.character(df[[col]]) == p, , drop = FALSE])
}

# ---- Writer ----
write_source_data <- function(sheets, figure, dir = "tables") {
  bad <- names(sheets)[nchar(names(sheets)) > 31]
  if (length(bad)) stop("Excel sheet names must be <= 31 characters: ", toString(bad))
  dir.create(dir, showWarnings = FALSE)

  sheets <- lapply(sheets, function(x) as.data.frame(dplyr::ungroup(x)))
  name <- paste0("Source_Data_", figure)
  writexl::write_xlsx(sheets, file.path(dir, paste0(name, ".xlsx")))
  for (sh in names(sheets)) {
    csv_name <- sub("_$", "", paste0(name, "_", gsub("[^A-Za-z0-9]+", "_", sh)))
    readr::write_excel_csv(sheets[[sh]], file.path(dir, paste0(csv_name, ".csv")))
  }
  invisible(sheets)
}
