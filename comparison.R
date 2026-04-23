library(airway)
library(tidySummarizedExperiment)
library(plyxp)
data(airway)
se <- airway
rowRanges(se) <- NULL

native_version <- function(se) {
  rowData(se)$length <- rowData(se)$gene_seq_end - rowData(se)$gene_seq_start
  assays(se)$counts_per_bp <- assay(se, "counts") / rowData(se)$length
  se
}

forced_conversion <- function(se) {
  tib <- tidySummarizedExperiment:::as_tibble.SummarizedExperiment(se)
  res <- tib |>
    dplyr::mutate(
      length = gene_seq_end - gene_seq_start,
      counts_per_bp = counts / length
    )
  tidySummarizedExperiment:::update_SE_from_tibble(res, se)
}

system.time({
  se_update <- forced_conversion(se)
})

xp <- new_plyxp(se)

plyxp_version <- function(xp) {
  xp |>
    mutate(
      rows(length = gene_seq_end - gene_seq_start),
      counts_per_bp = counts / .rows$length
    )
}

system.time({
  xp_update <- plyxp_version(xp)
})

bench::mark(
  # bioconductor native calls
  native_version(se),
  # naive eager tibble conversion
  forced_conversion(se),
  # tidySummarizedExperiment implementation
  mutate(se,
    length = gene_seq_end - gene_seq_start,
    counts_per_bp = counts / length
  ),
  # plyxp implementation
  plyxp_version(xp),
  check = FALSE,
  # small bug related to S7::S7_dispatch calls and how the
  # Rprofiler parsing works within `bench`.
  memory = FALSE,
  # protect against cold starts
  min_iterations = 10
)

all.equal(
  assay(se_update, "counts_per_bp"),
  assay(xp_update, "counts_per_bp")
)
