library(airway)
library(tidySummarizedExperiment)
library(plyxp)
data(airway)
se <- airway
rowRanges(se) <- NULL

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
  forced_conversion(se),
  plyxp_version(xp),
  check = FALSE,
  memory = FALSE
)

all.equal(
  assay(se_update, "counts_per_bp"),
  assay(xp_update, "counts_per_bp")
)
