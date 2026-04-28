library(airway) # the dataset
library(tidySummarizedExperiment)
library(plyxp)
library(purrr)
library(bench) # for benchmarking
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
  tidySE = mutate(
    se,
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

## Comparison of efficient asis computation vs plyxp reshaping

bench::mark(
  asis = mutate(xp, rows(ave = rowMeans(.assays_asis$counts))),
  reshape = mutate(xp, rows(ave = map_dbl(.assays$counts, mean))),
  min_iterations = 10,
  check = FALSE
)

all.equal(
  mutate(xp, rows(ave = rowMeans(.assays_asis$counts))) |>
    pull(rows(ave)) |>
    unname(),
  mutate(xp, rows(ave = map_dbl(.assays$counts, mean))) |>
    pull(rows(ave))
)


## Comparison of Complex Grouop Scenarios

# With increasing complexity of the expression,
# especially in grouping scenarios, plyxp begins to overtake
# native Bioconductor code.
#
# NOTE: We do not include comparisons to tidySummarizedExperiment
#  as summarize will return a tibble and not a SummarizedExperiment
#
# Consider the below where we aggregate counts by gene_biotype.

native_group_summarize <- function() {
  biotypes <- unique(rowData(se)$gene_biotype)
  se_list <- lapply(biotypes, function(bt) {
    idx <- rowData(se)$gene_biotype == bt
    mat <- matrix(
      colSums(assay(se)[idx, , drop = FALSE]),
      nrow = 1,
      dimnames = list(bt, colnames(se))
    )
    SummarizedExperiment(
      assays = list(counts = mat),
      colData = colData(se)
    )
  })
  do.call(rbind, se_list)
}

bench::mark(
  native_group_summarize(),
  plyxp_summarize = xp |>
    group_by(rows(gene_biotype)) |>
    summarize(
      rows(.features = unique(gene_biotype)),
      counts = colSums(counts)
    ),
  check = FALSE,
  memory = FALSE,
  min_iterations = 5
)

# when grouping both rows and columns, plyxp's performance is
# slightly hindered, but native approaches median runtime nearly
# doubles. Below, we use aggregate from the S4Vectors package,
# and precompute RLEs to attempt to save time.

native_group_summarize_both <- function() {
  ro <- order(rowData(se)$gene_biotype)
  co <- order(colData(se)$cell)
  .se <- se[ro, co]

  row_rle <- Rle(rowData(.se)$gene_biotype)
  col_rle <- Rle(colData(.se)$cell)
  # col_ind <- lapply(col_rle, \(rle) start(rle):end(rle))
  out <- aggregate(
    .se,
    by = row_rle,
    FUN = \(x) {
      SummarizedExperiment(
        assays = list(
          counts = matrix(
            tapply(t(assay(x, "counts")), INDEX = col_rle, sum),
            nrow = 1
          ),
          n = matrix(rep(prod(dim(x)), length(runValue(col_rle))), nrow = 1)
        ),
        rowData = DataFrame(
          gene_biotype = rowData(x)$gene_biotype[1],
          n = nrow(x)
        ),
        colData = DataFrame(
          cell = runValue(col_rle),
          n = runLength(col_rle)
        )
      )
    }
  ) |>
    do.call(what = rbind)
  out[order(-rowSums(assay(out, "counts"))), ]
}

# Below is an more optimized version in which indices are cached
# and the entire se object is NOT subsetted per row group.

native_group_summarize_both_optim <- function() {
  ## below code is commented out, techinically reordering the entire
  ## object takes more time too. Instead just get indices with
  ## vctrs::vec_group_loc

  # ro <- order(rowData(se)$gene_biotype)
  # co <- order(colData(se)$cell)
  # .se <- se[ro, co]

  # row_rle <- Rle(rowData(.se)$gene_biotype)
  # col_rle <- Rle(colData(.se)$cell)
  # index <- \(obj) Map(\(s, e) s:e, s = start(obj), e = end(obj))
  # row_rle_ind <- index(row_rle)
  # col_rle_ind <- index(col_rle)
  row_rle_ind <- vctrs::vec_group_loc(rowData(se)$gene_biotype)$loc
  col_rle_ind <- vctrs::vec_group_loc(colData(se)$cell)$loc
  out <- lapply(
    row_rle_ind,
    \(row_ind) {
      row_subset <- assay(.se, "counts")[row_ind, , drop = FALSE]
      mat <- vapply(
        col_rle_ind,
        \(col_ind) row_subset[, col_ind, drop = FALSE] |> sum(),
        numeric(1)
      ) |>
        matrix(nrow = 1)
      n_ <- vapply(col_rle_ind, \(ind) length(ind), numeric(1))
      SummarizedExperiment(
        assays = list(
          counts = mat,
          n = matrix(n_ * length(row_ind), nrow = 1)
        ),
        rowData = DataFrame(
          gene_biotype = rowData(.se)$gene_biotype[row_ind][1],
          n = length(row_ind)
        ),
        colData = DataFrame(
          cell = runValue(col_rle),
          n = runLength(col_rle)
        )
      )
    }
  ) |>
    do.call(what = rbind)
  out[order(-rowSums(assay(out, "counts"))), ]
}

bench::mark(
  native_group_summarize_both(),
  native_group_summarize_both_optim(),
  plyxp_summarize = xp |>
    group_by(rows(gene_biotype), cols(cell)) |>
    summarize(rows(n = n()), n = n(), cols(n = n()), counts = sum(counts)) |>
    arrange(rows(-rowSums(.assays_asis$counts))),
  check = FALSE,
  memory = FALSE,
  min_iterations = 5
)
