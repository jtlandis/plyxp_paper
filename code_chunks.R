library(airway)
library(plyxp)
library(purrr)

data("airway")

se <- airway
rowRanges(se) <- NULL

# setup
se <- new_plyxp(se) |>
  mutate(
    rows(gene_length = gene_seq_end - gene_seq_start),
    cols(mil_of_reads = colSums(.assays_asis$counts) / 1e6)
  ) |>
  # keep only the columns we care about for the rest of the code chunks
  select(
    counts,
    rows(gene_biotype, gene_length),
    cols(dex, mil_of_reads)
  ) |>
  se()

# Code Chunk Listing 1
xp <- new_plyxp(se)
xp |>
  mutate(
    log_counts = log10(counts + 1),
    cols(treated = dex == "trt"),
    cpbp = counts / .rows$gene_length,
    cpm = counts / .cols$mil_of_reads
  ) |>
  filter(rows(gene_biotype == "lincRNA"))

# Code Chunk Listing 2
xp |>
  mutate(rows(
    ave = rowMeans(.assays_asis$counts),
    ave2 = purrr::map_dbl(.assays$counts, mean)
  ))
