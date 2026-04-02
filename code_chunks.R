library(airway)
library(plyxp)
library(purrr)

se <- airway
rowRanges(se) <- NULL

# setup
data <- new_plyxp(se) |>
  mutate(
    rows(gene_length = gene_seq_end - gene_seq_start),
    cols(size_factor = colSums(.assays_asis$counts) / 1e6)
  ) |>
  # keep only the columns we care about for the rest of the code chunks
  select(
    counts,
    rows(gene_biotype, gene_length),
    cols(dex, avgLength, size_factor)
  )
print(data)

data |>
  mutate(
    log10_counts = log10(counts + 1),
    cols(treated = (dex == "trt")),
  )

data |>
  filter(
    rows(gene_biotype == "lincRNA"),
    cols(avgLength >= 126)
  )

data |>
  mutate(
    norm_cts = counts / .cols$size_factor,
    cts_bp = counts / .rows$gene_length
  )

data |>
  mutate(
    rows(
      gene_ave = rowMeans(.assays_asis$counts),
      gene_ave2 = map_dbl(.assays$counts, mean)
    )
  )
