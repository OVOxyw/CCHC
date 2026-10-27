library(rtracklayer)
library(dplyr)
library(readr)

annotation_file <- as.character(snakemake@input[["annotation"]])
output_file <- as.character(snakemake@output[["tx2gene"]])
gtf <- rtracklayer::import(annotation_file, feature.type = "transcript")
tx2gene <- as.data.frame(gtf) |>
  filter(!is.na(transcript_id), !is.na(gene_id)) |>
  select(transcript_id, gene_id) |>
  distinct()
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write_tsv(tx2gene, output_file)
cat("tx2gene written:", nrow(tx2gene), "transcript-gene pairs\n")
