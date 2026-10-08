library(tidyverse)
library(jsonlite)
library(ngsReports)

# Inputs from Snakemake
samples_file <- as.character(
  snakemake@input[["samples"]]
)
fastp_files <- unlist(
  snakemake@input[["fastp"]],
  use.names = FALSE
)
salmon_files <- unlist(
  snakemake@input[["salmon"]],
  use.names = FALSE
)
rrna_files <- unlist(
  snakemake@input[["rrna"]],
  use.names = FALSE
)
hbrna_files <- unlist(
  snakemake@input[["hbrna"]],
  use.names = FALSE
)
output_file <- as.character(
  snakemake@output[["summary"]]
)


# Sample metadata
samples <- read_tsv(
  samples_file,
  show_col_types = FALSE
)

# fastp QC
fastp_qc <- map_dfr(
  fastp_files,
  function(file) {
    x <- fromJSON(file)
    tibble(
      run = basename(file) |>
        str_remove("\\.json$"),
      raw_reads_total =
        x$summary$before_filtering$total_reads,
      trimmed_reads_total =
        x$summary$after_filtering$total_reads,
      retained_pct =
        x$summary$after_filtering$total_reads /
        x$summary$before_filtering$total_reads * 100,
      q30_before =
        x$summary$before_filtering$q30_rate * 100,
      q30_after =
        x$summary$after_filtering$q30_rate * 100,
      gc_before =
        x$summary$before_filtering$gc_content * 100,
      gc_after =
        x$summary$after_filtering$gc_content * 100,
      duplication_pct =
        x$duplication$rate * 100
    )
  }
)

# Salmon QC
salmon_qc <- map_dfr(
  salmon_files,
  function(file) {
    run <- basename(dirname(file))
    meta_file <- file.path(
      dirname(file),
      "aux_info",
      "meta_info.json"
    )
    x <- fromJSON(meta_file)
    tibble(
      run = run,
      salmon_processed_fragments =
        x$num_processed,
      salmon_mapped_fragments =
        x$num_mapped,
      mapping_rate =
        x$percent_mapped,
      library_type =
        x$detected_library_type
    )
  }
)

# flagstat QC
summarise_flagstat <- function(files, prefix) {
  x <- importNgsLogs(
    files,
    type = "flagstat"
  ) |>
    filter(
      flag %in% c(
        "primary",
        "primary mapped"
      )
    ) |>
    mutate(
      run = str_remove(
        Filename,
        "_[12]\\.[^.]+\\.flagstat$"
      ),
      read = str_match(
        Filename,
        "_([12])\\."
      )[, 2]
    ) |>
    select(
      run,
      read,
      flag,
      passed = `QC-passed`
    ) |>
    pivot_wider(
      names_from = flag,
      values_from = passed
    ) |>
    rename(
      primary_reads = primary,
      primary_mapped = `primary mapped`
    ) |>
    mutate(
      pct =
        primary_mapped /
        primary_reads * 100
    ) |>
    select(
      run,
      read,
      primary_reads,
      primary_mapped,
      pct
    ) |>
    pivot_wider(
      names_from = read,
      values_from = c(
        primary_reads,
        primary_mapped,
        pct
      ),
      names_glue = "r{read}_{.value}"
    ) |>
    rename_with(
      ~ paste0(prefix, "_", .x),
      -run
    )

  x
}

# rRNA QC
rrna_qc <- summarise_flagstat(
  rrna_files,
  "rrna"
)

# HbRNA QC
hbrna_qc <- summarise_flagstat(
  hbrna_files,
  "hbrna"
)

# Check that every sample has QC
check_runs <- function(x, label) {
  missing_runs <- setdiff(
    samples$run,
    x$run
  )
  if (length(missing_runs) > 0) {
    stop(
      label,
      " QC missing for: ",
      paste(
        missing_runs,
        collapse = ", "
      )
    )
  }
}

check_runs(fastp_qc, "fastp")
check_runs(salmon_qc, "Salmon")
check_runs(rrna_qc, "rRNA")
check_runs(hbrna_qc, "HbRNA")

# Final QC table
qc_summary <- samples |>
  left_join(
    fastp_qc,
    by = "run"
  ) |>
  left_join(
    salmon_qc,
    by = "run"
  ) |>
  left_join(
    rrna_qc,
    by = "run"
  ) |>
  left_join(
    hbrna_qc,
    by = "run"
  ) |>
  mutate(
    across(
      c(
        retained_pct,
        q30_before,
        q30_after,
        gc_before,
        gc_after,
        duplication_pct,
        mapping_rate,
        rrna_r1_pct,
        rrna_r2_pct,
        hbrna_r1_pct,
        hbrna_r2_pct
      ),
      ~ round(.x, 2)
    )
  )

# Write output
dir.create(
  dirname(output_file),
  recursive = TRUE,
  showWarnings = FALSE
)
write_tsv(
  qc_summary,
  output_file
)
cat(
  "QC summary written for",
  nrow(qc_summary),
  "samples\n"
)