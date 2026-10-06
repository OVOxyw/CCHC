library(tidyverse)
library(jsonlite)

# Load sample metadata
samples <- read_tsv(
  "config/samples.tsv",
  show_col_types = FALSE
)

# Extract fastp and Salmon QC metrics for all samples
qc_metrics <- map_dfr(samples$run, function(run) {

  fastp_file <- file.path(
    "output", "fastp", paste0(run, ".json")
  )

  salmon_file <- file.path(
    "data", "quant", run, "aux_info", "meta_info.json"
  )

  fastp <- fromJSON(fastp_file)
  salmon <- fromJSON(salmon_file)

  tibble(
    run = run,

    # fastp
    raw_reads_total =
      fastp$summary$before_filtering$total_reads,

    trimmed_reads_total =
      fastp$summary$after_filtering$total_reads,

    retained_pct =
      fastp$summary$after_filtering$total_reads /
      fastp$summary$before_filtering$total_reads * 100,

    q30_before =
      fastp$summary$before_filtering$q30_rate * 100,

    q30_after =
      fastp$summary$after_filtering$q30_rate * 100,

    gc_before =
      fastp$summary$before_filtering$gc_content * 100,

    gc_after =
      fastp$summary$after_filtering$gc_content * 100,

    duplication_pct =
      fastp$duplication$rate * 100,

    # Salmon
    salmon_processed_fragments =
      salmon$num_processed,

    salmon_mapped_fragments =
      salmon$num_mapped,

    mapping_rate =
      salmon$percent_mapped,

    library_type =
      salmon$detected_library_type
  )
})

# Add sample metadata and round QC percentages
qc_summary <- samples %>%
  left_join(qc_metrics, by = "run") %>%
  mutate(
    across(
      c(
        retained_pct,
        q30_before,
        q30_after,
        gc_before,
        gc_after,
        duplication_pct,
        mapping_rate
      ),
      ~ round(.x, 2)
    )
  )

# Save summary table
write_tsv(
  qc_summary,
  "output/qc/sample_qc_summary.tsv"
)

cat(
  "QC summary written for",
  nrow(qc_summary),
  "samples\n"
)