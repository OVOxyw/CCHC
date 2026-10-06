rule multiqc:
    input:
        fastp = expand(
            os.path.join(
                fastp_report_path,
                "{accession}.json"
            ),
            accession = runs
        ),
        salmon_meta = expand(
            os.path.join(
                quant_path,
                "{accession}",
                "aux_info",
                "meta_info.json"
            ),
            accession = runs
        ),
        salmon_lib = expand(
            os.path.join(
                quant_path,
                "{accession}",
                "lib_format_counts.json"
            ),
            accession = runs
        )

    output:
        os.path.join(qc_path, "multiqc.html")

    params:
        extra = config["params"]["multiqc"]["extra"],
        use_input_files_only = config["params"]["multiqc"]["use_input_files_only"]

    log:
        "workflow/logs/multiqc/multiqc.log"

    threads:
        1

    wrapper:
        "v9.15.0/bio/multiqc"