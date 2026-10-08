rule multiqc:
    input:
        fastp = expand(os.path.join(fastp_report_path, "{accession}.json"), accession=runs),
        salmon = expand(os.path.join(quant_path, "{accession}", "quant.sf"), accession=runs)
    output:
        os.path.join(qc_path, "multiqc.html")
    params:
        extra = config["params"]["multiqc"]["extra"],
        use_input_files_only = config["params"]["multiqc"]["use_input_files_only"]
    log:
        "workflow/logs/multiqc/multiqc.log"
    wrapper:
        "v9.15.0/bio/multiqc"