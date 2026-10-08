rule qc_summary:
    input:
        samples = config["samples"]["file"],
        fastp = expand(os.path.join(fastp_report_path, "{accession}.json"), accession=runs),
        salmon = expand(os.path.join(quant_path, "{accession}", "quant.sf"), accession=runs),
        rrna = expand(os.path.join(contamination_qc_path, "{accession}_{read}.rrna.flagstat"), accession=runs, read=["1", "2"]),
        hbrna = expand(os.path.join(contamination_qc_path, "{accession}_{read}.hbrna.flagstat"), accession=runs, read=["1", "2"])
    output:
        summary = os.path.join(qc_summary_path, "sample_qc_summary.tsv"),
        whole_blood = os.path.join(qc_summary_path, "whole_blood_qc_mqc.tsv")
    conda:
        "../envs/qc_summary.yml"
    log:
        "workflow/logs/qc_summary/qc_summary.log"
    threads:
        1
    script:
        "../scripts/build_qc_summary.R"