rule fastp:
    input:
      sample = [
        os.path.join(raw_path, "{accession}_1.fastq.gz"),
        os.path.join(raw_path, "{accession}_2.fastq.gz")
      ]

    output:
        trimmed = [
            os.path.join(trimmed_path, "{accession}_1.fastq.gz"),
            os.path.join(trimmed_path, "{accession}_2.fastq.gz")
        ],
        html = os.path.join(fastp_report_path, "{accession}.html"),
        json = os.path.join(fastp_report_path, "{accession}.json")
    
    log:
        "workflow/logs/fastp/{accession}.log"

    conda:
        "../envs/fastp.yml"

    params:
        extra = "" ,
        adapters = ""

    threads:
        4

    script:
        "../scripts/fastp.py"