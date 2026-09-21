rule sra_download:
    output:
        r1 = os.path.join(raw_path, "{accession}_1.fastq.gz"),
        r2 = os.path.join(raw_path, "{accession}_2.fastq.gz")

    log:
        "workflow/logs/sra_download/{accession}.log"

    conda:
        "../envs/sra.yml"

    params:
        extra = "--split-files"

    threads:
        4

    resources:
        tmpdir = lambda wildcards: f"data/tmp/sra_download/{wildcards.accession}"

    script:
        "../scripts/fasterq-dump.py"