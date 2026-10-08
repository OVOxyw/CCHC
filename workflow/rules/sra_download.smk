rule sra_download:
    output:
        r1 = protected(os.path.join(raw_path, "{accession}_1.fastq.gz")),
        r2 = protected(os.path.join(raw_path, "{accession}_2.fastq.gz"))
    log:
        "workflow/logs/sra_download/{accession}.log"
    conda:
        "../envs/sra.yml"
    params:
        extra = config["params"]["sra_download"]["extra"]  
    threads:
        2
    localrule: True
    resources:
        tmpdir = lambda wildcards: f"data/tmp/sra_download/{wildcards.accession}"
    script:
        "../scripts/fasterq-dump.py"