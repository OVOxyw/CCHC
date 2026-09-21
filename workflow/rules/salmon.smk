rule salmon_quant:
    input:
        r1 = os.path.join(trimmed_path, "{accession}_1.fastq.gz"),
        r2 = os.path.join(trimmed_path, "{accession}_2.fastq.gz"),
        index = salmon_index

    output:
        quant = os.path.join(quant_path, "{accession}", "quant.sf"),
        lib = os.path.join(quant_path, "{accession}", "lib_format_counts.json")

    log:
        "workflow/logs/salmon/{accession}.log"

    threads:
        8

    params:
        libtype = "A",
        extra = ""

    conda:
        "../envs/salmon.yml"

    wrapper:
        "v9.16.0/bio/salmon/quant"
