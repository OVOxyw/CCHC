rule bwa_index_rrna:
    input:
        rrna_ref
    output:
        rrna_ref + ".amb",
        rrna_ref + ".ann",
        rrna_ref + ".bwt",
        rrna_ref + ".pac",
        rrna_ref + ".sa"
    conda:
        "../envs/bwa_samtools.yml"
    threads:
        1
    log:
        "workflow/logs/contamination_qc/bwa_index_rrna.log"
    shell:
        "bwa index {input} > {log} 2>&1"

rule bwa_index_hbrna:
    input:
        hbrna_ref
    output:
        hbrna_ref + ".amb",
        hbrna_ref + ".ann",
        hbrna_ref + ".bwt",
        hbrna_ref + ".pac",
        hbrna_ref + ".sa"
    conda:
        "../envs/bwa_samtools.yml"
    threads:
        1
    log:
        "workflow/logs/contamination_qc/bwa_index_hbrna.log"
    shell:
        "bwa index {input} > {log} 2>&1"

rule sample_first_million:
    input:
        os.path.join(raw_path, "{accession}_{read}.fastq.gz")

    output:
        temp(
            os.path.join(
                contamination_tmp_path,
                "{accession}_{read}.firstmillion.fastq"
            )
        )
    wildcard_constraints:
        read = "1|2"
    threads:
        1
    log:
        "workflow/logs/contamination_qc/{accession}_{read}_sample.log"
    shell:
        r"""
        set +o pipefail
        gzip -dc {input} | head -n 4000000 > {output}
        set -o pipefail
        test "$(wc -l < {output})" -eq 4000000
        """
    
rule rrna_flagstat:
    input:
        reads = os.path.join(
            contamination_tmp_path,
            "{accession}_{read}.firstmillion.fastq"
        ),
        ref = rrna_ref,
        index = [
            rrna_ref + ".amb",
            rrna_ref + ".ann",
            rrna_ref + ".bwt",
            rrna_ref + ".pac",
            rrna_ref + ".sa"
        ]
    output:
        os.path.join(
            contamination_qc_path,
            "{accession}_{read}.rrna.flagstat"
        )
    conda:
        "../envs/bwa_samtools.yml"
    threads:
        1
    log:
        "workflow/logs/contamination_qc/{accession}_{read}_rrna.log"
    shell:
        """
        bwa mem -t {threads} {input.ref} {input.reads} 2>> {log} \
        | samtools view -hu - 2>> {log} \
        | samtools flagstat - > {output} 2>> {log}
        """
        
rule hbrna_flagstat:
    input:
        reads = os.path.join(
            contamination_tmp_path,
            "{accession}_{read}.firstmillion.fastq"
        ),
        ref = hbrna_ref,
        index = [
            hbrna_ref + ".amb",
            hbrna_ref + ".ann",
            hbrna_ref + ".bwt",
            hbrna_ref + ".pac",
            hbrna_ref + ".sa"
        ]
    output:
        os.path.join(
            contamination_qc_path,
            "{accession}_{read}.hbrna.flagstat"
        )
    conda:
        "../envs/bwa_samtools.yml"
    threads:
        1
    log:
        "workflow/logs/contamination_qc/{accession}_{read}_hbrna.log"
    shell:
        """
        bwa mem -t {threads} {input.ref} {input.reads} 2>> {log} \
        | samtools view -hu - 2>> {log} \
        | samtools flagstat - > {output} 2>> {log}
        """