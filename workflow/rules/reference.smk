rule get_genome:
    output:
        genome
    params:
        species = ref_species,
        datatype = "dna",
        build = ref_build,
        release = ref_release
    log:
        "workflow/logs/reference/get_genome.log"
    wrapper:
        "v9.12.0/bio/reference/ensembl-sequence"


rule get_transcriptome:
    output:
        transcriptome
    params:
        species = ref_species,
        datatype = "cdna",
        build = ref_build,
        release = ref_release
    log:
        "workflow/logs/reference/get_transcriptome.log"
    wrapper:
        "v9.12.0/bio/reference/ensembl-sequence"


rule get_annotation:
    output:
        annotation
    params:
        species = ref_species,
        build = ref_build,
        release = ref_release,
        flavor = ""
    log:
        "workflow/logs/reference/get_annotation.log"
    wrapper:
        "v9.4.2/bio/reference/ensembl-annotation"

rule make_decoys:
    input:
        genome = genome

    output:
        decoys = decoys

    shell:
        """
        gzip -dc {input.genome} \
        | grep "^>" \
        | cut -d " " -f 1 \
        | sed 's/>//' \
        > {output.decoys}
        """

rule make_gentrome:
    input:
        transcriptome = transcriptome,
        genome = genome

    output:
        gentrome = gentrome

    log:
        "workflow/logs/reference/make_gentrome.log"

    shell:
        "cat {input.transcriptome} {input.genome} > {output.gentrome} 2> {log}"

rule salmon_index:
    input:
        sequences = gentrome,
        decoys = decoys

    output:
        directory(salmon_index)

    conda:
        "../envs/salmon.yml"

    threads:
        8

    params:
        extra = ""

    resources:
        tmpdir = "data/tmp/salmon_index"

    log:
        "workflow/logs/reference/salmon_index.log"

    wrapper:
        "v9.17.0/bio/salmon/index"