rule tximport_gene:
    input:
        quant = expand(os.path.join(quant_path, "{accession}", "quant.sf"), accession=runs),
        tx_to_gene = tx2gene
    output:
        txi = os.path.join(tximport_path, "tximport_gene.rds")
    params:
        extra = "type='salmon', ignoreTxVersion=TRUE, dropInfReps=TRUE"
    resources:
        mem_mb = 16000,
        runtime = 20
    wrapper:
        "v9.3.0/bio/tximport"

rule tximport_gene_lengthscaled:
    input:
        quant = expand(os.path.join(quant_path, "{accession}", "quant.sf"), accession=runs),
        tx_to_gene = tx2gene
    output:
        txi = os.path.join(tximport_path, "tximport_gene_lengthscaled.rds")
    params:
        extra = "type='salmon', ignoreTxVersion=TRUE, dropInfReps=TRUE, countsFromAbundance='lengthScaledTPM'"
    resources:
        mem_mb = 16000,
        runtime = 20
    wrapper:
        "v9.3.0/bio/tximport"

rule tximport_transcript:
    input:
        quant = expand(os.path.join(quant_path, "{accession}", "quant.sf"), accession=runs)
    output:
        txi = os.path.join(tximport_path, "tximport_transcript.rds")
    params:
        extra = "type='salmon', txOut=TRUE, dropInfReps=TRUE"
    resources:
        mem_mb = 16000,
        runtime = 20
    wrapper:
        "v9.3.0/bio/tximport"