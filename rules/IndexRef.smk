REF_RAW = config["ref"].rsplit(".", 1)[0]
REF_FOR_MAPPING = get_reference_for_mapping()

rule bwa_index:
    input:
        ref = REF_FOR_MAPPING
    output:
        idx = multiext(REF_FOR_MAPPING, ".amb", ".ann", ".bwt", ".pac", ".sa")
    log:
        "bwa_index/bwa_index.log"
    conda:
        config["dir"] + "envs/NGS.yml"
    resources:
        mem_mb = config["mapping"]["indexing_mem_mb"]
    shell:
        """
        (bwa index {input}) &> {log}
        """

rule samtools_faidx:
    input:
        config["ref"]
    output:
        config["ref"] + ".fai"
    log:
        "bwa_index/faidx.log"
    conda:
        config["dir"] + "envs/NGS.yml"
    resources:
        mem_mb = 500
    shell:
        """
        (samtools faidx {input}) &> {log}
        """

rule create_dict:
    input:
        config["ref"]
    output:
        REF_RAW + ".dict"
    log:
        "bwa_index/create_dict.log"
    conda:
        config["dir"] + "envs/NGS.yml"
    resources:
        mem_mb = 500
    shell:
        """
        (picard CreateSequenceDictionary \
            R={input} \
            O={output}) &> {log}
        """
