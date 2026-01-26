REF_RAW = config["ref"].rsplit(".", 1)[0]


rule bwa_index:
    input:
        config["ref"]
    output:
        idx=multiext(config["ref"], ".amb", ".ann", ".bwt", ".pac", ".sa")
    log:
        "bwa_index/bwa_index.log"
    conda:
        NGS_SCRIPTS + "envs/NGS.yml"
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
        NGS_SCRIPTS + "envs/NGS.yml"
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
        NGS_SCRIPTS + "envs/gatk3.yml"
    shell:
        """
        (picard CreateSequenceDictionary \
            R={input} \
            O={output}) &> {log}
        """
