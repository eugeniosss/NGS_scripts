rule bwa_aln:
    input:
        fastq=config['adapterremoval2']['output_dir']+"/{SEQ}.concatenated.gz",
        idx=multiext(config["ref"], ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        config['mapping']['output_dir']+"/{SEQ}.sai"
    conda:
        config["dir"] + "envs/NGS.yml"
    params:
        extra=config["mapping"]["params"],
        ref=config["ref"]
    log:
        config['mapping']['output_dir']+"/{SEQ}_mapping.log",
    threads: config["mapping"]["threads"]
    shell:
        """
        (bwa aln \
        {params.extra} \
        -t {threads} \
        {params.ref} \
        {input.fastq} \
        > {output} )&> {log}
        """

rule bwa_samse:
    input:
        fastq=config['adapterremoval2']['output_dir']+"/{SEQ}.concatenated.gz",
        sai=config['mapping']['output_dir']+"/{SEQ}.sai",
        idx=multiext(config["ref"], ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        config['mapping']['output_dir']+"/{SEQ}.bam"
    conda:
        config["dir"] + "envs/NGS.yml"
    params:
        ref=config["ref"]
    log:
        config['mapping']['output_dir']+"/{SEQ}_sampe.log",
    threads: 1
    shell:
        """
        (bwa samse \
        {params.ref} \
        {input.sai} \
        {input.fastq} \
        | samtools view -Sb - \
        | samtools sort -o {output} -)&> {log}
        """
