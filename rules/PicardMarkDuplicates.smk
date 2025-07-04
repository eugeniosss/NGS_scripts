rule PicardRemoveDups:
    input:
        config['mapping']['output_dir']+"/{SEQ}.bam"
    output:
        bam=config['dedup']['output_dir']+"/{SEQ}.bam",
        metrics=config['dedup']['output_dir']+"/{SEQ}.metrics.txt"
    conda:
        config["dir"] + "envs/NGS.yml"
    params:
        extra=config['dedup']['params']
    resources:
        mem_mb=11000,
    log:
        config['dedup']['output_dir']+"/{SEQ}.log",
    threads: 1
    shell:
        """
        (picard MarkDuplicates \
        {params.extra} \
        --INPUT {input} \
        --OUTPUT {output.bam} \
        --METRICS_FILE {output.metrics})&> {log}
        """

rule index_bams:
    input:
        config['dedup']['output_dir']+"/{SEQ}.bam",
    output:
        config['dedup']['output_dir']+"/{SEQ}.bam.bai",
    log:
        config['dedup']['output_dir']+"/{SEQ}_indexing.log",
    threads: 1
    conda:
        config["dir"] + "envs/NGS.yml"
    shell:
        """
        (samtools index {input})&> {log}
        """