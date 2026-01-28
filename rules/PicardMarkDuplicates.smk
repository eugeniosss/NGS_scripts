rule PicardRemoveDups:
    input:
        config['mapping']['output_dir']+"/{SEQ}.bam"
    output:
        bam=temp_if_not_final(config['dedup']['output_dir']+"/{SEQ}.bam"),
        metrics=config['dedup']['output_dir']+"/{SEQ}.metrics.txt"
    conda:
        config["dir"] + "envs/NGS.yml"
    params:
        extra=config['dedup']['params'],
        mem_to_use=lambda wc: int(config["dedup"]["mem_mb"] / 1024)
    resources:
        mem_mb=int(config["dedup"]["mem_mb"]*config["mem_overhead"]),
    log:
        config['dedup']['output_dir']+"/{SEQ}.log",
    threads: 1
    shell:
        """
        (picard MarkDuplicates \
        -Xmx{params.mem_to_use}G \
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
    resources:
        mem_mb = 500
    conda:
        config["dir"] + "envs/NGS.yml"
    shell:
        """
        (samtools index {input})&> {log}
        """