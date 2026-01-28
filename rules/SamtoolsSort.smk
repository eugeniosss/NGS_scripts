rule sort_bam:
    wildcard_constraints:
        SEQ="[^_]+.*(?<!_unsorted)"
    input:
        bam=config['mapping']['output_dir']+"/{SEQ}_unsorted.bam"
    output:
        bam=maybe_temp(config['mapping']['output_dir']+"/{SEQ}.bam")
    threads: config["sorting"]["threads"]
    resources:
        mem_mb=int(config["sorting"]["mem_mb"]*config["mem_overhead"])
    params:
        mem_per_thread=lambda wildcards, threads: int(config["sorting"]["mem_mb"] / config["sorting"]["threads"])
    conda:
        config["dir"] + "envs/NGS.yml"
    log:
        config['mapping']['output_dir']+"/{SEQ}_sort.log"
    shell:
        """
        samtools sort \
            -@ {threads} \
            -o {output.bam} \
            -m {params.mem_per_thread}M \
            {input.bam} &> {log}
        """