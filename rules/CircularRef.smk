REF_RAW = config["ref"].rsplit(".", 1)[0]
REF_EXT = config["ref"].rsplit(".", 1)[1]
REF_FOR_MAPPING = get_reference_for_mapping()

rule circularize_reference:
    input:
        ref = config["ref"]
    output:
        REF_FOR_MAPPING
    params:
        chr = config["circular_mapper"]["chr"],
        ext = config["circular_mapper"]["elongation"]
    conda:
        config["dir"] + "envs/NGS.yml"
    log:
        "bwa_index/circ.log"
    resources:
        mem_mb = 500
    shell:
        r"""
        (circulargenerator \
            -i {input.ref} \
            -s {params.chr} \
            -e {params.ext}) &> {log}
        """

rule circularmapper_realign_seqs:
    input:
        bam = config['dedup']['output_dir']+"/{SEQ}.bam",
        ref = config["ref"]
    output:
        bam = config["circular_mapper"]["output_dir"] + "/{SEQ}.bam",
        bai = config["circular_mapper"]["output_dir"] + "/{SEQ}.bam.bai"
    params:
        ext = config["circular_mapper"]["elongation"]
    threads: 1
    resources:
        mem_mb=config["circular_mapper"]["mem_mb"]
    conda:
        config["dir"] + "envs/NGS.yml"
    log:
        config["circular_mapper"]["output_dir"] + "/{SEQ}.log"
    shell:
        r"""
        (realignsamfile \
            -i {input.bam} \
            -r {input.ref} \
            -e {params.ext}

        INP="{input.bam}"
        mv "${{INP%.bam}}_realigned.bam" {output.bam}

        samtools index {output.bam}) &> {log}
        """

rule circularmapper_realign_merged:
    input:
        bam = lambda wc: f"{config['merge_same_sample_runs']['output_dedup_dir']}/{wc.unit}.bam",
        ref = config["ref"],
    output:
        bam = config["merge_same_sample_runs"]["output_circular_dir"] + "/{unit}.bam",
        bai = config["merge_same_sample_runs"]["output_circular_dir"] + "/{unit}.bam.bai",
    params:
        ext = config["circular_mapper"]["elongation"]
    threads: 1
    resources:
        mem_mb=config["circular_mapper"]["mem_mb"]
    conda:
        config["dir"] + "envs/NGS.yml"
    log:
        config["merge_same_sample_runs"]["output_circular_dir"] + "/{unit}.log"
    shell:
        r"""
        (realignsamfile \
            -i {input.bam} \
            -r {input.ref} \
            -e {params.ext}

        INP="{input.bam}"
        mv "${{INP%.bam}}_realigned.bam" {output.bam}

        samtools index {output.bam}) &> {log}
        """