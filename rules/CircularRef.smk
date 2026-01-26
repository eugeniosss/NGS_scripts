REF_RAW = config["ref"].rsplit(".", 1)[0]
REF_EXT = config["ref"].rsplit(".", 1)[1]

rule circularize_reference:
    input:
        ref = config["ref"]
    output:
        REF_RAW + config["circular_mapper"]["elongation"] + REF_EXT
    params:
        chr = config["circular_mapper"]["chr"],
        ext = config["circular_mapper"]["elongation"]
    conda:
        config["dir"] + "envs/NGS.yml"
    log:
        config["config"]["mapping"]["output_dir"]+"circ.log"
    shell:
        r"""
        (circulargenerator fasta \
            -i {input.ref} \
            -s {params.chr} \
            -e {params.ext}) &> {log}
        """

rule bwa_index_circ:
    input:
        REF_RAW + config["circular_mapper"]["elongation"] + REF_EXT
    output:
        idx=multiext(REF_RAW + config["circular_mapper"]["elongation"] + REF_EXT, ".amb", ".ann", ".bwt", ".pac", ".sa")
    log:
        config["config"]["mapping"]["output_dir"]+"bwa_index_circ.log"
    conda:
        config["dir"]+ "envs/NGS.yml"
    threads: 1
    shell:
        """
        (bwa index {input}) &> {log}
        """

rule circularmapper_realign:
    input:
        #bam = lambda wc: f"{config['dedup']['output_dir']}/{wc.unit}.bam",
        #bai = lambda wc: f"{config['dedup']['output_dir']}/{wc.unit}.bam.bai",
        bam = lambda wc: f"{config['dedup']['output_dir']}/{wc.unit}.bam",
        bai = lambda wc: f"{config['dedup']['output_dir']}/{wc.unit}.bam.bai",
        ref = config["ref"],
    output:
        bam = config['circular_mapper']['output_dir']+"/{SEQ}.bam",
        bai = config['circular_mapper']['output_dir']+"/{SEQ}.bam.bai"
    params:
        ext = config["circular_mapper"]["extension"]
    threads: 1
    conda:
        config["dir"] + "envs/NGS.yml"
    log:
        config["config"]["circular_mapper"]["output_dir"]+"{{unit}}.log"
    shell:
        r"""
        realignsamfile \
            -i {input.bam} \
            -r {input.ref} \
            -e {params.ext}

        mv {input.bam.rsplit(".bam",1)[0]}_realigned.bam {output.bam}

        samtools index {output.bam}
        """