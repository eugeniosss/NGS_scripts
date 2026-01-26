rule bwa_aln_SE:
    input:
        fastq = lambda wc: get_final_fastqs(wc.SEQ),
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

rule bwa_samse_SE:
    input:
        fastq = lambda wc: get_final_fastqs(wc.SEQ),
        sai=config['mapping']['output_dir']+"/{SEQ}.sai",
        idx=multiext(config["ref"], ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        config['mapping']['output_dir']+"/{SEQ}.bam",
    conda:
        config["dir"] + "envs/NGS.yml",
    params:
        ref=config["ref"],
        parser=config["mapping"]["reads_group_parser"]
    log:
        config['mapping']['output_dir']+"/{SEQ}_sampe.log",
    threads: 1
    shell:
        r"""
        # Build read group
        if [ "{params.parser}" != "False" ]; then
            parts=($(echo {wildcards.SEQ} | tr '{params.parser}' ' '))
            SAMPLE=${{parts[0]}}
            LA=${{parts[1]:-NA}}
            DATE=${{parts[2]:-NA}}
            RG="@RG\\tID:${{SAMPLE}}_${{LA}}_${{DATE}}\\tSM:${{SAMPLE}}\\tLB:${{SAMPLE}}_${{LA}}"
        else
            RG="@RG\\tID:{wildcards.SEQ}\\tSM:{wildcards.SEQ}\\tLB:{wildcards.SEQ}"
        fi

        (
            bwa samse \
                {params.ref} \
                {input.sai} \
                {input.fastq} \
                -r "$RG" \
            | samtools view -Sb - \
            | samtools sort -o {output} -
        ) &> {log}
        """