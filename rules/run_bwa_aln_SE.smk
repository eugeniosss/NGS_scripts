rule bwa_aln_SE:
    input:
        fastq = lambda wc: get_final_fastqs(wc.SEQ),
        idx=lambda wc: multiext(get_reference_for_mapping(), ".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        maybe_temp(config['mapping']['output_dir']+"/{SEQ}.sai")
    conda:
        config["dir"] + "envs/NGS.yml"
    params:
        extra=config["mapping"]["params"],
        ref=lambda wc: get_reference_for_mapping(),
    log:
        config['mapping']['output_dir']+"/{SEQ}_mapping.log",
    threads: config["mapping"]["threads"]
    resources:
        mem_mb=config["mapping"]["mapping_mem_mb"]
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
        idx=lambda wc: multiext(get_reference_for_mapping(), ".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        bam=maybe_temp(config['mapping']['output_dir']+"/{SEQ}_unsorted.bam")
    conda:
        config["dir"] + "envs/NGS.yml",
    params:
        ref=lambda wc: get_reference_for_mapping(),
        parser=config["mapping"]["reads_group_parser"]
    log:
        config['mapping']['output_dir']+"/{SEQ}_sampe.log",
    threads: 1
    resources:
        mem_mb=config["mapping"]["mapping_mem_mb"]
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
        > {output.bam}) &> {log}
        """