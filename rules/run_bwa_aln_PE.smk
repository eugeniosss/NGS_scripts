rule bwa_aln_pe_r1:
    input:
        r1=lambda wc: get_final_fastqs(wc.SEQ)[0],
        idx=lambda wc: multiext(get_reference_for_mapping(), ".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        r1_sai=maybe_temp(config['mapping']['output_dir']+"/{SEQ}_R1.sai"),
    threads: config["mapping"]["threads"]
    resources:
        mem_mb=config["mapping"]["mapping_mem_mb"]
    conda:
        config["dir"] + "envs/NGS.yml"
    params:
        extra=config["mapping"]["params"],
        ref=lambda wc: get_reference_for_mapping()
    log:
        config['mapping']['output_dir']+"/{SEQ}_alnr1.log"
    shell:
        r"""
        (bwa aln {params.extra} -t {threads} {params.ref} {input.r1} > {output.r1_sai}) &> {log}
        """

rule bwa_aln_pe_r2:
    input:
        r2=lambda wc: get_final_fastqs(wc.SEQ)[1],
        idx=lambda wc: multiext(get_reference_for_mapping(), ".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        r2_sai=maybe_temp(config['mapping']['output_dir']+"/{SEQ}_R2.sai")
    threads: config["mapping"]["threads"]
    resources:
        mem_mb=config["mapping"]["mapping_mem_mb"]
    conda:
        config["dir"] + "envs/NGS.yml"
    params:
        extra=config["mapping"]["params"],
        ref=lambda wc: get_reference_for_mapping()
    log:
        config['mapping']['output_dir']+"/{SEQ}_alnr2.log"
    shell:
        r"""
        (bwa aln {params.extra} -t {threads} {params.ref} {input.r2} > {output.r2_sai}) &> {log}
        """

rule bwa_sampe:
    input:
        r1=lambda wc: get_final_fastqs(wc.SEQ)[0],
        r2=lambda wc: get_final_fastqs(wc.SEQ)[1],
        r1_sai=config['mapping']['output_dir']+"/{SEQ}_R1.sai",
        r2_sai=config['mapping']['output_dir']+"/{SEQ}_R2.sai",
        idx=lambda wc: multiext(get_reference_for_mapping(), ".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        bam=maybe_temp(config['mapping']['output_dir']+"/{SEQ}_unsorted.bam")
    threads: 1
    resources:
        mem_mb=config["mapping"]["mapping_mem_mb"]
    conda:
        config["dir"] + "envs/NGS.yml"
    params:
        ref=lambda wc: get_reference_for_mapping(),
        parser=config["mapping"]["reads_group_parser"]
    log:
        config['mapping']['output_dir']+"/{SEQ}_sampe.log"
    shell:
        r"""
        # Build read group
        (if [ "{params.parser}" != "False" ]; then
            parts=($(echo {wildcards.SEQ} | tr '{params.parser}' ' '))
            SAMPLE=${{parts[0]}}
            LA=${{parts[1]:-NA}}
            DATE=${{parts[2]:-NA}}
            RG="@RG\\tID:${{SAMPLE}}_${{LA}}_${{DATE}}\\tSM:${{SAMPLE}}\\tLB:${{SAMPLE}}_${{LA}}"
        else
            RG="@RG\\tID:{wildcards.SEQ}\\tSM:{wildcards.SEQ}\\tLB:{wildcards.SEQ}"
        fi

        bwa sampe \
            -r "$RG" \
            {params.ref} \
            {input.r1_sai} {input.r2_sai} \
            {input.r1} {input.r2} \
        | samtools view -Sb - \
        > {output.bam}) &> {log}
        """
