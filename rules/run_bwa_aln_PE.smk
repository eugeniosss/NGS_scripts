rule bwa_aln_pe:
    input:
        r1=lambda wc: get_final_fastqs(wc.SEQ)[0],
        r2=lambda wc: get_final_fastqs(wc.SEQ)[1],
        idx=multiext(config["ref"], ".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        r1_sai=config['mapping']['output_dir']+"/{SEQ}_R1.sai",
        r2_sai=config['mapping']['output_dir']+"/{SEQ}_R2.sai"
    threads: config["mapping"]["threads"]
    conda:
        config["dir"] + "envs/NGS.yml"
    params:
        extra=config["mapping"]["params"],
        ref=config["ref"]
    log:
        config['mapping']['output_dir']+"/{SEQ}_aln.log"
    shell:
        r"""
        (bwa aln {params.extra} -t {threads} {params.ref} {input.r1} > {output.r1_sai} &&
        bwa aln {params.extra} -t {threads} {params.ref} {input.r2} > {output.r2_sai}) &> {log}
        """

rule bwa_sampe:
    input:
        r1=lambda wc: get_final_fastqs(wc.SEQ)[0],
        r2=lambda wc: get_final_fastqs(wc.SEQ)[1],
        r1_sai=config['mapping']['output_dir']+"/{SEQ}_R1.sai",
        r2_sai=config['mapping']['output_dir']+"/{SEQ}_R2.sai",
        idx=multiext(config["ref"], ".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        bam=config['mapping']['output_dir']+"/{SEQ}.bam"
    threads: 1
    conda:
        config["dir"] + "envs/NGS.yml"
    params:
        ref=config["ref"],
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
        | samtools sort -o {output.bam} - ) &> {log}
        """
