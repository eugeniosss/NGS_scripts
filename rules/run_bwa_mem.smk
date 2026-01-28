rule bwa_mem:
    input:
        fastq = lambda wc: get_final_fastqs(wc.SEQ),
        idx=lambda wc: multiext(get_reference_for_mapping(), ".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        bam=maybe_temp(config['mapping']['output_dir']+"/{SEQ}_unsorted.bam")
    conda:
        config["dir"] + "envs/NGS.yml"
    params:
        ref=lambda wc: get_reference_for_mapping(),
        extra=config["mapping"]["params"],
        parser=config["mapping"]["reads_group_parser"]
    log:
        config['mapping']['output_dir'] + "/{SEQ}_mapping.log"
    threads: config["mapping"]["threads"]
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
            bwa mem \
                -t {threads} \
                {params.extra} \
                -R "$RG" \
                {params.ref} \
                {input.fastq} \
            | samtools view -Sb - \
        > {output.bam}) &> {log}
        """