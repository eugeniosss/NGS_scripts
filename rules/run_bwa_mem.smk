rule bwa_mem:
    input:
        fastq = lambda wc: get_final_fastqs(wc.SEQ),
        idx=lambda wc: multiext(get_reference_for_mapping(), ".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        bam=config['mapping']['output_dir'] + "/{SEQ}.bam"
    conda:
        config["dir"] + "envs/NGS.yml"
    params:
        ref=lambda wc: get_reference_for_mapping(),
        extra=config["mapping"]["params"],
        parser=config["mapping"]["reads_group_parser"]
    log:
        config['mapping']['output_dir'] + "/{SEQ}_mapping.log"
    threads: config["mapping"]["threads"]
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
            | samtools sort -@{threads} -o {output.bam} -
        ) &> {log}
        """