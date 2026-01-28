rule ConsensusHtsbox:
    input:
        bam = lambda wc: get_final_bam(wc.unit),
        bai = lambda wc: get_final_bai(wc.unit),
        ref = config["ref"]
    output:
        fasta = config["consensus"]["output_dir_prefix"] + "{chr}/{unit}.fasta",
        missing = config["consensus"]["output_dir_prefix"] + "{chr}/{unit}.missing.txt"
    params:
        chr = "{chr}",
        extra = config["consensus"]["params"]
    resources:
        mem_mb=config["consensus"]["mem_mb"]
    log:
        config["consensus"]["output_dir_prefix"] + "{chr}/{unit}.log"
    conda:
        config["dir"] + "envs/NGS.yml"
    shell:
        r"""
        mkdir -p {config[consensus][output_dir_prefix]}{params.chr}

        # 1. Generate consensus FASTA for this chromosome

        htsbox pileup \
            -f {input.ref} \
            -r {params.chr} \
            {input.bam} \
            {params.extra} \
            > {output.fasta} 2> {log}

        sed -i "1s/.*/>{wildcards.unit}/" {output.fasta}

        # 2. Calculate missing data (Ns)

        TOTAL=$(grep -v "^>" {output.fasta} | tr -d '\n' | wc -c)
        NS=$(grep -v "^>" {output.fasta} | tr -d '\n' | tr -cd 'Nn' | wc -c)

        if [ $TOTAL -eq 0 ]; then
            PCT="NA"
        else
            PCT=$(echo "scale=4; ($NS / $TOTAL) " | bc)
        fi

        echo -e "unit\tchr\ttotal_bases\tn_bases\tmissing_pct" > {output.missing}
        echo -e "{wildcards.unit}\t{wildcards.chr}\t$TOTAL\t$NS\t$PCT" >> {output.missing}
        """