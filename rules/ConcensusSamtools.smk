rule ConsensusSamtools:
    input:
        bam = lambda wc: get_final_bam(wc.unit),
        bai = lambda wc: get_final_bai(wc.unit),
    output:
        fasta = config["consensus"]["output_dir_prefix"] + "{chr}/{unit}.fasta",
        missing = config["consensus"]["output_dir_prefix"] + "{chr}/{unit}.missing.txt"
    params:
        extra = config["consensus"]["params"]
    log:
        config["consensus"]["output_dir_prefix"] + "{chr}/{unit}.log"
    conda:
        config["dir"] + "envs/NGS.yml"
    shell:
        r"""
        # Ensure output directory exists
        mkdir -p $(dirname {output.fasta})

        # 1. Generate consensus FASTA for this chromosome
        samtools consensus \
            -f fasta \
            -a \
            -r {wildcards.chr} \
            {params.extra} \
            -o {output.fasta} \
            {input.bam} 2> {log}

        sed -i "1s/.*/>{wildcards.unit}/" {output.fasta}

        # 2. Calculate missing data (Ns)
        SEQ=$(grep -v "^>" {output.fasta} | tr -d '\n')
        TOTAL=$(echo -n "$SEQ" | wc -c)
        NS=$(echo -n "$SEQ" | tr -cd 'Nn' | wc -c)

        if [ "$TOTAL" -eq 0 ]; then
            PCT="NA"
        else
            PCT=$(echo "scale=4; $NS / $TOTAL" | bc)
        fi

        echo -e "unit\tchr\ttotal_bases\tn_bases\tmissing_pct" > {output.missing}
        echo -e "{wildcards.unit}\t{wildcards.chr}\t$TOTAL\t$NS\t$PCT" >> {output.missing}
        """
