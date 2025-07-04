rule adapterremoval:
    input:
        r1= lambda wildcards: SEQS[wildcards.SEQ].split()[0],
        r2= lambda wildcards: SEQS[wildcards.SEQ].split()[1]
    output:
        collapsed="{output_dir}/{SEQ}.pair1.truncated.gz"
    params:
        outbase=lambda wildcards: f"{config['adapterremoval2']['output_dir']}/{wildcards.SEQ}",
        options=config["adapterremoval2"]["options"]
    threads: config["adapterremoval2"]["threads"]
    log:
        "{output_dir}/{SEQ}.log"
    conda:
        config["dir"] + "envs/NGS.yml"
    shell:
        """
        AdapterRemoval \
            --file1 {input.r1} \
            --file2 {input.r2} \
            --basename {params.outbase} \
            --threads {threads} \
            {params.options} \
            &> {log}
        """
