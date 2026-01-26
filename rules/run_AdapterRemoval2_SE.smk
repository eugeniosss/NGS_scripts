rule adapterremoval_se:
    input:
        r1=lambda wc: SEQS[wc.SEQ].split()[0],
        r2=lambda wc: SEQS[wc.SEQ].split()[1]
    output:
        collapsed=config['adapterremoval2']['output_dir'] + "/{SEQ}.collapsed.gz",
        pair1=config['adapterremoval2']['output_dir'] + "/{SEQ}.pair1.truncated.gz",
        pair2=config['adapterremoval2']['output_dir'] + "/{SEQ}.pair2.truncated.gz",
        settings=config['adapterremoval2']['output_dir'] + "/{SEQ}.settings"
    params:
        outbase=lambda wc: f"{config['adapterremoval2']['output_dir']}/{wc.SEQ}",
        options=config["adapterremoval2"]["options"]
    threads: config["adapterremoval2"]["threads"]
    log:
        config['adapterremoval2']['output_dir'] + "/{SEQ}.log"
    conda:
        config["dir"] + "envs/NGS.yml"
    shell:
        """
        AdapterRemoval \
            --file1 {input.r1} \
            --file2 {input.r2} \
            --basename {params.outbase} \
            --threads {threads} \
            --gzip \
            {params.options} \
            &> {log}
        """