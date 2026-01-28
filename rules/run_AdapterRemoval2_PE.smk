rule adapterremoval_pe:
    input:
        r1=lambda wc: SEQS[wc.SEQ].split()[0],
        r2=lambda wc: SEQS[wc.SEQ].split()[1]
    output:
        pair1=maybe_temp(config['adapterremoval2']['output_dir'] + "/{SEQ}.pair1.truncated.gz"),
        pair2=maybe_temp(config['adapterremoval2']['output_dir'] + "/{SEQ}.pair2.truncated.gz"),
        settings=config['adapterremoval2']['output_dir'] + "/{SEQ}.settings",
        singleton = temp("adapterremoval2/{SEQ}.singleton.truncated.gz"),
        discarded = temp("adapterremoval2/{SEQ}.discarded.gz"),
        collapsed_trunc = temp("adapterremoval2/{SEQ}.collapsed.truncated.gz")
    params:
        outbase=lambda wc: f"{config['adapterremoval2']['output_dir']}/{wc.SEQ}",
        options=config["adapterremoval2"]["options"]
    threads: config["adapterremoval2"]["threads"]
    resources:
        mem_mb=config["adapterremoval2"]["mem_mb"]
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

        touch {output.singleton}
        touch {output.discarded}
        touch {output.collapsed_trunc}
        """