rule concatenate:
    input:
        fq1=config['adapterremoval2']['output_dir']+"/{SEQ}.pair1.truncated.gz",  
        #fq2=config['adapterremoval2']['output_dir']+"/{SEQ}.pair2.truncated.gz",
        #collapsed=config['adapterremoval2']['output_dir']+"/{SEQ}.collapsed.gz"  
    output:
        config['adapterremoval2']['output_dir']+"/{SEQ}.concatenated.gz"
    params:
        prefix=config['adapterremoval2']['output_dir']+"/{SEQ}"
    shell: "cat {params.prefix}.pair1.truncated.gz {params.prefix}.pair2.truncated.gz {params.prefix}.collapsed.gz > {output}"
