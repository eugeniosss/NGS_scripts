rule concatenate:
    input:
        fq1=config['adapterremoval2']['output_dir']+"/{SEQ}.pair1.truncated.gz",  
        fq2=config['adapterremoval2']['output_dir']+"/{SEQ}.pair2.truncated.gz",
        collapsed=config['adapterremoval2']['output_dir']+"/{SEQ}.collapsed.gz"  
    output:
        maybe_temp(config['adapterremoval2']['output_dir']+"/{SEQ}.concatenated.gz")
    resources:
        mem_mb = 500
    shell: "cat {input.fq1} {input.fq2} {input.collapsed} > {output}"


##rule concatenate:
##    input:
##        # Only include files that actually exist
##        #dummy=config['adapterremoval2']['output_dir']+"/{SEQ}.settings",
##        files=lambda wc: [f for f in [
##            f"{config['adapterremoval2']['output_dir']}/{wc.SEQ}.pair1.truncated.gz",
##            f"{config['adapterremoval2']['output_dir']}/{wc.SEQ}.pair2.truncated.gz",
##            f"{config['adapterremoval2']['output_dir']}/{wc.SEQ}.collapsed.gz"
##        ] if os.path.exists(f)]
##    output:
##        f"{config['adapterremoval2']['output_dir']}/{{SEQ}}.concatenated.gz"
##    shell:
##        "cat {input.files} > {output}"