rule coverage_gatk3:
    input:
        bam = lambda wc: get_final_bam(wc.unit),
        bai = lambda wc: get_final_bai(wc.unit),
        # Corrected lambda: takes 'w' (wildcards) and accesses w.bedtype
        bed = lambda w: config["coverages"]["beds"][w.bedtype],
        ref = config["ref"],
        fai = config["ref"] + ".fai",
        dict = lambda wc: config["ref"].rsplit(".", 1)[0] + ".dict"
    output:
        multiext(
            config["coverages"]["output_dir_prefix"] + "{bedtype}/{unit}",
            ".sample_interval_summary",
            ".sample_cumulative_coverage_counts",
            ".sample_cumulative_coverage_proportions",
            ".sample_interval_statistics",
            ".sample_statistics",
            ".sample_summary"
        )
    params:
        gatk_params = config["coverages"].get("params", ""),
        out=config["coverages"]["output_dir_prefix"] + "{bedtype}/{unit}",
        mem_to_use=lambda wc : config["coverages"]["mem_mb"]
    log:
        config["coverages"]["output_dir_prefix"] + "{bedtype}/{unit}_coverage.log"
    threads: 1
    resources:
        mem_mb=int(config["coverages"]["mem_mb"]*config["mem_overhead"])
    conda:
        config["dir"] + "envs/gatk3.yml"
    shell:
        r"""
        gatk3 -T DepthOfCoverage \
            -Xmx{params.mem_to_use}M \
            -R {input.ref} \
            -I {input.bam} \
            -L {input.bed} \
            {params.gatk_params} \
            -o {params.out} 2> {log}
        """
