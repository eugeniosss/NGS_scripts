rule picard_merge_sample_bams:
    input:
        lambda wc: DIC_MERGE[wc.sample]
    output:
        bam = config['merge_same_sample_runs']['output_merge_dir'] + "/{sample}.bam"
    log:
        config['merge_same_sample_runs']['output_merge_dir'] + "/{sample}.log"
    params:
        extra = config["merge_same_sample_runs"]["params"],
        input_list = lambda wildcards, input: " ".join(f"--INPUT {f}" for f in input)
    conda:
        config["dir"] + "envs/NGS.yml"
    threads: 1
    shell:
        r"""
        picard MergeSamFiles \
            {params.extra} \
            {params.input_list} \
            --OUTPUT {output.bam} \
            --TMP_DIR tmp_{wildcards.sample} \
            &> {log}
        """

rule markduplicates_merged:
    input:
        bam = config['merge_same_sample_runs']['output_merge_dir'] + "/{sample}.bam"
    output:
        bam = config['merge_same_sample_runs']['output_dedup_dir'] + "/{sample}.bam",
        metrics = config['merge_same_sample_runs']['output_dedup_dir'] + "/{sample}.metrics.txt"
    conda:
        config["dir"] + "envs/NGS.yml"
    params:
        extra=config['dedup']['params']
    resources:
        mem_mb=11000,
    threads: 1
    log:
        config['merge_same_sample_runs']['output_dedup_dir'] + "/{sample}.markdup.log"
    shell:
        r"""
        (picard MarkDuplicates \
        {params.extra} \
        --INPUT {input} \
        --OUTPUT {output.bam} \
        --METRICS_FILE {output.metrics})&> {log}
        """

rule index_bams_samples:
    input:
        config['merge_same_sample_runs']['output_dedup_dir'] + "/{sample}.bam",
    output:
        config['merge_same_sample_runs']['output_dedup_dir'] + "/{sample}.bam.bai",
    log:
        config['merge_same_sample_runs']['output_dedup_dir']+"/{sample}_indexing.log",
    threads: 1
    conda:
        config["dir"] + "envs/NGS.yml"
    shell:
        r"""
        (samtools index {input}) &> {log}
        """


rule summarize_samples_metrics:
    input:
        basic_stats = expand(
            config["stats"]["output_dir"] + "/{sample}_basic_stats.txt",
            sample=MERGE_SAMPLES
        ),
        coverage = lambda wildcards: (
            expand(
                config["coverage"]["output_dir"] + "/{sample}.seq_summary",
                sample=MERGE_SAMPLES
            ) if config.get("coverage", {}).get("run", False) else []
        ),
        dups = expand(
            config["merge_same_sample_runs"]["output_dedup_dir"] + "/{sample}.metrics.txt",
            sample=MERGE_SAMPLES
        ),
        consensus = lambda wildcards: (
            expand(
                config["consensus"]["output_dir_prefix"] + "{chr}/{sample}.missing.txt",
                sample=MERGE_SAMPLES,
                chr=config["consensus"]["chrs"]
            ) if config.get("consensus", {}).get("run", False) else []
        )
    output:
        summary = "samples_summary.tsv"
    threads: 1
    run:
        import pandas as pd
        import os

        OUTFILE = output.summary
        all_data = []

        for sample in MERGE_SAMPLES:
            row = {"Sample": sample}

            # --- Basic stats ---
            basic_file = os.path.join(config["stats"]["output_dir"], f"{sample}_basic_stats.txt")
            if os.path.exists(basic_file):
                with open(basic_file) as f:
                    lines = [line.strip().split() for line in f.readlines()]
                    if len(lines) >= 2:
                        for k, v in zip(lines[0], lines[1]):
                            if not k.startswith("Raw_Reads"):  # skip R1/R2
                                row[k] = v

            # --- Duplicates ---
            dups_file = os.path.join(config["dedup"]["output_dir"], f"{sample}.metrics.txt")
            if os.path.exists(dups_file):
                with open(dups_file) as f:
                    lines = [l.strip() for l in f if l.strip() != "" and not l.startswith("#")]
                    data_line = None
                    for i, line in enumerate(lines):
                        if line.startswith("LIBRARY") and i+1 < len(lines):
                            data_line = lines[i+1]
                            break
                    if data_line:
                        values = data_line.split("\t")
                        row["Dups_amount"] = values[5]
                        row["Dups_rate"] = values[8]
                    else:
                        row["Dups_amount"] = "NA"
                        row["Dups_rate"] = "NA"
            else:
                row["Dups_amount"] = "NA"
                row["Dups_rate"] = "NA"

            # --- Coverage metrics ---
            if config.get("coverage", {}).get("run", False):
                for bed_name, bed_path in config["coverage"]["beds"].items():
                    cov_file = os.path.join(config["coverage"]["output_dir"], f"{sample}.seq_summary")
                    colname = f"Mean_coverage_{bed_name}"
                    if os.path.exists(cov_file):
                        with open(cov_file) as f:
                            lines = [line.strip().split() for line in f if line.strip()]
                            if len(lines) >= 2 and len(lines[1]) > 2:
                                row[colname] = lines[1][2]
                            else:
                                row[colname] = "NA"
                    else:
                        row[colname] = "NA"

            # --- Consensus metrics ---
            if config.get("consensus", {}).get("run", False):
                for chr_name in config["consensus"]["chrs"]:
                    missing_file = os.path.join(config["consensus"]["output_dir_prefix"], f"{chr_name}/{sample}.missing.txt")
                    col_missing = f"Missing_{chr_name}"
                    if os.path.exists(missing_file):
                        with open(missing_file) as f:
                            lines = [line.strip().split("\t") for line in f if line.strip()]
                            if len(lines) >= 2 and len(lines[1]) >= 5:
                                row[col_missing] = lines[1][4]
                            else:
                                row[col_missing] = "NA"
                    else:
                        row[col_missing] = "NA"

            all_data.append(row)

        df = pd.DataFrame(all_data).fillna("NA")

        # Reorder columns: Sample first
        cols = ["Sample"] + [c for c in df.columns if c != "Sample"]
        df = df[cols]

        df.to_csv(OUTFILE, sep="\t", index=False)

