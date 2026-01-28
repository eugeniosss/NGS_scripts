rule summarize_run_metrics:
    input:
        basic_stats = expand(
            config["stats"]["output_dir"] + "/{seq}_basic_stats.txt",
            seq=SEQS.keys()
        ),
        coverage = lambda wildcards: (
            expand(
                config["coverage"]["output_dir"] + "/{seq}.seq_summary",
                seq=SEQS.keys()
            ) if config.get("coverage", {}).get("run", False) else []
        ),
        dups = expand(
            config["dedup"]["output_dir"] + "/{seq}.metrics.txt",
            seq=SEQS.keys()
        ),
        adapter = lambda wildcards: (
            expand(
                config["adapterremoval2"]["output_dir"] + "/{seq}.settings",
                seq=SEQS.keys()
            ) if config.get("adapterremoval2", {}).get("run", False) else []
        ),
        consensus = lambda wildcards: (
            expand(
                config["consensus"]["output_dir_prefix"] + "{chr}/{seq}.missing.txt",
                seq=SEQS.keys(),
                chr=config["consensus"]["chrs"]
            ) if config.get("consensus", {}).get("run", False) else []
        )
    output:
        summary = "runs_summary.tsv"
    threads: 1
    resources:
        mem_mb = 500
    run:
        import pandas as pd
        import os

        OUTFILE = output.summary

        all_data = []
    
        for seq in SEQS.keys():
            row = {"Run": seq}
            
            # --- AdapterRemoval2 stats ---
            if config.get("adapterremoval2", {}).get("run", False):
                adapter_file = config["adapterremoval2"]["output_dir"] + f"/{seq}.settings"
                if os.path.exists(adapter_file):
                    metrics_of_interest = OrderedDict([
                        ("Number of retained reads", "Retained_reads"),
                        ("Number of full-length collapsed pairs", "Full_length_collapsed"),
                        ("Average length of retained reads", "Avg_length_retained")
                    ])
                    
                    # Read all lines once
                    with open(adapter_file) as f:
                        lines = [line.strip() for line in f if line.strip()]
                    
                    for key_text, col_name in metrics_of_interest.items():
                        # Search for the line starting with key_text
                        for line in lines:
                            if line.startswith(key_text):
                                row[col_name] = line.split()[-1]
                                break
                        else:
                            # If not found, mark as NA
                            row[col_name] = "NA"
    
    
            # --- Duplicates ---
            dups_file = config["dedup"]["output_dir"] + f"/{seq}.metrics.txt"
            if os.path.exists(dups_file):
                with open(dups_file) as f:
                    lines = [l.strip() for l in f if l.strip() != "" and not l.startswith("#")]
                    data_line = None
                    for i, line in enumerate(lines):
                        if line.startswith("LIBRARY"):
                            if i + 1 < len(lines):
                                data_line = lines[i + 1]
                            break
                    if data_line:
                        values = data_line.split("\t")
                        row["Dups_amount"] = values[5]           # PERCENT_DUPLICATION
                        row["Dups_rate"] = values[8]       # ESTIMATED_LIBRARY_SIZE
                    else:
                        row["dup_rate"] = "NA"
                        row["library_size"] = "NA"
            else:
                row["dup_rate"] = "NA"
                row["library_size"] = "NA"
    
            # --- Basic stats ---
            basic_file = config["stats"]["output_dir"] + f"/{seq}_basic_stats.txt"
            if os.path.exists(basic_file):
                with open(basic_file) as f:
                    lines = [line.strip().split() for line in f.readlines()]
                    if len(lines) >= 2:
                        for k, v in zip(lines[0], lines[1]):
                            row[k] = v
        
            # --- Coverage metrics ---
            if config.get("coverages", {}).get("run", False):
        
                # Each BED has its own summary file
                for bed_name, bed_path in config["coverages"]["beds"].items():
        
                    cov_file = (
                        config["coverages"]["output_dir_prefix"]
                        + f"{bed_name}/"
                        + f"{seq}.sample_summary"
                    )
        
                        # Column name becomes: Mean_coverage_all, Mean_coverage_exome …
                    colname = f"Mean_coverage_{bed_name}"
        
                    if os.path.exists(cov_file):
        
                        # ordered: line 2, column 3 → mean coverage (your example)
                        metrics_of_interest = OrderedDict([
                            (2, (colname, 2))
                        ])
        
                        with open(cov_file) as f:
                            lines = [line.strip().split() for line in f if line.strip()]
        
                        for line_num, (col_name, col_idx) in metrics_of_interest.items():
                            if len(lines) >= line_num and len(lines[line_num-1]) > col_idx:
                                row[col_name] = lines[line_num-1][col_idx]
                            else:
                                row[col_name] = "NA"
                    else:
                        row[colname] = "NA"
            
            # --- Consensus stats ---
            if config.get("consensus", {}).get("run", False):
                
                for chr_name in config["consensus"]["chrs"]:
                
                    missing_file = (
                        config["consensus"]["output_dir_prefix"]
                        + f"{chr_name}/"
                        + f"{seq}.missing.txt"
                    )
                
                    col_missing = f"Missing_{chr_name}"
                
                    if os.path.exists(missing_file):
                
                        with open(missing_file) as f:
                            header = f.readline().strip().split("\t")
                            values = f.readline().strip().split("\t")
                
                            if len(values) >= 5:  # ensure missing_pct is present
                                missing_pct = values[4]
                                row[col_missing] = missing_pct
                            else:
                                row[col_missing] = "NA"
                
                    else:
                        row[col_missing] = "NA"
            all_data.append(row)
    
        # Combine all rows
        df = pd.DataFrame(all_data).fillna("NA")

        # Desired first columns
        first_cols = ["Run", "Raw_Reads_R1", "Raw_Reads_R2"]

        # Columns that actually exist in the DataFrame
        existing_first = [c for c in first_cols if c in df.columns]

        # All remaining columns
        other_cols = [c for c in df.columns if c not in existing_first]

        # Reorder
        df = df[existing_first + other_cols]
        df.to_csv(OUTFILE, sep="\t", index=False)