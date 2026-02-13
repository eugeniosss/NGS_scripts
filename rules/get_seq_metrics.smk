rule basic_stats:
    input:
        bam = lambda wc: get_final_bam(wc.unit),
        bai = lambda wc: get_final_bai(wc.unit)
    output:
        txt = config['stats']['output_dir'] + "/{unit}_basic_stats.txt"
    params:
        chrs = lambda wc: " ".join(config["stats"].get("add_chrs", []) if config["stats"].get("add_chrs") else []),
        bed_names = lambda wc: " ".join(config["stats"].get("add_beds", {}).keys() if config["stats"].get("add_beds") else []),
        bed_files = lambda wc: " ".join(config["stats"].get("add_beds", {}).values() if config["stats"].get("add_beds") else []),
        raw_reads_r1 = lambda wc: SEQS[wc.unit].split()[0] if wc.unit in SEQS else "",
        raw_reads_r2 = lambda wc: SEQS[wc.unit].split()[1] if wc.unit in SEQS and " " in SEQS[wc.unit] else "",
        count_pp = config["stats"]["count_properly_paired"]
    conda:
        config["dir"] + "envs/NGS.yml"
    log:
        config['stats']['output_dir'] + "/{unit}.basic_stats.log"
    resources:
        mem_mb = 1000
    shell:
        r"""
        BAM="{input.bam}"
        OUTFILE="{output.txt}"
    
        FQ1="{params.raw_reads_r1}"
        FQ2="{params.raw_reads_r2}"
    
        LIST_CHRS="{params.chrs}"
        LIST_BED_NAMES="{params.bed_names}"
        LIST_BED_FILES="{params.bed_files}"
    
        TMP="tmp_{wildcards.unit}.txt"
        count_cols=0
        
        # --- Header & RAW READS ---
        if [ -f "$FQ1" ]; then
            echo Raw_Reads_R1 > $TMP
            if [[ "$FQ1" == *.gz ]]; then
                RAW1=$(zcat "$FQ1" | wc -l)
            else
                RAW1=$(cat "$FQ1" | wc -l)
            fi
            RAW1=$((RAW1/4))
            count_cols=$((count_cols+1))
        fi
        
        if [ -n "$FQ2" ] && [ -f "$FQ2" ]; then
            echo Raw_Reads_R2 >> $TMP
            if [[ "$FQ2" == *.gz ]]; then
                RAW2=$(zcat "$FQ2" | wc -l)
            else
                RAW2=$(cat "$FQ2" | wc -l)
            fi
            RAW2=$((RAW2/4))
            count_cols=$((count_cols+1))
        fi
        
        # --- BAM header columns ---
        echo Total_Reads >> $TMP
        echo Mapped_reads >> $TMP
        echo Mapped_reads_percent >> $TMP
        echo MQ30_reads >> $TMP
        echo MQ30_reads_percent >> $TMP
        count_cols=$((count_cols+5))

        if [ "{params.count_pp}" != "False" ]; then
            echo Properly_paired_reads >> $TMP
            echo Properly_paired_reads_percent >> $TMP
            echo Properly_paired_MQ30 >> $TMP
            echo Properly_paired_MQ30_percent >> $TMP
            count_cols=$((count_cols+4))
        fi
    
        # --- CHR-based stats ---
        for CHR in $LIST_CHRS; do
            echo Mapped_${{CHR}} >> $TMP
            echo Mapped_${{CHR}}_percent >> $TMP
            echo MQ30_${{CHR}} >> $TMP
            echo MQ30_${{CHR}}_percent >> $TMP
            count_cols=$((count_cols+4))
            if [ "{params.count_pp}" != "False" ]; then
                echo Properly_paired_${{CHR}} >> $TMP
                echo Properly_paired_${{CHR}}_percent >> $TMP
                echo Properly_paired_MQ30_${{CHR}} >> $TMP
                echo Properly_paired_MQ30_${{CHR}}_percent >> $TMP
                count_cols=$((count_cols+4))
            fi
        done
    
        # --- BED-based stats ---
        BED_NAMES_ARR=($LIST_BED_NAMES)
        BED_FILES_ARR=($LIST_BED_FILES)
    
        for i in $(seq 0 $((${{#BED_NAMES_ARR[@]}} - 1))); do
            PANEL=${{BED_NAMES_ARR[i]}}
            echo Mapped_${{PANEL}} >> $TMP
            echo Mapped_${{PANEL}}_percent >> $TMP
            echo MQ30_${{PANEL}} >> $TMP
            echo MQ30_${{PANEL}}_percent >> $TMP
            count_cols=$((count_cols+4))
            if [ "{params.count_pp}" != "False" ]; then
                echo Properly_paired_${{PANEL}} >> $TMP
                echo Properly_paired_${{PANEL}}_percent >> $TMP
                echo Properly_paired_MQ30_${{PANEL}} >> $TMP
                echo Properly_paired_MQ30_${{PANEL}}_percent >> $TMP
                count_cols=$((count_cols+4))
            fi
        done
    
        # --- BAM counts ---
        if [ -f "$FQ1" ]; then
            echo $RAW1 >> $TMP
        fi
        
        if [ -n "$FQ2" ] && [ -f "$FQ2" ]; then
            echo $RAW2 >> $TMP
        fi

        TR=$(samtools view -c $BAM)
        MAP=$(samtools view -c -F 260 $BAM)
        MQ30=$(samtools view -c -q 30 $BAM)
    
        echo $TR >> $TMP
        echo $MAP >> $TMP
        if [ "$TR" -eq 0 ]; then echo NA >> $TMP; else echo "scale=4; $MAP/$TR" | bc >> $TMP; fi
        echo $MQ30 >> $TMP
        if [ "$TR" -eq 0 ]; then echo NA >> $TMP; else echo "scale=4; $MQ30/$TR" | bc >> $TMP; fi

        if [ "{params.count_pp}" != "False" ]; then
            PP=$(samtools view -c -f 2 $BAM)
            PP_MQ30=$(samtools view -c -f 2 -q 30 $BAM)
            echo $PP >> $TMP
            if [ "$MAP" -eq 0 ]; then echo NA >> $TMP; else echo "scale=4; $PP/$MAP" | bc >> $TMP; fi
            echo $PP_MQ30 >> $TMP
            if [ "$MAP" -eq 0 ]; then echo NA >> $TMP; else echo "scale=4; $PP_MQ30/$MAP" | bc >> $TMP; fi            
        fi
    
        # --- Per-CHR percentages ---
        for CHR in $LIST_CHRS; do
            MAPCHR=$(samtools view -c -F 260 $BAM $CHR)
            MQ30CHR=$(samtools view -c -q 30 $BAM $CHR)
            echo $MAPCHR >> $TMP
            if [ "$MAP" -eq 0 ]; then echo NA >> $TMP; else echo "scale=4; $MAPCHR/$MAP" | bc >> $TMP; fi
            echo $MQ30CHR >> $TMP
            if [ "$MAP" -eq 0 ]; then echo NA >> $TMP; else echo "scale=4; $MQ30CHR/$MAP" | bc >> $TMP; fi

            if [ "{params.count_pp}" != "False" ]; then
                PPCHR=$(samtools view -c -f 2 $BAM $CHR)
                PPCHR_MQ30=$(samtools view -c -f 2 -q 30 $BAM $CHR)

                echo $PPCHR >> $TMP
                if [ "$MAP" -eq 0 ]; then echo NA >> $TMP; else echo "scale=4; $PPCHR/$MAP" | bc >> $TMP; fi
                echo $PPCHR_MQ30 >> $TMP
                if [ "$MAP" -eq 0 ]; then echo NA >> $TMP; else echo "scale=4; $PPCHR_MQ30/$MAP" | bc >> $TMP; fi
            fi

        done
    
        # --- Per-BED percentages ---
        for i in $(seq 0 $((${{#BED_NAMES_ARR[@]}} - 1))); do
            PANEL=${{BED_NAMES_ARR[i]}}
            BED=${{BED_FILES_ARR[i]}}
            MAPPAN=$(samtools view -c -F 260 -L $BED $BAM)
            MQ30PAN=$(samtools view -c -q 30 -L $BED $BAM)
            echo $MAPPAN >> $TMP
            if [ "$MAP" -eq 0 ]; then echo NA >> $TMP; else echo "scale=4; $MAPPAN/$MAP" | bc >> $TMP; fi
            echo $MQ30PAN >> $TMP
            if [ "$MAP" -eq 0 ]; then echo NA >> $TMP; else echo "scale=4; $MQ30PAN/$MAP" | bc >> $TMP; fi

            if [ "{params.count_pp}" != "False" ]; then
                PPPAN=$(samtools view -c -f 2 -L $BED $BAM)
                PPPAN_MQ30=$(samtools view -c -f 2 -q 30 -L $BED $BAM)

                echo $PPPAN >> $TMP
                if [ "$MAP" -eq 0 ]; then echo NA >> $TMP; else echo "scale=4; $PPPAN/$MAP" | bc >> $TMP; fi
                echo $PPPAN_MQ30 >> $TMP
                if [ "$MAP" -eq 0 ]; then echo NA >> $TMP; else echo "scale=4; $PPPAN_MQ30/$MAP" | bc >> $TMP; fi
            fi

        done
    
        # --- Final formatting ---
        COLS=""
        for i in $(seq $count_cols); do
            COLS="${{COLS}} -"
        done
        paste -d " " $COLS < $TMP > $OUTFILE
        rm $TMP
        """
