rule basic_stats:
    input:
        bam = lambda wc: get_final_bam(wc.unit),
        bai = lambda wc: get_final_bai(wc.unit)
    output:
        txt = config['stats']['output_dir'] + "/{unit}_basic_stats.txt"
    params:
        chrs = " ".join(config["stats"].get("add_chrs", [])), 
        bed_names = " ".join(config["stats"].get("add_beds", {}).keys()), 
        bed_files = " ".join(config["stats"].get("add_beds", {}).values()),
        raw_reads_r1 = lambda wc: SEQS[wc.unit].split()[0] if wc.unit in SEQS else "",
        raw_reads_r2 = lambda wc: SEQS[wc.unit].split()[1] if wc.unit in SEQS and " " in SEQS[wc.unit] else ""
    conda:
        config["dir"] + "envs/NGS.yml"
    log:
        config['stats']['output_dir'] + "/{unit}.basic_stats.log"
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
    
        # --- CHR-based stats ---
        for CHR in $LIST_CHRS; do
            echo Mapped_${{CHR}} >> $TMP
            echo Mapped_${{CHR}}_percent >> $TMP
            echo MQ30_${{CHR}} >> $TMP
            echo MQ30_${{CHR}}_percent >> $TMP
            count_cols=$((count_cols+4))
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
    
        # --- Per-CHR percentages ---
        for CHR in $LIST_CHRS; do
            MAPCHR=$(samtools view -c -F 260 $BAM $CHR)
            MQ30CHR=$(samtools view -c -q 30 $BAM $CHR)
            echo $MAPCHR >> $TMP
            if [ "$MAP" -eq 0 ]; then echo NA >> $TMP; else echo $((MAPCHR*100/MAP)) >> $TMP; fi
            echo $MQ30CHR >> $TMP
            if [ "$MAP" -eq 0 ]; then echo NA >> $TMP; else echo $((MQ30CHR*100/MAP)) >> $TMP; fi
        done
    
        # --- Per-BED percentages ---
        for i in $(seq 0 $((${{#BED_NAMES_ARR[@]}} - 1))); do
            PANEL=${{BED_NAMES_ARR[i]}}
            BED=${{BED_FILES_ARR[i]}}
            MAPPAN=$(samtools view -c -F 260 -L $BED $BAM)
            MQ30PAN=$(samtools view -c -q 30 -L $BED $BAM)
            echo $MAPPAN >> $TMP
            if [ "$MAP" -eq 0 ]; then echo NA >> $TMP; else echo $((MAPPAN*100/MAP)) >> $TMP; fi
            echo $MQ30PAN >> $TMP
            if [ "$MAP" -eq 0 ]; then echo NA >> $TMP; else echo $((MQ30PAN*100/MAP)) >> $TMP; fi
        done
    
        # --- Final formatting ---
        COLS=""
        for i in $(seq $count_cols); do
            COLS="${{COLS}} -"
        done
        paste -d " " $COLS < $TMP > $OUTFILE
        rm $TMP
        """