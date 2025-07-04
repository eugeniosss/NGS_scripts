rule basic_stats:
    input:
        bam = config['dedup']['output_dir']+"/{SEQ}.bam",
        bai = config['dedup']['output_dir']+"/{SEQ}.bam.bai"
    output:
        config['stats']['output_dir']+"/{SEQ}_basic_stats.txt"
    params:
        MT = config["MT_stats"]["run"],
        MT_bed = config["MT_stats"]["bed"],
        target = config["Target_region_stats"]["run"],
        target_bed = config["Target_region_stats"]["bed"]
    threads: 1
    log:
        config['stats']['output_dir']+"/{SEQ}.basic_stats.log"
    conda:
        config["dir"] + "envs/NGS.yml"
    shell:
        """
        (TMP=tmp_{wildcards.SEQ}.txt
        OUTFILE={output}

        # Header
        echo Total_Reads > $TMP
        echo Mapped_reads >> $TMP
        echo Mapped_reads_\(\%\) >> $TMP
        echo MQ_30 >> $TMP
        echo MQ_30_\(\%\) >> $TMP

        count_cols=5
        
        if [ {params.MT} == True ]; then
            echo Mapped_reads_on_MT >> $TMP
            echo Mapped_reads_on_MT\(\%\) >> $TMP
            echo MQ_30_on_MT >> $TMP
            echo MQ_30_on_MT\(\%\) >> $TMP
            count_cols=$((count_cols+4))
        fi

        if [ {params.target} == True ]; then
            echo Mapped_reads_on_Target >> $TMP
            echo Mapped_reads_on_Target\(\%\) >> $TMP
            echo MQ_30_on_Target >> $TMP
            echo MQ_30_on_Target\(\%\) >> $TMP
            count_cols=$((count_cols+4))
        fi

        BAM={input.bam}
        TR=$(samtools view -c $BAM)
        echo $TR >> $TMP

        MAP=$(samtools view -c -F 260 $BAM)
        echo $MAP >> $TMP


        if [ "$TR" -eq 0 ]; then
            echo "NA" >> "$TMP"
        else
            echo "scale=4 ; $MAP / $TR" | bc >> "$TMP"
        fi

        MQ30=$(samtools view -h -b -q 30 -c $BAM)
        echo $MQ30 >> $TMP


        if [ "$TR" -eq 0 ]; then
            echo "NA" >> "$TMP"
        else
            echo "scale=4 ; $MQ30 / $TR" | bc >> "$TMP"
        fi

       if [ {params.MT} == True ]; then
            MAPMT=$(samtools view -h -b -c -F 260 $BAM -L {params.MT_bed})
            echo $MAPMT >> $TMP

            if [ "$MAP" -eq 0 ]; then
                echo "NA" >> "$TMP"
            else
                echo "scale=4 ; $MAPMT / $MAP" | bc >> "$TMP"
            fi

            MAPMT30=$(samtools view -h -b -c -q 30 $BAM -L {params.MT_bed})
            echo $MAPMT30 >> $TMP

            if [ "$MAP" -eq 0 ]; then
                echo "NA" >> "$TMP"
            else
                echo "scale=4 ; $MAPMT30 / $MAP" | bc >> "$TMP"
            fi

        fi

        if [ {params.target} == True ]; then
            MAPTARGET=$(samtools view -h -b -c -F 260 $BAM -L {params.target_bed})
            echo $MAPTARGET >> $TMP

            if [ "$MAP" -eq 0 ]; then
                echo "NA" >> "$TMP"
            else
                echo "scale=4 ; $MAPTARGET / $MAP" | bc >> "$TMP"
            fi

            MAPTARGET30=$(samtools view -h -b -c -q 30 $BAM -L {params.target_bed})
            echo $MAPTARGET30 >> $TMP

            if [ "$MAP" -eq 0 ]; then
                echo "NA" >> "$TMP"
            else
                echo "scale=4 ; $MAPTARGET30 / $MAP" | bc >> "$TMP"
            fi

        fi

        COLS=""
        for i in $(seq $count_cols); do
            COLS="${{COLS}} -"
        done

        paste -d ' ' $COLS < $TMP > $OUTFILE

        #rm $TMP 
        )&> {log}
        """