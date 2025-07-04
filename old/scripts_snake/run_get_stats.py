import os
from snakemake.shell import shell
from snakemake_wrapper_utils.samtools import get_samtools_opts
from pathlib import Path

samtools_opts = get_samtools_opts(snakemake)
extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
code=Path(snakemake.input.bam).stem

if snakemake.params.bed =="False":
    shell("touch tmp_{code}.txt && "
    "echo Total_Reads >> tmp_{code}.txt && "
    "echo Mapped_reads >> tmp_{code}.txt && "
    "echo Mapped_reads_\(\%\) >> tmp_{code}.txt && "
    "echo MQ_30 >> tmp_{code}.txt && "
    "echo MQ_30_\(\%\) >> tmp_{code}.txt && "
    "echo Mapped_reads_on_MT >> tmp_{code}.txt && "
    "echo Mapped_reads_on_MT\(\%\) >> tmp_{code}.txt && "
    "echo MQ_30_on_MT >> tmp_{code}.txt && "
    "echo MQ_30_on_MT\(\%\) >> tmp_{code}.txt && "

    "TR=$(samtools view -c {snakemake.input.bam}) && "
    "echo $TR >> tmp_{code}.txt && "

    "MAP=$(samtools view -c -F 260 {snakemake.input.bam}) && "
    "echo $MAP >> tmp_{code}.txt && "

    'bc <<< "scale=4 ; $MAP / $TR" >> tmp_{code}.txt && '

    "MQ30=$(samtools view -h -b -q 30 -c {snakemake.input.bam}) && "
    "echo $MQ30 >> tmp_{code}.txt && "

    'bc <<< "scale=4 ; $MQ30 / $TR" >> tmp_{code}.txt && '

    "MAPMT=$(samtools view -h -b -c -F 260 {snakemake.input.bam} MT) && "
    "echo $MAPMT >> tmp_{code}.txt && "

    'bc <<< "scale=4 ; $MAPMT / $MAP" >> tmp_{code}.txt && '

    "MAPMT30=$(samtools view -h -b -c -q 30 {snakemake.input.bam} MT) && "
    "echo $MAPMT30 >> tmp_{code}.txt && "

    'bc <<< "scale=4 ; $MAPMT30 / $MAP" >> tmp_{code}.txt && '

    "paste -d '  ' - - - - - - - - - < tmp_{code}.txt > stats/{code}_basic_stats.txt && "
    "rm tmp_{code}.txt")

else: 
    shell(
    "touch tmp_{code}.txt && "
    "echo Total_Reads >> tmp_{code}.txt && "
    "echo Mapped_reads >> tmp_{code}.txt && "
    "echo Mapped_reads_\(\%\) >> tmp_{code}.txt && "
    "echo MQ_30 >> tmp_{code}.txt && "
    "echo MQ_30_\(\%\) >> tmp_{code}.txt && "
    "echo Mapped_reads_on_MT >> tmp_{code}.txt && "
    "echo Mapped_reads_on_MT\(\%\) >> tmp_{code}.txt && "
    "echo MQ_30_on_MT >> tmp_{code}.txt && "
    "echo MQ_30_on_MT\(\%\) >> tmp_{code}.txt && "
    "echo Mapped_reads_on_Target >> tmp_{code}.txt && "
    "echo Mapped_reads_on_Target\(\%\) >> tmp_{code}.txt && "
    "echo MQ_30_on_Target >> tmp_{code}.txt && "
    "echo MQ_30_on_Target\(\%\) >> tmp_{code}.txt && "

    "TR=$(samtools view -c {snakemake.input.bam}) && "
    "echo $TR >> tmp_{code}.txt && "

    "MAP=$(samtools view -c -F 260 {snakemake.input.bam}) && "
    "echo $MAP >> tmp_{code}.txt && "

    'bc <<< "scale=4 ; $MAP / $TR" >> tmp_{code}.txt && '

    "MQ30=$(samtools view -h -b -q 30 -c {snakemake.input.bam}) && "
    "echo $MQ30 >> tmp_{code}.txt && "

    'bc <<< "scale=4 ; $MQ30 / $TR" >> tmp_{code}.txt && '

    "MAPMT=$(samtools view -h -b -c -F 260 {snakemake.input.bam} MT) && "
    "echo $MAPMT >> tmp_{code}.txt && "

    'bc <<< "scale=4 ; $MAPMT / $MAP" >> tmp_{code}.txt && '

    "MAPMT30=$(samtools view -h -b -c -q 30 {snakemake.input.bam} MT) && "
    "echo $MAPMT30 >> tmp_{code}.txt && "

    'bc <<< "scale=4 ; $MAPMT30 / $MAP" >> tmp_{code}.txt && '

    "MAPTARGET=$(samtools view -h -b -c -F 260 {snakemake.input.bam} -L {snakemake.params.bed} ) && "
    "echo $MAPTARGET >> tmp_{code}.txt && "

    'bc <<< "scale=4 ; $MAPTARGET / $MAP" >> tmp_{code}.txt && '

    "MAPTARGET30=$(samtools view -h -b -c -q 30 {snakemake.input.bam} -L {snakemake.params.bed} ) && "
    "echo $MAPTARGET30 >> tmp_{code}.txt && "

    'bc <<< "scale=4 ; $MAPTARGET30 / $MAP" >> tmp_{code}.txt && '

    "paste -d '  ' - - - - - - - - - - - - - < tmp_{code}.txt > stats/{code}_basic_stats.txt && "
    "rm tmp_{code}.txt")
