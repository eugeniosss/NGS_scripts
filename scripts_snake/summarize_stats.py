import os
from snakemake.shell import shell
from snakemake_wrapper_utils.samtools import get_samtools_opts
from pathlib import Path
import numpy as np

shell("echo ID >> {snakemake.params.tmp} ")
count_cols=1

if snakemake.params.settings !="False":
    shell("echo Raw_Reads >> {snakemake.params.tmp} "
          "&& echo Reads_surviving_trimming >> {snakemake.params.tmp} "
          "&& echo Merged_reads >> {snakemake.params.tmp} "
          "&& echo Reads_length_average >> {snakemake.params.tmp} ")
    count_cols+=4


shell(
    "echo Duplicates  >> {snakemake.params.tmp} " 
    "&& echo Duplicates_\(\%\)  >> {snakemake.params.tmp} "
    "&& echo Total_Reads >> {snakemake.params.tmp} " 
    "&& echo Mapped_reads >> {snakemake.params.tmp} "
    "&& echo Mapped_reads_\(\%\) >> {snakemake.params.tmp} "
    "&& echo MQ_30 >> {snakemake.params.tmp} "
    "&& echo MQ_30_\(\%\) >> {snakemake.params.tmp} "
    "&& echo Mapped_reads_on_MT >> {snakemake.params.tmp} "
    "&& echo Mapped_reads_on_MT\(\%\) >> {snakemake.params.tmp} "
    "&& echo MQ_30_on_MT >> {snakemake.params.tmp} "
    "&& echo MQ_30_on_MT\(\%\) >> {snakemake.params.tmp} ")
count_cols+=11

if snakemake.params.reads_on_target=="True":
    shell("echo Mapped_reads_on_Target >> {snakemake.params.tmp} "
          "&& echo Mapped_reads_on_Target_\(\%\) >> {snakemake.params.tmp} "
          "&& echo MQ_30_on_Target >> {snakemake.params.tmp} "
          "&& echo MQ_30_on_Target_\(\%\) >> {snakemake.params.tmp} ")
    count_cols+=4

if snakemake.params.fastas !="False":
    shell(
    "echo Coverage_on_MT >> {snakemake.params.tmp} "
    "&& echo Missing_data_on_MT >> {snakemake.params.tmp} ")
    count_cols+=2

if snakemake.params.nuclear_coverage!="False":
    shell(
    "echo Coverage_on_Nuclear >> {snakemake.params.tmp} ")	
    count_cols+=1


for code in snakemake.params.codes:
    shell ("echo {code} >> {snakemake.params.tmp} ")
    if snakemake.params.settings !="False":
        shell ("sed -n -e '/Total number of read pairs:/p' {snakemake.params.settings}/{code}.settings | awk '{{print $NF}}' >> {snakemake.params.tmp} "
               "&& sed -n -e '/Number of retained reads:/p' {snakemake.params.settings}/{code}.settings | awk '{{print $NF}}' >> {snakemake.params.tmp} "
               "&& sed -n -e '/Number of full-length collapsed pairs:/p' {snakemake.params.settings}/{code}.settings | awk '{{print $NF}}' >> {snakemake.params.tmp} "
               "&& sed -n -e '/Average length of retained reads:/p' {snakemake.params.settings}/{code}.settings | awk '{{print $NF}}' >> {snakemake.params.tmp} ")
    shell("sed -n '8p' {snakemake.params.pic}/{code}.metrics.txt | awk -F' ' '{{print $6}}' >> {snakemake.params.tmp} "
          "&& sed -n '8p' {snakemake.params.pic}/{code}.metrics.txt | awk -F' ' '{{print $9}}' >> {snakemake.params.tmp} "
          "&& sed -n '2p' {snakemake.params.stats}/{code}_basic_stats.txt | awk -F' ' '{{print $1}}' >> {snakemake.params.tmp} "
          "&& sed -n '2p' {snakemake.params.stats}/{code}_basic_stats.txt | awk -F' ' '{{print $2}}' >> {snakemake.params.tmp} "
          "&& sed -n '2p' {snakemake.params.stats}/{code}_basic_stats.txt | awk -F' ' '{{print $3}}' >> {snakemake.params.tmp} "
          "&& sed -n '2p' {snakemake.params.stats}/{code}_basic_stats.txt | awk -F' ' '{{print $4}}' >> {snakemake.params.tmp} "
          "&& sed -n '2p' {snakemake.params.stats}/{code}_basic_stats.txt | awk -F' ' '{{print $5}}' >> {snakemake.params.tmp} "
          "&& sed -n '2p' {snakemake.params.stats}/{code}_basic_stats.txt | awk -F' ' '{{print $6}}' >> {snakemake.params.tmp} "
          "&& sed -n '2p' {snakemake.params.stats}/{code}_basic_stats.txt | awk -F' ' '{{print $7}}' >> {snakemake.params.tmp} "
          "&& sed -n '2p' {snakemake.params.stats}/{code}_basic_stats.txt | awk -F' ' '{{print $8}}' >> {snakemake.params.tmp} "
          "&& sed -n '2p' {snakemake.params.stats}/{code}_basic_stats.txt | awk -F' ' '{{print $9}}' >> {snakemake.params.tmp} ")

    if snakemake.params.reads_on_target=="True":
        shell("sed -n '2p' {snakemake.params.stats}/{code}_basic_stats.txt | awk -F' ' '{{print $10}}' >> {snakemake.params.tmp} "
          "&& sed -n '2p' {snakemake.params.stats}/{code}_basic_stats.txt | awk -F' ' '{{print $11}}' >> {snakemake.params.tmp} "
          "&& sed -n '2p' {snakemake.params.stats}/{code}_basic_stats.txt | awk -F' ' '{{print $12}}' >> {snakemake.params.tmp} "
          "&& sed -n '2p' {snakemake.params.stats}/{code}_basic_stats.txt | awk -F' ' '{{print $13}}' >> {snakemake.params.tmp} ")


    if snakemake.params.fastas!="False":
         shell("sed -n '2p' {snakemake.params.MT_coverage}/{code}_MT.sample_summary | awk -F' ' '{{print $3}}' >> {snakemake.params.tmp} "
               "&& sed -n '2p' {snakemake.params.fastas}/{code}.log | awk '{{print $NF}}' >> {snakemake.params.tmp} ")
    if snakemake.params.nuclear_coverage!="False":
         shell("sed -n '2p' {snakemake.params.nuclear_coverage}/{code}_nuclear.sample_summary | awk -F' ' '{{print $3}}' >> {snakemake.params.tmp} ")

cols=" ".join(np.repeat("-", count_cols))

shell("paste -d '  ' {cols} < {snakemake.params.tmp} > {snakemake.output.final} && rm {snakemake.params.tmp}")

