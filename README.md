# 🧬 NGS Mapping & Analysis Pipeline (Snakemake)

This repository contains a **Snakemake-based NGS pipeline** designed for flexible processing of short-read sequencing data, with particular support for ancient DNA–style workflows (adapter trimming, collapsing, circular genomes, deduplication, consensus, statistics, and coverage).

The pipeline is fully configurable via a single YAML file.

📥 Input
input: dic.txt

dic.txt maps sample IDs to FASTQ files

Format example:

sample1   reads_R1.fastq.gz reads_R2.fastq.gz
sample2   reads.fastq.gz

📁 General Settings
dir: /media/jbod2/eugenio/NGS_scripts/
ref: /media/jbod2/eugenio/ref_genomes/myotragus_balearicus_NC_042943.fasta


dir: Base directory of the pipeline (used to locate rules, envs, scripts)

ref: Reference genome FASTA

mem_overhead: 1.2


Memory safety factor applied to memory-intensive rules

Useful for samtools sort, Picard, and GATK

Example: mem_mb: 8000 → actually requests 9600 MB

clean_intermediates: true


Automatically deletes intermediate files once downstream steps are completed

Helps reduce disk usage

✂️ Adapter Removal (AdapterRemoval2)
adapterremoval2:
  run: True
  output_dir: adapterremoval2
  threads: 4
  options: "--trimns --trimqualities --minlength 30"
  concatenate: False
  mem_mb: 2000


Uses AdapterRemoval2

Supports:

Adapter trimming

Quality trimming

Read collapsing (optional)

Concatenation of PE reads (optional)

⚠️ Important:

If --collapse or concatenate: True is used, proper pairing information is lost

This affects downstream statistics (see stats.count_properly_paired)

🧭 Mapping
mapping:
  soft: bwa-aln
  output_dir: mapping
  threads: 4
  params: "-l 1024"
  reads_group_parser: "_"
  mapping_mem_mb: 4000
  indexing_mem_mb: 10000


Supports BWA (bwa aln by default)

Read group names are parsed using reads_group_parser

Reference indexing and mapping memory can be tuned separately

🔃 Sorting
sorting:
  threads: 4
  mem_mb: 16000


Coordinate sorting using samtools sort

Memory automatically scaled using mem_overhead

🔁 Circular Genome Mapping (Optional)
circular_mapper:
  run: True
  output_dir: circular
  elongation: 500
  chr: MT
  mem_mb: 2000


Designed for circular genomes (e.g. mitochondria)

Artificially elongates reference to recover edge-spanning reads

Typically used for mitochondrial DNA (chr: MT)

🧹 Deduplication
dedup:
  output_dir: dedup
  params: "--REMOVE_DUPLICATES true --VALIDATION_STRINGENCY LENIENT --ASSUME_SORT_ORDER coordinate"
  mem_mb: 10000


Uses Picard MarkDuplicates

Designed to work with both linear and circular mappings

🔀 Merging Multiple Runs (Optional)
merge_same_sample_runs:
  run: False
  params: "--VALIDATION_STRINGENCY LENIENT"
  output_merge_dir: merged
  output_dedup_dir: merged_dedup
  output_circular_dir: merged_circular
  mem_mb: 10000


Merges BAMs from multiple sequencing runs of the same sample

Can be applied before or after deduplication

🧬 Consensus Generation
consensus:
  run: True
  soft: samtools
  chrs: ["MT"]
  output_dir_prefix: fastas_
  params: "--min-MQ 30 --min-BQ 30 -d 5"
  mem_mb: 2000


Generates consensus FASTA sequences

Typically used for mitochondrial or targeted regions

Fully configurable quality thresholds

📊 Basic Statistics
stats:
  output_dir: stats
  count_properly_paired: True
  add_chrs: ["MT"]
  add_beds:
    MT2: /media/jbod2/eugenio/ref_genomes/myotragus_balearicus_NC_042943.fasta.bed


The statistics module computes:

Total reads

Mapped reads and percentages

MQ30 reads and percentages

Optional properly paired statistics

Optional per-chromosome and per-BED statistics

⚠️ Important constraints:

count_properly_paired: True is only valid if:

Data are paired-end

Reads were not collapsed

Reads were not concatenated

The pipeline will fail early if this condition is violated.

📈 Coverage Calculation
coverages:
  run: True
  threads: 1
  output_dir_prefix: coverage_
  params: "--minMappingQuality 30 --omitDepthOutputAtEachBase"
  beds:
    all: /media/jbod2/eugenio/ref_genomes/myotragus_balearicus_NC_042943.fasta.bed
    exome: /media/jbod2/eugenio/ref_genomes/myotragus_balearicus_NC_042943.fasta.bed
  mem_mb : 8000


Computes depth and coverage statistics

Supports multiple BED regions

Outputs are organized per region

✅ Design Philosophy

Fail early on invalid configurations

Highly modular: turn steps on/off easily

Ancient DNA–friendly

Disk-aware via automatic cleanup

Reproducible via Conda environments
