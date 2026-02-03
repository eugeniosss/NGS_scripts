# 🧬 NGS Mapping & Analysis Pipeline (Snakemake)

This repository contains a **Snakemake-based NGS pipeline** designed for flexible processing of short-read sequencing data, with particular support for ancient DNA–style workflows (adapter trimming, collapsing, circular genomes, deduplication, consensus, statistics, and coverage).

The pipeline is fully configurable via a single YAML file.

## 📥 Input
<pre>input: dic.txt</pre>

dic.txt maps sample IDs to FASTQ files

Format example:

sample1,reads_R1.fastq.gz reads_R2.fastq.gz

## 📁 General Settings
<pre>dir: Directory where NGS_scripts is.

ref: Path to fasta file.

mem_overhead: Memory safety factor applied to memory-intensive rules. Useful for samtools sort, Picard, and GATK. Example: mem_mb: 8000 → actually requests 9600 MB

clean_intermediates: Automatically deletes intermediate files once downstream steps are completed. Helps reduce disk usage</pre>

## ✂️ Adapter Removal (AdapterRemoval2)
<pre>adapterremoval2:
  run: True or False
  output_dir: directory to save output
  threads: number of threads
  options: additional options
  concatenate: concatenate collapsed with r1 and r2 that passed filters. True or False
  mem_mb: memmory to use at this step (in MB)</pre>

⚠️ Important:

If --collapse or concatenate: True is used, proper pairing information is lost

This affects downstream statistics (see stats.count_properly_paired)

## 🧭 Mapping
<pre>mapping:
  soft: software for mapping. bwa-aln or bwa-mem
  output_dir: directory to save output
  threads: number of threads
  params: additional options
  reads_group_parser: Read group names are parsed using reads_group_parser. One character or False
  mapping_mem_mb: memmory to use at this step (in MB)
  indexing_mem_mb: memmory to use at this step (in MB)</pre>

## 🔃 Sorting
<pre>sorting:
  threads: number of threads
  mem_mb:  memmory to use at this step (in MB)</pre>

## 🔁 Circular Genome Mapping (Optional)
<pre>circular_mapper:
  run: True or False
  output_dir: directory to save output
  elongation: bps to elongate
  chr: chr to elongate
  mem_mb: memmory to use at this step (in MB)</pre>

## 🧹 Deduplication
<pre>dedup:
  output_dir: directory to save output
  params: additional options
  mem_mb: memmory to use at this step (in MB)</pre>

## 🔀 Merging Multiple Runs (Optional)
<pre>merge_same_sample_runs:
  run: True or False
  params: additional options
  output_merge_dir: directory to save merged output
  output_dedup_dir: directory to save merged output after dedup
  output_circular_dir: if circularmapper, directory to save merged output after dedup and circularmapper
  mem_mb: memmory to use at merging step (in MB)</pre>

## 🧬 Consensus Generation
<pre>consensus:
  run: True or False
  soft: software for consensus. htsbox or samtools
  chrs: chrs to be conensensus called ex ["MT"]
  output_dir_prefix: prefix directory to save output
  params: additional options
  mem_mb: memmory to use at this step (in MB)</pre>

## 📊 Basic Statistics
<pre>stats:
  output_dir: directory to save output
  count_properly_paired: True or False
  add_chrs: chrs to count statistics. ex ["MT"]
  add_beds: beds to count statistics (python dictionary or False)
    Exome: path/to/exome.bed or False</pre>

⚠️ Important constraints:

count_properly_paired: True is only valid if:

Data are paired-end

Reads were not collapsed

Reads were not concatenated

The pipeline will fail early if this condition is violated.

## 📈 Coverage Calculation
<pre>coverages:
  run: True or False
  threads: number of threads
  output_dir_prefix: prefix directory to save output
  params: additional options
  beds: beds to calculate coverage (python dictionary)
    all: /media/jbod2/eugenio/ref_genomes/myotragus_balearicus_NC_042943.fasta.bed
    exome: /media/jbod2/eugenio/ref_genomes/myotragus_balearicus_NC_042943.fasta.bed
  mem_mb : 8000</pre>

## ✅ Design Philosophy

Fail early on invalid configurations

Highly modular: turn steps on/off easily

Ancient DNA–friendly

Disk-aware via automatic cleanup

Reproducible via Conda environments

# 🚀 Running the Pipeline

## 1. 📦 Clone the Project

<pre>git clone https://github.com/eugeniosss/NGS_scripts
cd NGS_scripts</pre>

## 2. 🧪 Create a Snakemake Environment

<pre>conda create -n snakemake -c conda-forge -c bioconda snakemake
conda activate snakemake</pre>

## 3. 🧬 Prepare Input Files

Prepare Input Files

<pre>sample1 reads_R1.fastq.gz reads_R2.fastq.gz</pre>

## 4. ⚙️ Configure the Pipeline

Edit the configuration file to match your data and system:

<pre>nano config.yaml</pre>

Make sure to update:
- Input sample dictionary (`input`)
- Reference genome (`ref`)
- Output directories
- Enabled/disabled modules
- Memory and thread settings

## 5. ▶️ Run with Snakemake (example using 20 cores)

<pre>snakemake --use-conda --cores 20</pre>

## 6. 🔍 (Optional) Dry-run to Preview Jobs

<pre>snakemake -n</pre>
