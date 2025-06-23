🧬 Ancient DNA Mapping and Analysis Pipeline

This Snakemake pipeline processes ancient DNA sequencing data: it performs read preprocessing, mapping to a reference genome, duplicate removal, mitochondrial genome assembly, coverage analysis, and generates summary statistics.

⚙️ Pipeline Overview

🔹 Preprocessing
Adapter trimming with AdapterRemoval2

Merge overlapping reads

🔹 Mapping
Index circularized reference genome

Align reads using BWA

Realign to account for circular MT genome structure

🔹 Post-mapping
Sort BAM files

Mark duplicates using Picard

Index BAMs

🔹 Statistics and QC
Extract coverage metrics (GATK3)

Create FASTA consensus (HTSBox)

Generate summary tables for libraries and samples

🚀 Running the Pipeline
1. 📦 Clone the Project
git clone https://github.com/eugeniosss/NGS_scripts

2. 🧪 Create a Snakemake Environment
conda create -n snakemake -c conda-forge -c bioconda snakemake
conda activate snakemake

3. ▶️ Run with Snakemake (e.g., using 20 cores)
snakemake --use-conda --cores 20
