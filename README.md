🧬 Ancient DNA Mapping and Analysis Pipeline
This Snakemake pipeline processes ancient DNA sequencing data, performs read preprocessing, mapping to a reference genome, duplicate removal, mitochondrial genome assembly, coverage analysis, and generates summary statistics.

⚙️ Pipeline Overview
Preprocessing:

  Adapter trimming with AdapterRemoval2

  Merge overlapping reads

Mapping:

  Index circularized reference genome

  Align reads using BWA

  Realign to account for circular MT genome structure

Post-mapping:

  Sort BAMs

  Mark duplicates (Picard)

  Index BAMs

Statistics and QC:

  Extract coverage metrics (GATK3)

  Create FASTA consensus (HTSBox)

  Generate summary tables for libraries and samples

🚀 Running the Pipeline
Clone/Prepare Project
git clone https://github.com/eugeniosss/NGS_scripts

Install snakemake in a new environment
conda create -c conda-forge -c bioconda -n snakemake snakemake


Run with Snakemake (for example 20 cores)
snakemake --use-conda --cores 20
