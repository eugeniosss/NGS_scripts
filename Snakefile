import csv
import re
import os

### INPUT ###

configfile:"config.yaml"
csvfile=csv.reader(open(config["input"]))
SEQS = dict(csvfile)
DIC_MERGE = {}

### HELPER VARIABLES ###

adapter_config = config.get("adapterremoval2", {})
DO_ADAPTER = adapter_config.get("run", False)
DO_COLLAPSE = DO_ADAPTER and ("--collapse" in adapter_config.get("options", ""))
DO_CONCAT = adapter_config.get("concatenate", False)
DO_CIRC = config["circular_mapper"]["run"]
CLEANUP = config["clean_intermediates"]
COUNT_PP = config["stats"]["count_properly_paired"]

### NON SENSE ERRORS ###


if DO_CONCAT and not DO_ADAPTER:
    raise ValueError(
        "adapterremoval2.concatenate=True requires adapterremoval2.run=True"
    )

if len(config["mapping"]["reads_group_parser"])!=1 and config["mapping"]["reads_group_parser"]!=False:
        raise ValueError(
            f"ERROR: mapping.reads_group_parser must be length=1 for bwa aln, "
            f"but got '{config['mapping']['reads_group_parser']}' (length={len(config['mapping']['reads_group_parser'])})."
        )

if config["mapping"]["soft"]=="bwa-mem":
    if DO_CIRC:
        raise ValueError(
            f"ERROR: mapper bwa-mem is not compatible with circular_mapper option.\n"
            f"Please either use bwa mem with no circular_mapper (fine for modern) or bwa aln with circular_mapper (suggested for ancient)."
        )
elif config["mapping"]["soft"]!="bwa-aln":
    raise ValueError(
        f"ERROR: mapper should be 'bwa-mem' or 'bwa-aln'.\n"
        f"Invalid option used: '{config['mapping']['soft']}"
    ) 

if COUNT_PP:
    if DO_COLLAPSE:
        raise ValueError(
            "stats.count_properly_paired=True is incompatible with AdapterRemoval --collapse "
            "(collapsed reads are single-end)."
        )

    if DO_CONCAT:
        raise ValueError(
            "stats.count_properly_paired=True is incompatible with concatenated FASTQs "
            "(concatenation destroys pairing information)."
        )

if not (isinstance(config["stats"]["add_chrs"], list) or config["stats"]["add_chrs"] is False):
    raise Exception("Error: config['stats']['add_chrs'] must be a list or False")

if not (isinstance(config["stats"]["add_beds"], dict) or config["stats"]["add_beds"] is False):
    raise Exception("Error: config['stats']['add_beds'] must be a dict or False")


### HELPER FUNCTIONS ###

def get_final_fastqs(seq):
    """
    Return the list of FASTQ files to map for a given SEQ.
    
    Hierarchy:
    1. Concatenate outputs if concatenate=True
        - If AdapterRemoval2 collapsed -> concatenate collapsed R1 and R2
        - If not collapsed -> concatenate surviving R1 and R2
    2. If no concatenation:
        - AdapterRemoval2 outputs (collapsed or surviving)
    3. If AdapterRemoval2 not run:
        - Raw fastqs from dic.txt
    """
    if DO_CONCAT:
        if DO_ADAPTER and DO_COLLAPSE:
            return f"{adapter_config['output_dir']}/{seq}.concatenated.gz"
        elif DO_ADAPTER and not DO_COLLAPSE:
            return [f"{adapter_config['output_dir']}/{seq}.pair1.truncated.gz",
                    f"{adapter_config['output_dir']}/{seq}.pair2.truncated.gz"]
        else:
            return SEQS[seq].split()
    else:
        if DO_ADAPTER:
            if DO_COLLAPSE:
                return f"{adapter_config['output_dir']}/{seq}.collapsed.gz"
            else:
                return [f"{adapter_config['output_dir']}/{seq}.pair1.truncated.gz",
                          f"{adapter_config['output_dir']}/{seq}.pair2.truncated.gz"]
        else:
            return SEQS[seq].split()

def extract_sample(seq):
    """
    Extract sample name from SEQ using the configured parser.
    If parser is False, return the SEQ unchanged.
    """
    parser = config["mapping"].get("reads_group_parser", False)

    if not parser:
        raise ValueError(
            f"ERROR: to merge runs of same bams, need to assign a reads_group_parser in mapping step!"
        )

    return seq.split(parser)[0]

def get_fastq_r1(unit):
    # Only return FASTQ if the unit is an individual SEQ (not merged)
    if unit in SEQS:
        fq = SEQS[unit].split(" ")[0]
        if os.path.exists(fq):
            return fq
    return ""  # merged samples: no fastq

def get_fastq_r2(unit):
    if unit in SEQS:
        fq2 = SEQS[unit].split(" ")[1] if " " in SEQS[unit] else ""
        if fq2 and os.path.exists(fq2):
            return fq2
    return ""

def get_final_bam(unit):
    """
    Return BAM to use for downstream steps (basic_stats, coverage, etc.)
    """
    if unit in DIC_MERGE:
        if DO_CIRC:
            return f"{config['merge_same_sample_runs']['output_circular_dir']}/{unit}.bam"
        else:
            return f"{config['merge_same_sample_runs']['output_dedup_dir']}/{unit}.bam"
    if DO_CIRC:
        return f"{config['circular_mapper']['output_dir']}/{unit}.bam"
    return f"{config['dedup']['output_dir']}/{unit}.bam"

def get_final_bai(unit):
    return get_final_bam(unit) + ".bai"

def get_reference_for_mapping():
    """
    Return the path to the reference FASTA that should be used for mapping.

    - If circular_mapper is enabled, return the circularized reference.
    - Otherwise, return the normal reference.
    """
    if DO_CIRC:
        REF_RAW = config["ref"].rsplit(".", 1)[0]
        REF_EXT = config["ref"].rsplit(".", 1)[1]  
        ext = str(config["circular_mapper"]["elongation"])
        return f"{REF_RAW}_{ext}.{REF_EXT}"
    else:
        return config["ref"]

def maybe_temp(path):
    return temp(path) if CLEANUP else path

def temp_if_not_final(path):
    """
    Return temp(path) if circular mapping is enabled,
    otherwise return path as-is (final file).
    """
    if DO_CIRC:
        return temp(path)
    return path

### Build processing units ###

# Individual SEQ units
SEQ_UNITS = list(SEQS.keys())
SAMPLE_UNITS = []
FINAL_OUTPUTS = []

if config.get("merge_same_sample_runs", {}).get("run", False):

    PARSER = config["mapping"]["reads_group_parser"]

    if not PARSER:
        raise ValueError("mapping.reads_group_parser must be set")

    for seq in SEQS.keys():
        parts = seq.split(PARSER)
        if len(parts) != 3:
            raise ValueError(
                f"ERROR: SEQS must be parsed in 3 substrings by '{PARSER}' "
                f"(SAMPLE{PARSER}LIBRARY{PARSER}ID). "
                f"Invalid SEQ: '{seq}'"
            )

    P = re.escape(PARSER)

    DIC_MERGE = {}
    for seq_key in SEQS.keys():
        sample = extract_sample(seq_key)
        bam = f"{config['dedup']['output_dir']}/{seq_key}.bam"

        DIC_MERGE.setdefault(sample, []).append(bam)

    # Only include samples for which bam files actually exist
    MERGE_SAMPLES = sorted(DIC_MERGE.keys())

    if DO_CIRC:
        FINAL_OUTPUTS += expand(
            f"{config['merge_same_sample_runs']['output_circular_dir']}" + "/{sample}.bam",
            sample=MERGE_SAMPLES
        )
        FINAL_OUTPUTS += expand(
            config['merge_same_sample_runs']['output_circular_dir'] + "/{sample}.bam.bai",
            sample=MERGE_SAMPLES
        )

    else:
        FINAL_OUTPUTS += expand(
            f"{config['merge_same_sample_runs']['output_dedup_dir']}" + "/{sample}.bam",
            sample=MERGE_SAMPLES
        )
        FINAL_OUTPUTS += expand(
            config['merge_same_sample_runs']['output_dedup_dir'] + "/{sample}.bam.bai",
            sample=MERGE_SAMPLES
        )

    FINAL_OUTPUTS += ["samples_summary.tsv"]

    SAMPLE_UNITS = list(DIC_MERGE.keys())

    include: config["dir"]+"/rules/PicardMerge.smk"

# Combined list: everything downstream will run on these
UNITS = SEQ_UNITS + SAMPLE_UNITS

### OUTPUTS GENERATION ###

#Final Stats
FINAL_OUTPUTS += expand(
    "{output_dir}/{unit}_basic_stats.txt",
    unit=UNITS,
    output_dir=config["stats"]["output_dir"]
)

#Final non merged bams
FINAL_OUTPUTS += [
    get_final_bam(unit)
    for unit in SEQS.keys()
]

#Final non merged bam indexes
FINAL_OUTPUTS += [
    get_final_bai(unit)
    for unit in SEQS.keys()
]

### IMPORT RULES

if DO_CIRC:
    include: config["dir"] + "/rules/CircularRef.smk"

include: config["dir"]+"/rules/IndexRef.smk"

include: config["dir"]+"/rules/SamtoolsSort.smk"

include: config["dir"]+"/rules/PicardMarkDuplicates.smk"

include: config["dir"]+"/rules/get_seq_metrics.smk"

include: config["dir"]+"/rules/summarize_metrics.smk"

if DO_ADAPTER:
    if DO_COLLAPSE:
        include: config["dir"] + "/rules/run_AdapterRemoval2_SE.smk"
    else:
        include: config["dir"] + "/rules/run_AdapterRemoval2_PE.smk"

if DO_CONCAT:
    include: config["dir"] + "/rules/concatenate.smk"

if config["mapping"]["soft"]=="bwa-aln":
    if DO_CONCAT or DO_COLLAPSE:
        include: config["dir"] + "/rules/run_bwa_aln_SE.smk"
    else:
        include: config["dir"] + "/rules/run_bwa_aln_PE.smk"

if config["mapping"]["soft"]=="bwa-mem":
    include: config["dir"]+"/rules/run_bwa_mem.smk"

if config["consensus"]["run"]==True:
    if config["consensus"]["soft"]=="htsbox":
        include: config["dir"]+"/rules/ConcensusHtsbox.smk"
        FINAL_OUTPUTS += expand(
            config["consensus"]["output_dir_prefix"] + "{chr}/{unit}.fasta",
            chr=config["consensus"]["chrs"],
            unit=UNITS
        )

    elif config["consensus"]["soft"]=="samtools":
        include: config["dir"]+"/rules/ConcensusSamtools.smk"
        FINAL_OUTPUTS += expand(
            config["consensus"]["output_dir_prefix"] + "{chr}/{unit}.fasta",
            chr=config["consensus"]["chrs"],
            unit=UNITS
        )

if config["coverages"]["run"]==True:
    include: config["dir"]+"/rules/CoverageGatk3.smk"
    FINAL_OUTPUTS += expand(
        # The structure of the output files
        config["coverages"]["output_dir_prefix"] + "{bedtype}/{unit}.{extension}",
        
        # 1. Wildcards for bed type and sequence
        bedtype=config["coverages"]["beds"].keys(),
        unit=UNITS,
        
        # 2. Wildcard for the file extension/suffix
        extension=[
            "sample_summary",
            "sample_interval_summary",
            "sample_cumulative_coverage_counts",
            "sample_cumulative_coverage_proportions",
            "sample_interval_statistics",
            "sample_statistics"
        ]
    )

rule all:
    input:
       FINAL_OUTPUTS,
       "runs_summary.tsv"
