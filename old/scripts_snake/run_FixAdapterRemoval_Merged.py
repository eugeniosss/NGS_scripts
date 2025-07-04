import os
from snakemake.shell import shell
from snakemake_wrapper_utils.samtools import get_samtools_opts
from pathlib import Path

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell("AdapterRemovalFixPrefix {snakemake.input} | sed '/^[@M_+]/ s/M_//g' > {snakemake.output}")
