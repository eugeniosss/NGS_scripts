k__email__ = "lauri.mesilaakso@gmail.com"
__license__ = "MIT"

import tempfile
from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts
from os import path

java_opts = get_java_opts(snakemake)
extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# Extract basename from the output file names
out_basename = path.commonprefix(snakemake.output).rstrip(".")

with tempfile.TemporaryDirectory() as tmpdir:
    shell(
        "gatk3 {java_opts}"
        " -T DepthOfCoverage"
        " -I {snakemake.input.bam}"
        " -L {snakemake.input.intervals}"
        " -R {snakemake.input.fasta}"
        " -o {out_basename}"
        " {extra}"
        " {log}"
    )
