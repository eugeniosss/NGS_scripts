__author__ = "Johannes Köster"
__copyright__ = "Copyright 2016, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


from snakemake.shell import shell
from snakemake_wrapper_utils.samtools import get_samtools_opts
import os


samtools_opts = get_samtools_opts(snakemake)
extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
sample=os.path.basename(str(snakemake.output))

if len(snakemake.input)==1:
        shell("cp {snakemake.input} {snakemake.output}")

elif snakemake.params.herit=="True":
        shell("touch rg_{sample}.txt && for bam in {snakemake.input};do samtools view -H $bam | grep '^@RG' | tail -n 1 >> rg_{sample}.txt; done && samtools merge -rh rg_{sample}.txt {snakemake.output} {snakemake.input} {log}")
else:
        shell("samtools merge {snakemake.output} {snakemake.input} {log}")


