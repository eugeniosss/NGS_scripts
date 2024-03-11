import os
from snakemake.shell import shell
from snakemake_wrapper_utils.samtools import get_samtools_opts


samtools_opts = get_samtools_opts(snakemake)
extra = snakemake.params.get("extra", "")

reference=str("-f "+snakemake.input.reference)
bed=str("-b "+snakemake.input.bed)

shell(
    "(htsbox pileup "
    "{reference} "
    "{bed} "
    "{extra} "
    "{snakemake.input.bam} > {snakemake.output} "
    ") > {snakemake.log} 2>&1 ")
#    "&& echo Missing data percentage: >> {snakemake.log} "
#    "&& N=$(seqtk comp {snakemake.output} | awk '{{x+=$9}}END{{print x}}')"
#    "&& T=$(seqtk comp {snakemake.output} | awk '{{x+=$2}}END{{print x}}') "
#    "&& echo 'scale =3; $N/$T' | bc -l >> {snakemake.log} ")


dnaSequence=[]
f1=open(str(snakemake.output))
for line in f1:
    line = line.rstrip()
    if line.startswith('>'):
        header=line
    else:
        dnaSequence += line
nCount = 0
for c in dnaSequence:
    if c == 'N' or c == "n":
        nCount = nCount + 1
sequenceLength = len(dnaSequence)

if nCount !=0:
    Ns=round(float(nCount) / sequenceLength,3)
    shell("echo Missing data percentage: {Ns} >> {snakemake.log} ")
else:
    shell("echo Missing data percentage: 1.0 >> {snakemake.log} ")
    shell("echo [M::write_fa] average depth for contig \'MT\': -nan >> {snakemake.log} ")

file=str(snakemake.output).split('/')[-1]
seq=str(file).split('.')[0]

shell("sed -i 's/>.*/>{seq}/g' {snakemake.output}")
