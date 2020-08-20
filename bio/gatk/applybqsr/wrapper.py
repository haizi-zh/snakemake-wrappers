__author__ = "Christopher Schröder"
__copyright__ = "Copyright 2020, Christopher Schröder"
__email__ = "christopher.schroeder@tu-dortmund.de"
__license__ = "MIT"


import os

extra = snakemake.params.get("extra", "")
java_opts = snakemake.params.get("java_opts", "")

log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)
os.system(
    f"gatk --java-options '{java_opts}' ApplyBQSR {extra} -R {snakemake.input.ref} -I {snakemake.input.bam} "
    f"--bqsr-recal-file {snakemake.input.recal_table} "
    f"-O {snakemake.output.bam} {log}"
)
