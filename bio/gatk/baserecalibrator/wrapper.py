__author__ = "Christopher Schröder"
__copyright__ = "Copyright 2020, Christopher Schröder"
__email__ = "christopher.schroeder@tu-dortmund.de"
__license__ = "MIT"


import os

extra = snakemake.params.get("extra", "")
java_opts = snakemake.params.get("java_opts", "")

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
known = snakemake.input.get("known", "")
if known:
    known = "--known-sites {}".format(known)

os.system(
    f"gatk --java-options '{java_opts}' BaseRecalibrator {extra} "
    f"-R {snakemake.input.ref} -I {snakemake.input.bam} "
    f"-O {snakemake.output.recal_table} {known} {log}"
)
