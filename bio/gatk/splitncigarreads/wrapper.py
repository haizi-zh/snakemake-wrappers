__author__ = "Jan Forster"
__copyright__ = "Copyright 2019, Jan Forster"
__email__ = "jan.forster@uk-essen.de"
__license__ = "MIT"

import os

import os

extra = snakemake.params.get("extra", "")
java_opts = snakemake.params.get("java_opts", "")

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
os.system(
    f"gatk --java-options '{java_opts}' SplitNCigarReads {extra} "
    f" -R {snakemake.input.ref} -I {snakemake.input.bam} "
    f"-O {snakemake.output} {log}"
)
