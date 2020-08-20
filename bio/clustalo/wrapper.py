"""Snakemake wrapper for clustal omega."""

__author__ = "Michael Hall"
__copyright__ = "Copyright 2019, Michael Hall"
__email__ = "mbhall88@gmail.com"
__license__ = "MIT"


import os

# Placeholder for optional parameters
extra = snakemake.params.get("extra", "")
# Formats the log redrection string
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# Executed shell command
os.system(
    f"clustalo {extra}"
    f" --threads={snakemake.threads}"
    f" --in {snakemake.input[0]}"
    f" --out {snakemake.output[0]} "
    f" {log}"
)
