"""Snakemake wrapper for picard SortSam."""

__author__ = "Julian de Ruiter"
__copyright__ = "Copyright 2017, Julian de Ruiter"
__email__ = "julianderuiter@gmail.com"
__license__ = "MIT"


import os


extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

os.system(
    f"picard"
    f" SortSam"
    f" {extra}"
    f" INPUT={snakemake.input[0]}"
    f" OUTPUT={snakemake.output[0]}"
    f" SORT_ORDER={snakemake.params.sort_order}"
    f" {log}"
)
