__author__ = "Patrik Smeds"
__copyright__ = "Copyright 2018, Patrik Smeds"
__email__ = "patrik.smeds@gmail.com"
__license__ = "MIT"


import os

input_flag = "--vcf"
if snakemake.input[0].endswith(".gz"):
    input_flag = "--gzvcf"

output = " > " + snakemake.output[0]
if output.endswith(".gz"):
    output = " | gzip" + output

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

extra = snakemake.params.get("extra", "")

os.system(
    f"vcftools "
    f"{input_flag} "
    f"{snakemake.input} "
    f"{extra} "
    f"--recode "
    f"--stdout "
    f"{output} "
    f"{log}"
)
