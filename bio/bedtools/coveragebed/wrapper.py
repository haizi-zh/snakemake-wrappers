__author__ = "Patrik Smeds"
__copyright__ = "Copyright 2019, Patrik Smeds"
__email__ = "patrik.smeds@gmail.com"
__license__ = "MIT"

import os


log = snakemake.log_fmt_shell(stdout=False, stderr=True)

extra_params = snakemake.params.get("extra", "")

input_a = snakemake.input.a
input_b = snakemake.input.b

output_file = snakemake.output[0]

if not isinstance(output_file, str) and len(snakemake.output) != 1:
    raise ValueError("Output should be one file: " + str(output_file) + "!")

os.system(
    f"coverageBed"
    f" -a {input_a}"
    f" -b {input_b}"
    f" {extra_params}"
    f" > {output_file}"
    f" {log}"
)
