__author__ = "Johannes Köster, Derek Croote"
__copyright__ = "Copyright 2020, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import os
import tempfile

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

outdir = os.path.dirname(snakemake.output[0])
if outdir:
    outdir = "--outdir {}".format(outdir)

extra = snakemake.params.get("extra", "")

with tempfile.TemporaryDirectory() as tmp:
    os.system(
        f"fasterq-dump --temp {tmp} --threads {snakemake.threads} "
        f"{extra} {outdir} {snakemake.wildcards.accession} {log}"
    )
