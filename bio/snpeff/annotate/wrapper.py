__author__ = "Bradford Powell"
__copyright__ = "Copyright 2018, Bradford Powell"
__email__ = "bpow@unc.edu"
__license__ = "BSD"


import os
from os import path
import shutil
import tempfile
from pathlib import Path

outcalls = snakemake.output.calls
if outcalls.endswith(".vcf.gz"):
    outprefix = "| bcftools view -Oz"
elif outcalls.endswith(".bcf"):
    outprefix = "| bcftools view -Ob"
else:
    outprefix = ""

incalls = snakemake.input[0]
if incalls.endswith(".bcf"):
    incalls = "< <(bcftools view {})".format(incalls)

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

extra = snakemake.params.get("extra", "")

data_dir = Path(snakemake.input.db).parent.resolve()

stats = snakemake.output.get("stats", "")
csvstats = snakemake.output.get("csvstats", "")
csvstats_opt = "" if not csvstats else "-csvStats {}".format(csvstats)
stats_opt = "-noStats" if not stats else "-stats {}".format(stats)

reference = path.basename(snakemake.input.db)

os.system(
    f"snpEff -dataDir {data_dir} {stats_opt} {csvstats_opt} {extra} "
    f"{reference} {incalls} "
    f"{outprefix} > {outcalls} {log}"
)
