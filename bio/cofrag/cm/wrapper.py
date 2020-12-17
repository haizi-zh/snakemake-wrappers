__author__ = "Patrik Smeds"
__copyright__ = "Copyright 2016, Patrik Smeds"
__email__ = "patrik.smeds@gmail.com"
__license__ = "MIT"

from os import path
from snakemake.shell import shell
from tempfile import TemporaryDirectory
from snakemake.logging import logger


with TemporaryDirectory() as tempdir:
    job_label = snakemake.params.get("label")
    if job_label:
        logger.info(f"Job started: {job_label}")

    sinput = snakemake.input
    sparams = snakemake.params

    sample_id = sparams.get("sample_id")
    resol = sparams.get("resol", "")
    metrics = sparams.get("metrics", "")
    bootstrap = sparams.get("bootstrap", "")
    subsample = sparams.get("subsample", "")
    chroms = sparams.get("chroms", "")
    exclude_chroms = sparams.get("exclude_chroms", "")
    seed = sparams.get("seed", "")
    min_mapq = sparams.get("min_mapq", "")
    min_fraglen = sparams.get("min_fraglen", "")
    max_fraglen = sparams.get("max_fraglen", "")

    def join_regions(regions):
        if regions is None:
            return ""
        elif isinstance(regions, str):
            return regions
        else:
            return ":".join(regions)

    intersect_region = join_regions(sinput.get("intersect_region"))
    exclude_region = join_regions(sinput.get("exclude_region"))
    chrom_sizes = sinput.get("chrom_sizes", "")

    cmd = (
        f"Rscript {sparams.get('script')} "
        + f"--input {sinput.frag} "
        + f"--output-dir {tempdir} "
        + f"--sample-id {sample_id} "
        + f"--ncores {snakemake.threads} "
        + (f"--res {resol} " if resol else "")
        + (f"--metrics {metrics} " if metrics else "")
        + (f"--bootstrap {bootstrap} " if bootstrap else "")
        + (f"--subsample {subsample} " if subsample else "")
        + (f"--chroms {chroms} " if chroms else "")
        + (f"--exclude-chroms {exclude_chroms} " if exclude_chroms else "")
        + (f"--seed {seed} " if seed else "")
        + (f"--min-mapq {min_mapq} " if min_mapq else "")
        + (f"--min-fraglen {min_fraglen} " if min_fraglen else "")
        + (f"--max-fraglen {max_fraglen} " if max_fraglen else "")
        + (f"--intersect-region {intersect_region} " if intersect_region else "")
        + (f"--exclude-region {exclude_region} " if exclude_region else "")
        + (f"--chrom-sizes {chrom_sizes} " if chrom_sizes else "")
    )
    logger.info(cmd)

    shell("{cmd} " + ("2>&1 | tee {snakemake.log}" if (snakemake.log) else "2>&1"))

    output_dir = path.dirname(snakemake.output.cm)
    logger.info(f"Copying to destination: {output_dir}")
    shell("mv {tempdir}/{sample_id}.cofrag_cm.bed.gz {snakemake.output.cm}")