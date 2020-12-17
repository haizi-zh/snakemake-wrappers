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

    def build_standards(ids, files):
        if not ids or not files:
            return ""

        if isinstance(ids, str):
            ids = [ids]
        if isinstance(files, str):
            files = [files]

        assert len(ids) == len(files)

        return ":".join([f"{v1}:{v2}" for v1, v2 in list(zip(ids, files))])

    standard_hic = build_standards(
        sparams.get("standard_hic"), sinput.get("standard_hic")
    )
    standard_compartment = build_standards(
        sparams.get("standard_compartment"), sinput.get("standard_compartment")
    )
    chrom_sizes = sinput.get("chrom_sizes")
    juicer = sinput.get("juicer")

    cmd = (
        f"Rscript {sparams.get('script')} "
        + f"--input {sinput.cm} "
        + f"--output-dir {tempdir} "
        + f"--sample-id {sample_id} "
        + (f"--standard-hic {standard_hic} " if standard_hic else "")
        + (
            f"--standard-compartment {standard_compartment} "
            if standard_compartment
            else ""
        )
        + (f"--juicer {juicer} " if juicer else "")
        + (f"--chrom-sizes {chrom_sizes} " if chrom_sizes else "")
    )
    logger.info(cmd)

    shell("{cmd} " + ("2>&1 | tee {snakemake.log}" if (snakemake.log) else "2>&1"))

    output_dir = path.dirname(snakemake.output.compartment)
    logger.info(f"Copying to destination: {output_dir}")
    shell(
        "mv {tempdir}/{sample_id}.compartment.bedGraph.gz {snakemake.output.compartment}"
    )
    if snakemake.output.get("compcor") and snakemake.output.get("compcor_chart"):
        shell(
            "mv {tempdir}/{sample_id}.compartment-correlation.pdf {snakemake.output.compcor_chart}"
        )
        shell(
            "mv {tempdir}/{sample_id}.compartment-correlation.txt {snakemake.output.compcor}"
        )
    if snakemake.output.get("hicrep") and snakemake.output.get("hicrep_chart"):
        shell("mv {tempdir}/{sample_id}.hicrep.pdf {snakemake.output.hicrep_chart}")
        shell("mv {tempdir}/{sample_id}.hicrep.txt {snakemake.output.hicrep}")
