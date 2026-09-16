__author__ = "Johannes Köster, Derek Croote"
__copyright__ = "Copyright 2020, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import os
import tempfile
from snakemake.shell import shell
from snakemake_wrapper_utils.snakemake import get_mem

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
extra = snakemake.params.get("extra", "")


# Parse memory
mem_mb = get_mem(snakemake, "MiB")


# Outdir
outdir = os.path.dirname(snakemake.output[0])
if outdir:
    outdir = f"--outdir {outdir}"


# Output compression
compress = ""
mem = f"-m{mem_mb}" if mem_mb else ""

for output in snakemake.output:
    out_name, out_ext = os.path.splitext(output)
    if out_ext == ".gz":
        compress += f"pigz -p {snakemake.threads} {out_name}; "
    elif out_ext == ".bz2":
        compress += f"pbzip2 -p{snakemake.threads} {mem} {out_name}; "


with tempfile.TemporaryDirectory() as tmpdir:
    mem = f"--mem {mem_mb}M" if mem_mb else ""

    accession = snakemake.wildcards.accession
    sra_dir = os.path.join(tmpdir, accession)
    fasterq_tmp = os.path.join(tmpdir, "fasterq_tmp")

    os.makedirs(fasterq_tmp)

    shell(
    "(prefetch {accession} -O {tmpdir} --max-size u || "
    "prefetch {accession} -O {tmpdir} --max-size u || "
    "prefetch {accession} -O {tmpdir} --max-size u; "
    "fasterq-dump {sra_dir} "
    "--temp {fasterq_tmp} "
    "--threads {snakemake.threads} {mem} "
    "{extra} {outdir}; "
    "{compress}"
    ") {log}"
)