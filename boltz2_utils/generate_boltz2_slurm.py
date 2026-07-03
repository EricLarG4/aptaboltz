from pathlib import Path

import yaml


SLURM_TEMPLATE = """\
#!/bin/bash

#############################
# SLURM OPTIONS
#############################
#SBATCH -J {job_name}
#SBATCH -t 23:59:59
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --gres gpu:1
#SBATCH -p gpu-l40,gpu-a40,gpu-rtx6000
#SBATCH --mem=64G
#############################

echo "#############################"
echo "User:" $USER
echo "Date:" $(date)
echo "Host:" $(hostname)
echo "PWD:" $(pwd)
echo "SLURM_JOBID:" $SLURM_JOBID
echo "#############################"


#############################
# IDENTIFIERS
#############################
SUFFIX=_J${{SLURM_JOBID}}

###############################
# DIRECTORIES (SAFE ISOLATION)
###############################
SCRATCHDIR=${{SLURM_SUBMIT_DIR}}/{cluster_path}
RUNROOT="${{SCRATCHDIR}}/${{SLURM_JOB_NAME}}/J${{SLURM_JOBID}}"
mkdir -p "${{RUNROOT}}"

#############################
# TEMP WORKDIR (COMPUTE NODE)
#############################
TMPBASE=${{TMPDIR:-/tmp}}
WORKDIR="${{TMPBASE}}/boltz${{SUFFIX}}"
mkdir -p "${{WORKDIR}}"
cd "${{WORKDIR}}"

echo "Workdir : $WORKDIR"
echo "Out root: $RUNROOT"

#############################
# STAGE INPUT FILES
#############################
echo "Staging input files"

cp ${{SCRATCHDIR}}/yaml/{job_name}.yaml input.yaml

#############################
# METADATA
#############################
echo "${{SLURM_JOBID}} $(hostname)" > run.info

#############################
# RUN
#############################

# --recycling_steps 10 --diffusion_samples 25 is the default for AlphaFold3
# --use_potentials: inference time potential that significantly improve the physical quality of the poses

boltz predict input.yaml \\
    --out_dir output/ \\
    --accelerator gpu \\
\t--recycling_steps 12 \\
\t--sampling_steps 800 \\
\t--diffusion_samples 25 \\
\t--max_parallel_samples 5 \\
\t--step_scale 1.5 \\
\t--output_format mmcif \\
\t--num_workers 1 \\
{msa_flags}\t--use_potentials \\
\t--write_full_pae \\
\t--write_full_pde \\

#############################
# COPY BACK RESULTS
#############################
echo "Copying results to $RUNROOT"
cp -r output/. "${{RUNROOT}}/"     # copy contents of output/ (incl. boltz_results_*/)

echo "Boltz is done"
"""


# Two MSA-related flags that are omitted when the YAML already supplies an MSA.
_MSA_FLAGS = "\t--max_msa_seqs 8192 \\\n\t--use_msa_server \\\n"


def _yaml_has_msa(yaml_path: Path) -> bool:
    """Return True if any DNA entry in *yaml_path* contains an ``msa`` field."""
    with open(yaml_path) as fh:
        data = yaml.safe_load(fh)
    return any(
        "msa" in entry["dna"]
        for entry in data.get("sequences", [])
        if "dna" in entry
    )


def generate_boltz2_slurm_files(
    yaml_source: str | Path | list[Path],
    cluster_path: str,
    output_dir: str | Path | None = None,
) -> list[Path]:
    """Generate SLURM submission scripts for all YAML files in a directory.

    Accepts either:
    - A local directory (``str`` or ``Path``): scans ``{yaml_source}/yaml/``
      for ``*.yaml`` files.
    - A list of ``Path`` objects (e.g. the return value of
      ``generate_boltz2_yamls``).

    The SLURM job name and staged input filename are derived from each YAML
    stem (e.g. ``CSS1_C0R.yaml`` → job ``CSS1_C0R``).

    Parameters
    ----------
    yaml_source:
        Local project root directory or list of ``.yaml`` ``Path`` objects.
    cluster_path:
        Path to the project directory on the cluster, relative to
        ``$SLURM_SUBMIT_DIR`` (i.e. where you call ``sbatch`` from).
        For example: ``"scratch/boltz/projects/CSS"``.
        Embedded in each script as ``SCRATCHDIR=${SLURM_SUBMIT_DIR}/{cluster_path}``.
    output_dir:
        Where to write the ``.slurm`` files locally.  Defaults to the project
        root inferred from ``yaml_source``.  Created automatically if needed.

    Returns
    -------
    list[Path]
        Paths of the written ``.slurm`` files.
    """
    if isinstance(yaml_source, list):
        yaml_files = sorted(yaml_source)
        if not yaml_files:
            raise ValueError("The list of yaml paths is empty.")
        local_root = yaml_files[0].parent.parent
    else:
        local_root = Path(yaml_source)
        yaml_dir = local_root / "yaml"
        if not yaml_dir.exists():
            raise FileNotFoundError(f"YAML directory not found: {yaml_dir}")
        yaml_files = sorted(yaml_dir.glob("*.yaml"))
        if not yaml_files:
            raise FileNotFoundError(f"No .yaml files found in {yaml_dir}")

    out_dir = Path(output_dir) if output_dir else local_root
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nGenerating SLURM files for {len(yaml_files)} YAML file(s) in {local_root}...")

    written: list[Path] = []
    for yaml_path in yaml_files:
        job_name = yaml_path.stem
        msa_flags = "" if _yaml_has_msa(yaml_path) else _MSA_FLAGS
        content = SLURM_TEMPLATE.format(
            job_name=job_name,
            cluster_path=cluster_path.strip("/"),
            msa_flags=msa_flags,
        )
        # If the YAML filename contains "_custom_msa", propagate that suffix
        # to the SLURM filename so the two files stay paired.
        slurm_filename = f"{job_name}.slurm"
        slurm_path = out_dir / slurm_filename
        slurm_path.write_text(content, newline="\n")
        written.append(slurm_path)
        print(f"  ✓ {slurm_path}")

    print(f"\n{len(written)} SLURM file(s) written to {out_dir}\n")
    return written
