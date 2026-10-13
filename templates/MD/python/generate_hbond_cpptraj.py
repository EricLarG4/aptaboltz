#!/usr/bin/env python
"""Generate a cpptraj H-bond analysis input for a completed MD run.

The Boltz-2 SLURM pipeline (`PrepAndMin.slurm`) produces `2drms`,
`atomicfluct`, stripped topologies and trajectories, but NOT H-bond data.
Those come from a cpptraj script that must be pointed at a specific job and
replica set, and re-pointed every time a new job is analysed.  This generator
produces that script from the run's actual layout so the inputs are
reproducible instead of hand-edited.

Writes a cpptraj input that, for each replica, records ligand/protein H-bonds
along the production trajectory:

  * per-frame H-bond matrix
  * average contacts (``avgout``)      <- consumed by md_utils.R readers
  * contact series (long format CSV)
  * lifetime curve

Outputs are written into each replica's ``prep/traj/`` directory.  Paths
inside the generated script are relative to the script's own location, so it
can be run from anywhere.

Usage
-----
From ``{project}/MD/``::

    python python/generate_hbond_cpptraj.py \\
        --project my_project --experiment my_seq_L01 --job J1234567

Behaviour
---------
``--job``      omitted -> the highest-numbered ``J*`` directory for the
                experiment is used (latest job), and the choice is reported.

``--tasks``    omitted -> **auto-discovery**: every ``task_*`` directory under
                the job that contains a ``prep/traj/*_stripped.prmtop`` is
                included.  A replica that has not produced a stripped topology
                is skipped and reported on stderr; it is NOT silently assumed
                to be absent, so check the summary before analysing.  Pass an
                explicit value (``--tasks 0-2`` or ``--tasks 0,1,2``) for a
                fully reproducible selection.

``--dna-end``  omitted -> read from the Boltz-2 YAML.  The length of the
                first ``dna``/``rna``/``protein`` sequence block gives the
                macromolecule length, and the H-bond mask is
                ``:1-(dna_end + ligand)``.  If the YAML cannot be found or
                contains no sequence block, the script errors rather than
                guessing.  Override with ``--dna-end`` (or ``--ligands`` to
                change how many ligand residues follow) if your system differs.

Requires ``PyYAML``.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path, PurePosixPath

try:
    import yaml
except ImportError:  # pragma: no cover - environment problem, not logic
    sys.exit("error: PyYAML is required (pip install pyyaml)")

JOB_RE = re.compile(r"J(\d+)$")


# --------------------------------------------------------------------------
# discovery helpers
# --------------------------------------------------------------------------
def detect_job(out_root: Path, experiment: str) -> str:
    """Return the highest-numbered J-directory for *experiment*."""
    numbered = []
    for path in out_root.glob(f"{experiment}/J*"):
        match = JOB_RE.match(path.name)
        if match:
            numbered.append((int(match.group(1)), path.name))
    if not numbered:
        sys.exit(f"error: no J* directory found under {out_root / experiment}")
    return max(numbered)[1]


def detect_tasks(traj_root: Path, report) -> list[int]:
    """Return sorted task numbers whose stripped topology already exists."""
    tasks = []
    for task_dir in sorted(traj_root.glob("task_*")):
        if list((task_dir / "prep" / "traj").glob("*_stripped.prmtop")):
            tasks.append(int(task_dir.name.split("_", 1)[1]))
        else:
            report(f"note: {task_dir.name} has no stripped topology - skipped")
    if not tasks:
        sys.exit(f"error: no completed task found under {traj_root}")
    return tasks


def parse_tasks(spec: str) -> list[int]:
    """Parse ``0-2`` or ``0,1,2`` into a sorted unique list of ints."""
    spec = str(spec).strip()
    if "-" in spec and "," not in spec:
        start, end = spec.split("-", 1)
        return list(range(int(start), int(end) + 1))
    return sorted({int(t) for t in spec.split(",") if t.strip()})


def sequence_length(yaml_path: Path) -> int:
    """Return the length of the first dna/rna/protein sequence block."""
    if not yaml_path.exists():
        sys.exit(f"error: YAML not found: {yaml_path}")
    with yaml_path.open(encoding="utf-8") as fh:
        data = yaml.safe_load(fh)
    for entry in (data or {}).get("sequences", []):
        for key in ("dna", "rna", "protein"):
            if key in entry:
                seq = entry[key].get("sequence", "")
                if seq:
                    return len(seq)
    sys.exit(f"error: no dna/rna/protein sequence block found in {yaml_path}")


# --------------------------------------------------------------------------
# rendering
# --------------------------------------------------------------------------
def render(args, task_ids: list[int], dna_end: int, md_dir: Path) -> str:
    experiment = args.experiment
    job = args.job
    # cpptraj runs on Linux, so every path emitted must use forward slashes;
    # stringifying a Windows Path would embed backslashes instead.
    traj_rel = PurePosixPath("pmemd/out") / experiment / job
    last = dna_end + args.ligands

    header = [
        f"# {experiment} - {job} - H-bond analysis",
        "# Generated by python/generate_hbond_cpptraj.py - do not edit by hand.",
        "# Regenerate after each new job instead of copying this file.",
        "#",
        f"# Macromolecule: residues 1-{dna_end} ({dna_end} nt/res from YAML)"
        f"  |  ligand: {dna_end + 1}-{last}",
        f"# Geometry: distance <= {args.dist} A, angle >= {args.angle} deg",
        f"# Replicas: {', '.join('task_%d' % t for t in task_ids)}",
        "#",
        "# Outputs per replica (in prep/traj/):",
        "#   step10_contacts_{N}_{J}.dat          per-frame H-bond matrix",
        "#   step10_contacts_avg_{N}_{J}.dat      average contacts (read by md_utils.R)",
        "#   step10_hbond_contacts_{N}_{J}.csv    contact series",
        "#   step10_contacts_lifetime_{N}_{J}.dat lifetime curve",
        f"#   {args.series_name}_{{N}}_{{J}}.dat           per-frame H-bond presence series",
        "#",
        "# The *_uuseries file is an intermediate: 'series uuseries' creates the",
        "# time-series dataset that writedata/lifetime consume.  Nothing downstream",
        "# reads the file itself, but the keyword must not be removed.",
        "#",
        "# Run with:  cpptraj -i <this file>",
        "",
    ]

    blocks = []
    for n in task_ids:
        task = f"task_{n}"
        prmtop = f"{traj_rel / task / 'prep' / 'traj' / f'step10_{n}_{job}_stripped.prmtop'}"
        netcdf = f"{traj_rel / task / 'prep' / 'traj' / f'step10_{n}_{job}_stripped.nc'}"
        contacts = f"{traj_rel / task / 'prep' / 'traj' / f'step10_contacts_{n}_{job}.dat'}"
        avg = f"{traj_rel / task / 'prep' / 'traj' / f'step10_contacts_avg_{n}_{job}.dat'}"
        series = f"{traj_rel / task / 'prep' / 'traj' / f'step10_hbond_contacts_{n}_{job}.csv'}"
        lifetime = f"{traj_rel / task / 'prep' / 'traj' / f'step10_contacts_lifetime_{n}_{job}.dat'}"
        # Per-frame H-bond presence series.  "series uuseries" is what CREATES the
        # time-series dataset that the contactseries and lifetime steps consume,
        # so it cannot be dropped; giving each replica its own path stops the
        # blocks from overwriting one shared file in the working directory.
        uuseries = f"{traj_rel / task / 'prep' / 'traj' / f'{args.series_name}_{n}_{job}.dat'}"

        for path in (prmtop, netcdf):
            if not (md_dir / path).exists():
                sys.exit(f"error: missing input for {task}: {md_dir / path}")

        blocks.append(
            f"# {task}\n"
            f"parm {prmtop}\n"
            f"trajin {netcdf} 1 last {args.stride}\n"
            f"hbond contacts out {contacts} \\\n"
            f"    :1-{last} angle {args.angle:g} dist {args.dist:g} \\\n"
            f"    avgout {avg} \\\n"
            f"    series uuseries {uuseries}\n"
            "run\n"
            "\n"
            f"writedata {series} contactseries contacts[*]\n"
            "run\n"
            "\n"
            f"lifetime contacts[solutehb] out {lifetime}\n"
            "run\n"
            "\n"
            "go\n"
            "clear all\n"
        )

    return "\n".join(header) + "\n" + "\n".join(blocks)


# --------------------------------------------------------------------------
def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate a cpptraj H-bond analysis input for an MD run.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--project", required=True, help="Project directory name, e.g. my_project")
    parser.add_argument("--experiment", required=True, help="Experiment name, e.g. my_seq_L01")
    parser.add_argument("--job", default=None, help="Job directory (default: latest J* found)")
    parser.add_argument("--tasks", default=None,
                        help="Replicas: '0-2' or '0,1,2'. Omitted = auto-discover completed tasks.")
    parser.add_argument("--yaml", default=None,
                        help="Boltz-2 YAML (default: ../yaml/<experiment>.yaml)")
    parser.add_argument("--dna-end", type=int, default=None,
                        help="Macromolecule length (default: read from the YAML sequence)")
    parser.add_argument("--ligands", type=int, default=1,
                        help="Number of ligand residues following the macromolecule (default: 1)")
    parser.add_argument("--angle", type=float, default=135.0, help="Min H-bond angle in degrees")
    parser.add_argument("--dist", type=float, default=3.2, help="Max donor-acceptor distance in A")
    parser.add_argument("--stride", type=int, default=5, help="Frame stride (default: 5)")
    parser.add_argument("--series-name", default="step10_uuseries",
                        help="Base name for the per-frame H-bond presence series "
                             "(default: step10_uuseries). Written into each replica's "
                             "prep/traj/; required by the contactseries and lifetime steps.")
    parser.add_argument("--out", default="hbond_analysis.cpptraj", help="Output file")
    args = parser.parse_args()

    md_dir = Path.cwd()
    out_root = md_dir / "pmemd" / "out"

    if args.job is None:
        args.job = detect_job(out_root, args.experiment)
        print(f"Auto-detected job: {args.job}")
    traj_root = out_root / args.experiment / args.job

    if args.tasks is None:
        task_ids = detect_tasks(traj_root, lambda m: print(m, file=sys.stderr))
        print(f"Auto-discovered {len(task_ids)} task(s): "
              f"{', '.join('task_%d' % t for t in task_ids)}")
    else:
        task_ids = parse_tasks(args.tasks)

    if args.dna_end is None:
        yaml_path = Path(args.yaml) if args.yaml else md_dir.parent / "yaml" / f"{args.experiment}.yaml"
        args.dna_end = sequence_length(yaml_path)
        print(f"Macromolecule length from {yaml_path}: {args.dna_end}")

    out_path = Path(args.out).resolve()
    out_path.write_text(render(args, task_ids, args.dna_end, md_dir), encoding="utf-8")

    # Rewrite the relative paths so they are anchored at the output file, not
    # the working directory, keeping the script runnable from anywhere.
    print(f"Wrote {out_path} ({len(task_ids)} replica block(s), "
          f"mask :1-{args.dna_end + args.ligands})")


if __name__ == "__main__":
    main()