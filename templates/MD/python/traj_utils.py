"""
Trajectory visualisation utilities for Amber MD outputs.

Provides functions for generating interactive Mol* viewer files (MVSJ)
from Amber MD trajectory data (prmtop + netcdf trajectory pairs).

Functions
---------
find_trajectory_pairs(job_dir)
    Discover (prmtop, nctraj) pairs in a job's task directories.

generate_trajectory_viewer(project, experiment, job, verbose=True)
    Generate a multi-snapshot .mvsj file for trajectory visualisation.
    Each snapshot is titled "model N" (the short identifier derived
    from the task directory name).
"""

import json
import os

import molviewspec as mvs
from molviewspec.nodes import States, GlobalMetadata


# =====================================================================
#  File discovery
# =====================================================================


def find_trajectory_pairs(job_dir):
    """Find (prmtop, nctraj) pairs in a job's task directories.

    Scans task directories for ``prep/traj/`` folders containing
    a ``.prmtop`` topology file and a corresponding ``.nc`` NetCDF
    trajectory file.  Tasks without a ``prep/traj/`` directory are
    skipped gracefully.

    Parameters
    ----------
    job_dir : str
        Path to the job directory (e.g.
        ``CSS/MD/pmemd/out/CSS1_free_constrained/J1129505``).

    Returns
    -------
    list of (str, str, str)
        ``(task_name, prmtop_path, nctraj_path)`` tuples, one per task
        that has a valid trajectory pair.
    """
    pairs = []
    if not os.path.isdir(job_dir):
        return pairs
    tasks = sorted(os.listdir(job_dir))
    for task in tasks:
        if not task.startswith("task_"):
            continue
        traj_dir = os.path.join(job_dir, task, "prep", "traj")
        if not os.path.isdir(traj_dir):
            continue
        # Find all prmtop files
        prmtop_files = [f for f in os.listdir(traj_dir) if f.endswith(".prmtop")]
        if not prmtop_files:
            continue
        # Find all .nc files
        nc_files = {f for f in os.listdir(traj_dir) if f.endswith(".nc")}
        if not nc_files:
            continue
        # Pair by matching the prefix before "_stripped"
        for prmtop_f in prmtop_files:
            prefix = prmtop_f.replace("_stripped.prmtop", "")
            match = [f for f in nc_files if prefix in f]
            if match:
                pairs.append((
                    task,
                    os.path.join(traj_dir, prmtop_f),
                    os.path.join(traj_dir, match[0]),
                ))
                break
    return pairs


# =====================================================================
#  MVSJ generation
# =====================================================================


def generate_trajectory_viewer(project, experiment, job, verbose=True, base_url=None):
    """Generate a multi-snapshot .mvsj file from Amber MD trajectory data.

    Each task (replica) is shown as a snapshot with the topology
    (prmtop) and trajectory (nctraj) loaded for interactive viewing.
    The viewer shows DNA/RNA polymer as a cartoon and ligands in
    ball-and-stick representation.

    The MVSJ file is written to the job directory with the name::

        {experiment}_{job}_trajectory.mvsj

    Parameters
    ----------
    project : str
        Project directory name (e.g. ``"CSS"``, ``"PQ4"``).
    experiment : str
        Experiment name (e.g. ``"CSS1_free_constrained"``).
    job : str
        Job identifier (e.g. ``"J1129505"``).
    verbose : bool, default True
        If False, suppresses status prints.
    base_url : str or None, default None
        Root URL for remote file serving (e.g. a Zenodo record or web
        server).  When set, data file paths in the MVSJ are absolute
        URLs built as ``{base_url}/{project}/MD/pmemd/out/{experiment}/{job}/{rel_path}``.
        When ``None`` (default), paths are relative to the MVSJ file,
        suitable for local viewing.
    """
    job_dir = os.path.join(project, "MD", "pmemd", "out", experiment, job)
    pairs = find_trajectory_pairs(job_dir)

    if not pairs:
        if verbose:
            print(f"  No trajectory data found in {job_dir}")
        return

    # Relative path prefix from job directory to data files
    # When base_url is set, prepend the full remote path
    snapshots = []
    for task_name, prmtop_path, nctraj_path in pairs:
        task_dir = os.path.relpath(os.path.dirname(nctraj_path), job_dir)
        prmtop_name = os.path.basename(prmtop_path)
        nctraj_name = os.path.basename(nctraj_path)
        rel_prmtop = os.path.join(task_dir, prmtop_name).replace("\\", "/")
        rel_nctraj = os.path.join(task_dir, nctraj_name).replace("\\", "/")

        if base_url:
            remote_dir = f"{project}/MD/pmemd/out/{experiment}/{job}"
            prmtop_url = f"{base_url}/{remote_dir}/{rel_prmtop}"
            nctraj_url = f"{base_url}/{remote_dir}/{rel_nctraj}"
        else:
            prmtop_url = rel_prmtop
            nctraj_url = rel_nctraj

        b = mvs.create_builder()

        b.canvas(custom={
            "molstar_postprocessing": {
                "enable_outline": True,
            },
        })

        ds_topo = b.download(url=prmtop_url)
        ps_topo = ds_topo.parse(format="prmtop", ref="topo")

        ds_traj = b.download(url=nctraj_url)
        ds_traj.parse(format="nctraj", ref="traj")

        ms = ps_topo.model_structure(coordinates_ref="traj")

        poly = ms.component(selector="polymer")
        poly_rep = poly.representation(
            type="cartoon",
            custom={"molstar_representation_params": {"ignoreLight": True}},
        )
        poly_rep.color(color="#e1e9ef")

        lig = ms.component(selector="ligand")
        lig_rep = lig.representation(
            type="ball_and_stick",
            size_factor=0.5,
            custom={"molstar_representation_params": {"ignoreLight": True}},
        )
        lig_rep.color(color="#ff6666")

        snap = b.get_snapshot(
            title=f"model {task_name.split('_')[-1]}",
            linger_duration_ms=3000,
            transition_duration_ms=500,
        )
        snapshots.append(snap)

    if not snapshots:
        return

    mvsj_path = os.path.join(job_dir, f"{experiment}_{job}_trajectory.mvsj")
    states = States(
        metadata=GlobalMetadata(title=f"{experiment} trajectory"),
        snapshots=snapshots,
    )
    states_dict = states.model_dump(exclude_none=True)
    with open(mvsj_path, "w") as fh:
        json.dump(states_dict, fh, indent=2)
    if verbose:
        print(f"  Wrote trajectory MVSJ: {mvsj_path}")
