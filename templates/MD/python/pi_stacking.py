#!/usr/bin/env python
"""Systematic ligand-DNA pi-pi stacking analysis.

Precompute step (run externally, like the cpptraj H-bond analysis):

    conda run -n boltz python <project>/MD/python/pi_stacking.py

Copy this script into your project's MD/python/ directory and edit the
Configuration block below (PROJECT, EXPERIMENT, LIGAND_RESNAME, LIGAND_SMILES).
For every aromatic ring of the ligand and every aromatic base ring of the DNA,
it evaluates the standard geometric pi-stacking criteria for every analysed
trajectory frame (stride 5, matching the H-bond analysis), without any prior
knowledge of which interactions are formed. Aromatic ligand rings are
enumerated automatically from the ligand SMILES (in the protonation state of
the MD species) via RDKit, mapped onto the topology atom names; base rings use
the standard DNA base ring atom sets.

Criteria (tunable via command line):
  - parallel/displaced (sandwich): centroid distance <= dist_cut (5.5 A) and
    acute interplanar angle <= angle_cut (30 deg) and lateral offset
    <= offset_cut (2.0 A)
  - T-shaped (edge-to-face):      centroid distance <= dist_cut (5.5 A) and
    acute interplanar angle in the 60-90 deg band

Outputs (written to <project>/MD/pi_stacking/):
  - {experiment}_{model}.csv  one row per (frame, ring pair) within a buffer
    distance, with columns: Frame, pair, lig_ring, base_residue, present,
    mode, d, angle, offset.
  - ring_defs.json  aromatic ring atom definitions (atom names as used in the
    Amber topology/PDB) for the ligand and the DNA bases, read by md_utils.R
    for the in-situ evaluation of the final minimised structures.

NOTE: re-run this script after feeding any new MD data into the report.
"""
import argparse
import json
import os
import sys

import numpy as np
import MDAnalysis as mda
from rdkit import Chem

# ---------------------------------------------------------------------------
# Configuration — EDIT FOR YOUR PROJECT
# ---------------------------------------------------------------------------
# DNA base ring atom sets (standard Amber/PDB names). Unchanged between
# projects unless you use a non-canonical base.
BASE_RINGS = {
    "DG": [["N1", "C2", "N3", "C4", "C5", "C6"], ["C4", "C5", "N7", "C8", "N9"]],
    "DA": [["N1", "C2", "N3", "C4", "C5", "C6"], ["C4", "C5", "N7", "C8", "N9"]],
    "DC": [["N1", "C2", "N3", "C4", "C5", "C6"]],
    "DT": [["N1", "C2", "N3", "C4", "C5", "C6"]],
}
PROJECT = "PLACEHOLDER_PROJECT"        # project directory name (e.g. "PQ4")
EXPERIMENT = "PLACEHOLDER_EXPERIMENT"  # e.g. "PQ4_PQ_constrained"
LIGAND_RESNAME = "PLACEHOLDER_LIGAND"  # ligand residue name in the topology/PDB (e.g. "PQ")
LIGAND_BASE = f"{PROJECT}/MD/pmemd/out"
OUT_DIR = f"{PROJECT}/MD/pi_stacking"
BUFFER = 6.5  # only store ring pairs within this centroid distance (A)


# ---------------------------------------------------------------------------
# Ligand SMILES
# ---------------------------------------------------------------------------
# SMILES of the ligand in its MD protonation state. The heavy-atom graph must
# be isomorphic to the topology graph; hydrogens are ignored in the mapping.
# Example (piperaquine, neutral bis-quinoline):
#   "C(CN1CCN(CC1)c1ccnc2cc(Cl)ccc12)CN1CCN(CC1)c1ccnc2cc(Cl)ccc12"
LIGAND_SMILES = "PLACEHOLDER_SMILES"


def element_of_name(name):
    """Guess an element symbol from an Amber/GAFF atom name."""
    base = name.rstrip("0123456789'\"+-")
    if base in ("CL", "Cl", "cl"):
        return "Cl"
    return base[0].upper()


def build_topology_mol(atom_names, bonds):
    """RDKit molecule from the topology graph (element + connectivity)."""
    el = [element_of_name(n) for n in atom_names]
    m = Chem.RWMol()
    for i, name in enumerate(atom_names):
        m.AddAtom(Chem.Atom(el[i]))
    for i, j in bonds:
        m.AddBond(i, j, Chem.BondType.UNSPECIFIED)
    return m.GetMol()


def ligand_rings_from_topology(prmtop, ligand_resname, smiles):
    """Map aromatic rings from SMILES onto topology atom names.

    Uses the heavy-atom graph only: the neutral SMILES molecule and the
    topology heavy-atom connectivity must be isomorphic (the protonation
    state of the MD species only affects hydrogen atoms).

    Returns a dict {ring_name: [atom names in topology]}.
    """
    u = mda.Universe(prmtop)
    lig = u.select_atoms(f"resname {ligand_resname}")
    if lig.n_atoms == 0:
        raise RuntimeError(f"no {ligand_resname} residue in {prmtop}")

    heavy = lig[lig.elements != "H"]
    lig_pos = {idx: pos for pos, idx in enumerate(heavy.indices)}
    bonds = []
    for bond in u.bonds:
        i, j = bond.indices
        if i in lig_pos and j in lig_pos:
            bonds.append((lig_pos[i], lig_pos[j]))

    m1 = build_topology_mol(heavy.names, bonds)
    m2 = Chem.MolFromSmiles(smiles)
    if m2 is None:
        raise RuntimeError(f"could not parse SMILES: {smiles}")

    # match[m2_idx] = m1 index of the matching topology atom
    match = m1.GetSubstructMatch(m2)
    if not match or len(match) != m2.GetNumAtoms():
        raise RuntimeError("SMILES <-> topology graph mismatch for the ligand")

    rings = [
        list(r) for r in m2.GetRingInfo().AtomRings()
        if all(m2.GetAtomWithIdx(i).GetIsAromatic() for i in r)
    ]
    ring_names = {}
    for i, ring in enumerate(rings):
        topo_idx = [match[j] for j in ring]
        ring_names[f"R{i}"] = [heavy.names[k] for k in topo_idx]
    return ring_names


def base_ring_groups(u, base_rings_def):
    """Return list of (resid, resname, ring_name, [atom names]) for the DNA."""
    groups = []
    for res in u.residues:
        rn = res.resname
        if not rn.startswith(("DG", "DA", "DC", "DT")):
            continue
        key = rn[:2]
        if key not in base_rings_def:
            continue
        names = set(res.atoms.names)
        for k, ring in enumerate(base_rings_def[key]):
            if any(n not in names for n in ring):
                continue
            groups.append((res.resid, rn, f"base{k}", list(ring)))
    return groups


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------
def ring_centroid(coords):
    return coords.mean(axis=0)


def ring_normal(coords):
    """Unit normal to the best-fit plane (smallest singular vector)."""
    c = coords - coords.mean(axis=0)
    _, _, vt = np.linalg.svd(c)
    return vt[2]


def acute_angle(n1, n2):
    """Acute angle (deg) between two unit vectors."""
    cosv = abs(float(np.dot(n1, n2)))
    return float(np.degrees(np.arccos(min(1.0, cosv))))


# ---------------------------------------------------------------------------
# Per-trajectory analysis
# ---------------------------------------------------------------------------
def analyse_trajectory(prmtop, nc, lig_rings, base_groups, stride,
                       dist_cut, angle_cut, offset_cut):
    """Return a list of row dicts for one trajectory."""
    u = mda.Universe(prmtop, nc)
    lig = u.select_atoms(f"resname {LIGAND_RESNAME}")
    if lig.n_atoms == 0:
        return []

    # Precompute per-ring atom index arrays
    lig_ring_idx = {}
    for ring_name, names in lig_rings.items():
        sel = lig.select_atoms(f"name {' '.join(names)}")
        lig_ring_idx[ring_name] = sel.indices
    base_ring_idx = {}
    for resid, resname, ring_name, names in base_groups:
        key = (resid, resname, ring_name)
        sel = u.select_atoms(f"resid {resid} and name {' '.join(names)}")
        base_ring_idx[key] = sel.indices

    rows = []
    for frame_idx, ts in enumerate(u.trajectory[::stride], start=1):
        frame = frame_idx  # analysed-frame index (1-based), as in the cpptraj CSVs
        pos = u.atoms.positions

        lc, ln = {}, {}
        for rn, idx in lig_ring_idx.items():
            c = pos[idx]
            lc[rn] = ring_centroid(c)
            ln[rn] = ring_normal(c)
        bc, bn = {}, {}
        for key, idx in base_ring_idx.items():
            c = pos[idx]
            bc[key] = ring_centroid(c)
            bn[key] = ring_normal(c)

        for rn, c1 in lc.items():
            n1 = ln[rn]
            for key, c2 in bc.items():
                dvec = c2 - c1
                d = float(np.linalg.norm(dvec))
                if d > BUFFER:
                    continue
                resid, resname, bname = key
                ang = acute_angle(n1, bn[key])
                h = abs(float(np.dot(dvec, n1)))
                off = float(np.sqrt(max(0.0, d * d - h * h)))
                if d <= dist_cut and ang <= angle_cut and off <= offset_cut:
                    mode, present = "parallel", 1
                elif d <= dist_cut and 60.0 <= ang <= 90.0:
                    mode, present = "T-shaped", 1
                else:
                    mode, present = "", 0
                rows.append(dict(
                    frame=frame, pair=f"{rn}-{resname}{resid}",
                    lig_ring=rn, base_residue=f"{resname}{resid}",
                    present=present, mode=mode, d=round(d, 3),
                    angle=round(ang, 2), offset=round(off, 3)))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stride", type=int, default=5)
    ap.add_argument("--dist-cut", type=float, default=5.5)
    ap.add_argument("--angle-cut", type=float, default=30.0)
    ap.add_argument("--offset-cut", type=float, default=2.0)
    args = ap.parse_args()

    exp_dir = os.path.join(LIGAND_BASE, EXPERIMENT)
    job_dirs = [d for d in os.listdir(exp_dir)
                if d.startswith("J") and os.path.isdir(os.path.join(exp_dir, d))]
    if not job_dirs:
        sys.exit(f"no job dirs under {exp_dir}")
    last_job = sorted(job_dirs, key=lambda x: int(x[1:]))[-1]

    # Ring definitions (from topology, which defines the atom names)
    prmtop_ref = os.path.join(exp_dir, last_job, "task_0", "prep", "traj",
                              f"step10_0_{last_job}_stripped.prmtop")
    lig_rings = ligand_rings_from_topology(prmtop_ref, LIGAND_RESNAME,
                                           LIGAND_SMILES)
    print("Ligand aromatic rings (topology atom names):")
    for rn, names in lig_rings.items():
        print("  ", rn, names)

    os.makedirs(OUT_DIR, exist_ok=True)
    # ring_defs.json for the R in-situ evaluation (final minimised structures)
    base_defs = {k: [list(r) for r in v] for k, v in BASE_RINGS.items()}
    with open(os.path.join(OUT_DIR, "ring_defs.json"), "w") as fh:
        json.dump({"ligand": {LIGAND_RESNAME: lig_rings},
                   "bases": base_defs}, fh, indent=2)

    # Per-model base rings come from a representative topology (same sequence)
    u_ref = mda.Universe(prmtop_ref)
    base_groups = base_ring_groups(u_ref, BASE_RINGS)
    print(f"{len(base_groups)} DNA base rings")

    task_dirs = sorted(
        [d for d in os.listdir(os.path.join(exp_dir, last_job))
         if d.startswith("task_")],
        key=lambda x: int(x.split("_")[1]))
    for task in task_dirs:
        t_idx = task.split("_")[1]
        base = os.path.join(exp_dir, last_job, task, "prep", "traj")
        prmtop = os.path.join(base, f"step10_{t_idx}_{last_job}_stripped.prmtop")
        nc = os.path.join(base, f"step10_{t_idx}_{last_job}_stripped.nc")
        if not (os.path.exists(prmtop) and os.path.exists(nc)):
            print(f"skip {task}: missing prmtop/nc")
            continue
        print(f"analysing {EXPERIMENT}/{task} ...")
        rows = analyse_trajectory(prmtop, nc, lig_rings, base_groups,
                                  args.stride, args.dist_cut,
                                  args.angle_cut, args.offset_cut)
        out = os.path.join(OUT_DIR, f"{EXPERIMENT}_{task}.csv")
        if rows:
            cols = ["frame", "pair", "lig_ring", "base_residue", "present",
                    "mode", "d", "angle", "offset"]
            header = "Frame,pair,lig_ring,base_residue,present,mode,d,angle,offset"
            with open(out, "w") as fh:
                fh.write(header + "\n")
                for r in sorted(rows, key=lambda r: (r["frame"], r["pair"])):
                    fh.write(",".join(str(r[c]) for c in cols) + "\n")
            print(f"  wrote {out} ({len(rows)} rows)")
        else:
            print(f"  {task}: no candidate pairs (empty)")

    print("done")


if __name__ == "__main__":
    main()
