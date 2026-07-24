"""
CLI entry point for Boltz-2 input file generation.

Usage:
    boltz-generate-input --config config.json
    boltz-generate-input --project my_project --seq Seq1 "ACG" --free

For full documentation see README.md §3.
"""

import argparse
import json
import sys
from pathlib import Path

from boltz2_utils.generate_boltz2_yaml import generate_boltz2_yamls
from boltz2_utils.generate_boltz2_slurm import generate_boltz2_slurm_files


def _normalize_sequences(sequences):
    seqs = {}
    overrides = {}
    for name, value in sequences.items():
        if isinstance(value, str):
            seqs[name] = value
        elif isinstance(value, dict):
            seqs[name] = value["sequence"]
            overrides[name] = {k: v for k, v in value.items() if k != "sequence"}
        else:
            raise TypeError(f"Sequence '{name}': expected str or dict, got {type(value).__name__}")
    return seqs, overrides


def run(
    project_name,
    sequences,
    ligands,
    molecule_type="dna",
    stem_length=8,
    max_distance=4.0,
    additional_pairs=None,
    custom_msa=False,
    output_dir=None,
    seq_stem_length=None,
):
    seq_stem_length = seq_stem_length or {}
    out = output_dir or project_name

    seqs, seq_overrides = _normalize_sequences(sequences)
    for name, sl in seq_stem_length.items():
        seq_overrides.setdefault(name, {})["stem_length"] = sl

    entries = []
    for seq_name, sequence in seqs.items():
        over = seq_overrides.get(seq_name, {})
        for lig_name, lig in ligands.items():
            entry = {
                "seq_name": seq_name,
                "sequence": sequence,
                "ligand_name": lig_name,
                "ccd": lig["ccd"],
                "smiles": lig["smiles"],
            }
            if over:
                entry.update(over)
            entry.setdefault("stem_length", stem_length)
            entries.append(entry)

    if custom_msa:
        for entry in entries:
            msa_path_full = Path(f"{project_name}/msa/{entry['seq_name']}.a3m")
            msa_path_yaml = f"msa/{entry['seq_name']}.a3m"
            entry["msa"] = str(msa_path_yaml) if msa_path_full.exists() else None

    paths = generate_boltz2_yamls(
        entries,
        output_dir=f"{out}/yaml",
        stem_length=stem_length,
        max_distance=max_distance,
        additional_pairs=additional_pairs,
        molecule_type=molecule_type,
    )

    generate_boltz2_slurm_files(
        paths,
        cluster_path=f"scratch/boltz/projects/{project_name}",
        output_dir=out,
    )


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Generate Boltz-2 YAML and SLURM input files."
    )
    parser.add_argument(
        "--config",
        help="Path to JSON config file. Overrides all other flags.",
    )
    parser.add_argument(
        "--project",
        help="Project directory name.",
    )
    parser.add_argument(
        "--molecule-type",
        choices=["dna", "rna", "protein"],
        default="dna",
        help="Molecule type (default: dna).",
    )
    parser.add_argument(
        "--seq",
        action="append",
        nargs=2,
        metavar=("NAME", "SEQUENCE"),
        help="Sequence name and nucleotide sequence (repeatable).",
    )
    parser.add_argument(
        "--ligand",
        action="append",
        nargs=2,
        metavar=("NAME", "CCD"),
        help="Ligand name and CCD code (repeatable).",
    )
    parser.add_argument(
        "--ligand-smiles",
        action="append",
        nargs=2,
        metavar=("NAME", "SMILES"),
        help="Ligand name and SMILES string (repeatable).",
    )
    parser.add_argument(
        "--free",
        action="store_true",
        help="Include a ligand-free ('free') entry.",
    )
    parser.add_argument(
        "--stem-length",
        type=int,
        default=8,
        help="Number of stem base-pair constraints (default: 8).",
    )
    parser.add_argument(
        "--max-distance",
        type=float,
        default=4.0,
        help="Max distance (Å) for stem contacts (default: 4.0).",
    )
    parser.add_argument(
        "--additional-pairs",
        help="JSON string of additional contact constraints. See README §11.",
    )
    parser.add_argument(
        "--custom-msa",
        action="store_true",
        help="Enable custom MSA (.a3m) files in <project>/msa/.",
    )
    parser.add_argument(
        "--seq-stem-length",
        action="append",
        nargs=2,
        metavar=("NAME", "STEM_LENGTH"),
        help="Per-sequence stem_length override (repeatable, e.g. --seq-stem-length Seq1 0).",
    )
    parser.add_argument(
        "--output-dir",
        help="Output directory for generated files (default: <project>).",
    )
    return parser.parse_args(argv)


def load_config(path):
    text = Path(path).read_text(encoding="utf-8")
    if path.endswith((".yaml", ".yml")):
        import yaml
        return yaml.safe_load(text)
    return json.loads(text)


def main(argv=None):
    args = parse_args(argv)

    if args.config:
        cfg = load_config(args.config)
        cfg_seq_sl = cfg.get("seq_stem_length", {})
        if isinstance(cfg_seq_sl, list):
            cfg_seq_sl = dict(cfg_seq_sl)
        run(
            project_name=cfg["project"],
            sequences=cfg["sequences"],
            ligands=cfg["ligands"],
            molecule_type=cfg.get("molecule_type", "dna"),
            stem_length=cfg.get("stem_length", 8),
            max_distance=cfg.get("max_distance", 4.0),
            additional_pairs=cfg.get("additional_pairs"),
            custom_msa=cfg.get("custom_msa", False),
            output_dir=cfg.get("output_dir"),
            seq_stem_length=cfg_seq_sl,
        )
        return

    if not args.project:
        parser = argparse.ArgumentParser()
        parser.error("the following arguments are required: --project (or use --config)")

    sequences = dict(args.seq) if args.seq else {}
    seq_stem_length = {}
    if args.seq_stem_length:
        for name, sl in args.seq_stem_length:
            seq_stem_length[name] = int(sl)
    ligands = {}
    if args.free:
        ligands["free"] = {"ccd": None, "smiles": None}
    if args.ligand:
        for name, ccd in args.ligand:
            ligands[name] = {"ccd": ccd, "smiles": None}
    if args.ligand_smiles:
        for name, smiles in args.ligand_smiles:
            ligands[name] = {"ccd": None, "smiles": smiles}

    run(
        project_name=args.project,
        sequences=sequences,
        ligands=ligands,
        molecule_type=args.molecule_type,
        stem_length=args.stem_length,
        max_distance=args.max_distance,
        additional_pairs=json.loads(args.additional_pairs) if args.additional_pairs else None,
        custom_msa=args.custom_msa,
        output_dir=args.output_dir,
        seq_stem_length=seq_stem_length,
    )


if __name__ == "__main__":
    main()
