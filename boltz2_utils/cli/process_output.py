"""
CLI entry point for Boltz-2 output processing.

Usage:
    boltz-process-output --project CSS ^
        --model CSS1 --model CSS2 --model CSS3 ^
        --ligand-dict "{\"Corticosteron\": \"C0R\", \"free\": \"free\"}" ^
        --name-map "{\"C0R\": \"Corticosteron\"}"

For full documentation see README.md §5.
"""

import argparse
import json
from itertools import product

from boltz2_utils.process_boltz_results import process_single_experiment


def run(project_name, models, lgd_dict, name_map, suffixes, job_dirs="last", ligand_processing=True):
    for model, lgd, suffix in product(models, lgd_dict.values(), suffixes):
        print(f"Processing {model} with ligand {lgd} and suffix {suffix}")
        experiment_name = f"{model}_{lgd}{suffix}"

        process_single_experiment(
            project_name,
            experiment_name,
            LGD_DICT=lgd_dict,
            NAME_MAP=name_map,
            ligand_processing=ligand_processing,
            job_dirs=job_dirs,
        )


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Process Boltz-2 prediction outputs (PyMOL + RDKit)."
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
        "--model",
        action="append",
        metavar="NAME",
        help="Model/sequence name (repeatable).",
    )
    parser.add_argument(
        "--ligand-dict",
        help="JSON string mapping display names to CCD codes (e.g. "
             '{\\"Apt\\": \\"C0R\\", \\"free\\": \\"free\\"}).',
    )
    parser.add_argument(
        "--name-map",
        help="JSON string mapping CCD codes to display names (e.g. "
             '{\\"C0R\\": \\"Apt\\"}).',
    )
    parser.add_argument(
        "--suffix",
        nargs="*",
        default=["_constrained", ""],
        help="YAML/SLURM suffixes to process (default: _constrained '').",
    )
    parser.add_argument(
        "--job-dirs",
        default="last",
        help="Which J-directories to process: 'last', 'all', a JID, "
             "or a JSON list of JIDs (default: last).",
    )
    parser.add_argument(
        "--no-ligand-processing",
        action="store_true",
        help="Skip ligand extraction and chirality analysis.",
    )
    return parser.parse_args(argv)


def load_config(path):
    text = open(path, encoding="utf-8").read()
    if path.endswith((".yaml", ".yml")):
        import yaml
        return yaml.safe_load(text)
    return json.loads(text)


def resolve_job_dirs(val):
    if val is None:
        return None
    try:
        parsed = json.loads(val)
    except (json.JSONDecodeError, TypeError):
        return val
    if parsed == "all":
        return None
    return parsed


def main(argv=None):
    args = parse_args(argv)

    if args.config:
        cfg = load_config(args.config)
        job_dirs = resolve_job_dirs(cfg.get("job_dirs", "last"))
        run(
            project_name=cfg["project"],
            models=cfg["models"],
            lgd_dict=cfg.get("lgd_dict", cfg.get("LGD_DICT", {})),
            name_map=cfg.get("name_map", cfg.get("NAME_MAP", {})),
            suffixes=cfg.get("suffixes", ["_constrained", ""]),
            job_dirs=job_dirs,
            ligand_processing=cfg.get("ligand_processing", True),
        )
        return

    if not args.project:
        parser = argparse.ArgumentParser()
        parser.error("the following arguments are required: --project (or use --config)")

    lgd_dict = json.loads(args.ligand_dict) if args.ligand_dict else {}
    name_map = json.loads(args.name_map) if args.name_map else {}
    job_dirs = resolve_job_dirs(args.job_dirs)

    run(
        project_name=args.project,
        models=args.model,
        lgd_dict=lgd_dict,
        name_map=name_map,
        suffixes=args.suffix,
        job_dirs=job_dirs,
        ligand_processing=not args.no_ligand_processing,
    )


if __name__ == "__main__":
    main()
