from boltz2_utils.generate_boltz2_yaml import generate_boltz2_yamls
from boltz2_utils.generate_boltz2_slurm import generate_boltz2_slurm_files
from pathlib import Path

SEQUENCES = {
    "CSS1": "GGGACGACGCCCGCATGTTCCATGGATAGTCTTGACTAGTCGTCCC",
    "CSS2": "GGGACGACTAGCGTATGCGCCAGAAGTATACGAGGATAGTCGTCCC",
    "CSS3": "GGGACGACGCCAGAAGTTTACGAGGATATGGTAACATAGTCGTCCC",
}

LIGANDS = {
    "C0R": {"ccd": "C0R", "smiles": None},
    "HCY": {"ccd": "HCY", "smiles": None},
    "11DC": {
        "ccd": None,
        "smiles": (
            r"O=C4\C=C2/[C@]([C@H]1CC[C@@]3([C@@](O)(C(=O)CO)CC[C@H]3"
            r"[C@@H]1CC2)C)(C)CC4"
        ),
    },
}

# Build the flat list of entries
entries = [
    {
        "seq_name": seq_name,
        "sequence": sequence,
        "ligand_name": lig_name,
        "ccd": lig["ccd"],
        "smiles": lig["smiles"],
    }
    for seq_name, sequence in SEQUENCES.items()
    for lig_name, lig in LIGANDS.items()
]

# Toggle use of custom msa
custom_msa = False

# NOTE: if custom_msa is True, the msa files must be named as seq_name.a3m and placed in CSS/msa/ directory. If the file does not exist, the msa field will be set to None in the YAML files.

if custom_msa:
    # To entries, add msa: that is "msa/{seq_name}.a3m" if the file exists, else None
    for entry in entries:
        msa_path_full = Path(f"CSS/msa/{entry['seq_name']}.a3m")
        msa_path_yaml = f"msa/{entry['seq_name']}.a3m"
        entry["msa"] = str(msa_path_yaml) if msa_path_full.exists() else None

# Generate the YAML files
paths = generate_boltz2_yamls(entries, output_dir="CSS/yaml", stem_length=8, max_distance=4.0)

# Generate the SLURM submission scripts
# paths contains the last generated YAML files
# so custom_msa toggle will be reflected in the SLURM scripts as well
slurm_paths = generate_boltz2_slurm_files(
    paths,
    cluster_path="scratch/boltz/projects/CSS",
    output_dir="CSS",
)


# USAGE OF FILES ON CALI3 USER ROOT----
# Run all the generated SLURM scripts with
# for f in scratch/boltz/projects/CSS/*.slurm; do sbatch "$f"; sleep 2; done

# Run only those with custom msa
# for f in scratch/boltz/projects/CSS/*_custom_msa.slurm; do sbatch "$f"; sleep 2; done