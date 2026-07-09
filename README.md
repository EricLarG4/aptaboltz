# aptaBoltz — Nucleic-Acid / Aptamer Analysis Pipeline for Boltz-2

A modular pipeline for generating Boltz-2 input files, processing
prediction outputs, and visualising results.  Designed primarily for
**DNA and RNA aptamers** with small-molecule ligands.  The pipeline is
split into two layers:

1. **Shared utility packages** (imported / sourced by all projects)
   - `boltz2_utils/`  — Python package (YAML generation, SLURM script
     generation, output processing with PyMOL + RDKit)
   - `boltz_R_utils/` — R scripts (confidence table, PAE/PDE heatmaps,
     pLDDT per-residue plots)

2. **Project-specific template files** (copy from `templates/` to each project and edit)
   - `templates/input_file_generator.py` — define sequences, ligands, constraints
   - `templates/output_file_processing.py` — invoke the Python post-processing
   - `templates/process.R` — invoke the R visualisation

   The `templates/` directory contains these three files with `PLACEHOLDER_` values
   that you replace with your project's actual names and sequences.

> **Note on molecule types:** This pipeline works natively with **DNA** and **RNA**
> (nucleic-acid aptamers).  Basic **protein** support is available: the YAML
> generator accepts `molecule_type="protein"`, and PyMOL styling colours protein
> chains in a uniform muted colour.  However, the stem-contact logic (Watson-Crick
> base pairing) is skipped for proteins, and the R pLDDT visualiser labels residues
> using single-letter codes read directly from the YAML sequence.

---

## Table of Contents

1. [Prerequisites and Installation](#1-prerequisites-and-installation)
2. [Setting Up a New Project](#2-setting-up-a-new-project)
3. [Generating Input Files (YAML + SLURM)](#3-generating-input-files-yaml--slurm)
4. [Running Boltz-2 on a Cluster](#4-running-boltz-2-on-a-cluster)
5. [Post-Processing: Python (PyMOL + RDKit)](#5-post-processing-python-pymol--rdkit)
6. [Post-Processing: R Visualisation](#6-post-processing-r-visualisation)
7. [Molecule Type Support](#7-molecule-type-support)
8. [MD Force-Field Parameterisation Pipeline](#8-md-force-field-parameterisation-pipeline)
9. [File Reference](#9-file-reference)
10. [Appendix: Contact Constraints Format](#10-appendix-contact-constraints-format)

---

## 1. Prerequisites and Installation

### Python package (boltz2_utils)

Install the shared Python utilities in editable mode:

```bash
pip install -e ./boltz2_utils
```

This makes the following modules importable:

- `boltz2_utils.process_boltz_results` — prediction output processing
- `boltz2_utils.generate_boltz2_yaml` — YAML input file generation
- `boltz2_utils.generate_boltz2_slurm` — SLURM job script generation

### R environment

Run the one-time setup script to install required R packages and create
the `boltz` conda environment (defined in `environment.yml`):

```r
source("boltz_R_utils/install_packages.R")
```

This script is idempotent — safe to re-run.

**Required R packages:** `data.table`, `ggplot2`, `DT`, `khroma`,
`jsonlite`, `yaml`, `bslib`, `reticulate`

**Required Python (conda) packages:** defined in `environment.yml`
(automatically installed into the `boltz` conda environment).

---

## 2. Setting Up a New Project

For each new project (e.g. `my_project`):

1. Create a project directory:
   ```bash
   mkdir my_project
   mkdir my_project/yaml
   mkdir my_project/msa      # optional — only if using custom MSA
   ```

2. Copy the three template files from `templates/` into the project directory:
   ```bash
   cp templates/input_file_generator.py my_project/
   cp templates/output_file_processing.py my_project/
   cp templates/process.R my_project/
   ```

   Each file contains `PLACEHOLDER_` values that you replace with your
   project's actual names and sequences (see Steps 3–6 below).

3. The shared utilities (`boltz2_utils/`, `boltz_R_utils/`) remain in
   the repository root — do **not** copy them.  They are imported /
   sourced automatically.

---

## 3. Generating Input Files (YAML + SLURM)

Edit the project-local copy of `input_file_generator.py`.

### Step 1 — Set project name, molecule type, and sequences

```python
project = "my_project"
molecule_type = "dna"       # "dna", "rna", or "protein"
```

Define sequences in the `SEQUENCES` dict:

```python
SEQUENCES = {
    "Seq1": "GGGACGACGCCCGCATGTTCCATGGATAGTCTTGACTAGTCGTCCC",
    "Seq2": "GGGACGACTAGCGTATGCGCCAGAAGTATACGAGGATAGTCGTCCC",
}
```

Adjust stem constraints if needed (only applies to DNA and RNA):

```python
stem_length = 8
max_distance = 4.0
```

### Step 2 — (Optional) Define contact constraints

Set `additional_pairs` to add extra residue-pair contact constraints.
Three formats are accepted.  See [Appendix 10](#10-appendix-contact-constraints-format)
for full details and examples.

Set to `None` to omit all additional constraints.

### Step 3 — Define ligands

Each ligand is specified either by CCD code (for standard ligands from
the PDB) or by SMILES string (for novel / non-standard ligands).

```python
LIGANDS = {
    "free": {"ccd": None, "smiles": None},               # no ligand
    "L01":  {"ccd": "C0R", "smiles": None},              # CCD-based
    "L02": {                                              # SMILES-based
        "ccd": None,
        "smiles": "C1=CC=C2C(=C1)C(=O)C3=C(C2=O)C4=CC=CC=C4C3=O",
    },
}
```

### Step 4 — Build the flat list of (sequence × ligand) entries

The template includes the nested loop automatically.  No user action
required.

### Step 5 — Toggle custom MSA usage

```python
custom_msa = False     # set to True if using custom MSA files
```

If `True`, place `.a3m` MSA files in `{project}/msa/` and name them
`{seq_name}.a3m` (e.g. `Seq1.a3m`).

### Step 6—7 — Generate files

Run the script:

```bash
python my_project/input_file_generator.py
```

This produces:

- YAML files in `{project}/yaml/` — one per (sequence × ligand × variant)
- SLURM submission scripts in `{project}/` — one per YAML file

Naming conventions:

| Variant                     | Example filename                       |
|-----------------------------|----------------------------------------|
| Standard                    | `Seq1_L01.yaml` / `.slurm`            |
| With custom MSA             | `Seq1_L01_custom_msa.yaml` / `.slurm` |
| With contact constraints    | `Seq1_L01_constrained.yaml` / `.slurm`|
| MSA + constraints           | `Seq1_L01_custom_msa_constrained.yaml` / `.slurm` |

---

## 4. Running Boltz-2 on a Cluster

1. Transfer the following to the cluster:
   ```
   my_project/              → scratch/boltz/projects/
   my_project/yaml/         → scratch/boltz/projects/my_project/yaml/
   my_project/msa/          → scratch/boltz/projects/my_project/msa/   (if used)
   my_project/*.slurm       → scratch/boltz/projects/my_project/
   ```

2. Submit jobs.  Useful patterns:

   ```bash
   # Submit ALL jobs
   for f in scratch/boltz/projects/my_project/*.slurm; do
       sbatch "$f"; sleep 2
   done

   # Submit only custom_msa jobs
   for f in scratch/boltz/projects/my_project/*_custom_msa.slurm; do
       sbatch "$f"; sleep 2
   done

   # Submit only constrained jobs
   for f in scratch/boltz/projects/my_project/*_constrained.slurm; do
       sbatch "$f"; sleep 2
   done

   # Submit jobs that are NOT custom_msa
   for f in scratch/boltz/projects/my_project/*.slurm; do
       [[ "$f" != *_custom_msa.slurm ]] && sbatch "$f" && sleep 2
   done

   # Submit jobs that are NOT constrained
   for f in scratch/boltz/projects/my_project/*.slurm; do
       [[ "$f" != *_constrained.slurm ]] && sbatch "$f" && sleep 2
   done
   ```

---

## 5. Post-Processing: Python (PyMOL + RDKit)

Edit the project-local copy of `output_file_processing.py`.

### Edit per-project configuration

```python
LGD_DICT = {"Aptamer1": "L01", "Aptamer2": "L02"}
NAME_MAP = {"L01": "Aptamer1", "L02": "Aptamer2", ...}

project = "my_project"
models = ["Seq1", "Seq2"]
suffixes = ["_constrained", ""]
```

The `LGD_DICT` maps aptamer display names to their CCD ligand codes.
The `NAME_MAP` reverses this mapping.  Both control ligand extraction
and chirality visualisation.

### Run

```bash
python my_project/output_file_processing.py
```

### What this does

For each (model × ligand × suffix × job-directory) combination:

1. **Loads CIF predictions** into PyMOL (all models from the job).
2. **Aligns and saves** a multi-object CIF and a `.pse` session file.
3. **Colours** nucleic-acid chains by base (DNA/RNA), ligand by element,
   protein chains in a muted uniform colour, constrained residues by
   component (pastel palette), and all residues by pLDDT (red→blue).
4. **Extracts ligands** as SDF files from the predicted structures and
   fetches the corresponding reference ligand from the PDB.
5. **Analyses chirality** — loads SDFs into RDKit, detects chiral
   centres, and draws an annotated grid image (PNG) with chiral centres
   highlighted.
6. **Aggregates confidence / PAE / PDE / pLDDT** metrics from the
   prediction JSON/NPZ files and saves them as CSV files (one per
   experiment).

### Controlling which job directories are processed

The `job_dirs` parameter of `process_single_experiment()` accepts:

| Value                   | Behaviour                                    |
|-------------------------|----------------------------------------------|
| `None`                  | Process **all** J-directories                |
| `"last"`                | Process only the highest-numbered directory  |
| `"J1121077"`            | Process only that specific directory         |
| `["J1121077", "J1123593"]` | Process only the listed directories      |

---

## 6. Post-Processing: R Visualisation

Edit the project-local copy of `process.R`.

### Minimal setup

```r
source('boltz_R_utils/install_packages.R')
source('boltz_R_utils/processor.R')

project <- 'my_project'

# 1. Interactive confidence table
table_confidence(process_confidence(project))

# 2. PAE / PDE heatmaps
pae_dt <- process_pxe(project, type = "pae")
pde_dt <- process_pxe(project, type = "pde")

lapply(names(pae_dt), function(n) {
    plot_pxe(pae_dt[n], type = "pae", ligand_number = 47)
})
lapply(names(pde_dt), function(n) {
    plot_pxe(pde_dt[n], type = "pde", ligand_number = 47)
})

# 3. pLDDT per-residue plots
plot_plddt(project)
```

Adjust `ligand_number` to the residue index at which the ligand region
begins in your system (i.e. the number of nucleic-acid / protein residues).

### What each function does

| Function                | Input                                 | Output                                                           |
|-------------------------|---------------------------------------|------------------------------------------------------------------|
| `process_confidence()`  | Reads all `*_confidence.csv` files    | A `data.table` with one row per model, columns for each metric  |
| `table_confidence()`    | The returned data.table               | Interactive HTML table saved to `{project}/confidence_table.html`|
| `process_pxe()`         | Reads all `*_pae.csv` or `*_pde.csv`  | Named list of `data.table`s (one per experiment)                 |
| `plot_pxe()`            | One element from `process_pxe()`      | Faceted heatmap PNG per region (DNA / Ligand)                   |
| `plot_plddt()`          | Scans project for `*_plddt.csv` files | Faceted line-plot PNG per experiment                             |

---

## 7. Molecule Type Support

The `molecule_type` parameter (`input_file_generator.py`) controls how
the pipeline handles different biopolymers:

| Type      | YAML key | Stem contacts | PyMOL styling                    | pLDDT labels      |
|-----------|----------|---------------|----------------------------------|--------------------|
| `"dna"`   | `dna`    | ✓             | Base-coloured (A, G, T, C, U)    | Single-letter DNA  |
| `"rna"`   | `rna`    | ✓             | Base-coloured (A, G, U, C)       | Single-letter RNA  |
| `"protein"` | `protein` | ✗ (no base pairing) | Uniform muted colour (`warmpink`) | 1-letter AA codes  |

- **DNA and RNA** are fully supported with Watson-Crick stem contacts,
  base-identity colouring, and residue labelling.
- **Protein** support is basic: the YAML generator and PyMOL styling
  work, but stem contacts are omitted.  Protein-specific features
  (secondary-structure colouring, domain parsing) are not implemented.

The ligand-extraction and chirality-analysis components are
molecule-type agnostic and work for any system.

---


## 8. MD Force-Field Parameterisation Pipeline

In addition to the Boltz-2 prediction pipeline, the `templates/` directory
contains four scripts for generating AMBER-compatible force-field parameters
from quantum-mechanical (QM) calculations.  This is a separate workflow used
to prepare **ligand parameters** for classical molecular-dynamics simulations
that complement the Boltz-2 structure predictions.

### Pipeline Overview

```
QM optimisation          ──→  orca_steps_wsl.sh         (WSL)
   │                              │
   │                     (copy .molden files to Windows)
   │                              │
   │                     multiwfn_steps_windows.ps1     (Windows)
   │                              │
   │                     (copy .chg file to WSL)
   │                              │
   ├──→  ligand_prep.sh          (WSL)
   │         │
   │         ├──  antechamber/   (GAFF2 types + BCC charges)
   │         ├──  ff/            (prepin, frcmod, resp.pdb)
   │         └──  (auxiliary files)
   │
   └──→  leap.sh                 (WSL)
              │
              └──  leap/CSS1_<LIGAND>_constrained/
                     (prmtop, rst7, pdb per model)
```

### Step-by-Step Workflow

#### Step 1 — ORCA single-point calculations (`orca_steps_wsl.sh`)

Runs two consecutive single-point energy calculations with ORCA on an
optimised ligand geometry (XYZ format):

- **Gas phase** (B3LYP-D3/def2-TZVP with RIJCOSX)
- **Solvent phase** (same level, CPCM-SMD continuum solvation)

Both calculations produce Molden-format files with an Nval valence-electron
header (required for Multiwfn RESP fitting with def2 ECP basis sets).

**Deliverables:** `SP_gas.molden`, `SP_solv.molden`

#### Step 2 — RESP2 charge fitting (`multiwfn_steps_windows.ps1`)

Runs on Windows.  Reads the gas- and solvent-phase Molden files (copied from
WSL) and invokes Multiwfn to compute RESP charges for each phase.  It then
combines the two charge sets into RESP2 charges using the mixing parameter
delta (default: 0.5, i.e. RESP2-0.5).

**Deliverable:** `<ligand>.chg` — RESP2 charges for the ligand

#### Step 3 — Ligand parameterisation (`ligand_prep.sh`)

Takes the optimised `.mol2` and RESP2 `.chg` files from the `QM/` directory
and runs the full AmberTools pipeline:

1. **Antechamber** — assigns GAFF2 atom types and BCC charges.
2. **Charge substitution** — replaces BCC charges with RESP2 charges.
3. **Prepgen** — generates an Amber prepin library file.
4. **Parmchk2** — checks for missing force-field parameters (frcmod).

**Deliverables:**
- `ff/<ligand>_resp.prepin` — Amber prepin library
- `ff/<ligand>_resp.frcmod` — Force-field modifications
- `ff/<ligand>_resp.pdb` — PDB-like file for tLeap complex preparation

#### Step 4 — System solvation and ionisation (`leap.sh`)

Creates a combined tLEAP input script that:
- Loads force fields (GAFF2, DNA.OL21, OPC water)
- Loads the ligand prepin and frcmod
- Loads each PDB model from `pdb_for_md/CSS1_<LIGAND>_constrained/`
- Solvates in an OPC water box (14.0 Å buffer)
- Adds Na⁺/Mg²⁺/Cl⁻ ions with automatic Cl⁻ distribution (MgCl₂ + NaCl)

Runs tLEAP on the combined input, then converts each rst7 to PDB with ambpdb.

**Deliverables (per model):** `.prmtop`, `.rst7`, `.pdb`

### Dependencies

| Script | Platform | Dependencies |
|--------|----------|-------------|
| `orca_steps_wsl.sh` | WSL/Linux | ORCA (≥ 6.x), orca_2mkl |
| `multiwfn_steps_windows.ps1` | Windows | Multiwfn, PowerShell |
| `ligand_prep.sh` | WSL/Linux | AmberTools (antechamber, prepgen, parmchk2) |
| `leap.sh` | WSL/Linux | AmberTools (tleap, ambpdb) |

### Troubleshooting

| Issue | Likely Cause | Solution |
|-------|-------------|----------|
| `SP_gas.molden` not found | WSL→Windows copy not performed | Copy Molden files from WSL to `<MoldenDir>` before running the PS1 script |
| Multiwfn fails silently | Batch wrapper issue | Check host output for errors; verify Multiwfn path |
| `mol2` file not found in `QM/` | Wrong filename or directory | Name the file `<MOL>.mol2` and place in `QM/` |
| RESP2 charges look wrong | Atom ordering mismatch | Verify identical atom ordering between `.mol2` and `.chg` files |
| tLEAP reports missing parameters | frcmod has ATTN warnings | Review `ff/<ligand>_resp.frcmod` and adjust |
| OPC water model unrecognised | AmberTools version | Requires AmberTools ≥ 18 for `leaprc.water.opc` |

---

## 9. File Reference

### Shared utilities (imported, do not edit for each project)

| File | Purpose |
|------|---------|
| `boltz2_utils/__init__.py` | Package initialisation |
| `boltz2_utils/base.py` | PyMOL styling script (nucleic-acid base colours, protein detection, solvent/ion hiding) |
| `boltz2_utils/process_boltz_results.py` | PyMOL alignment + colouring, ligand extraction, chirality grid, CSV aggregation (confidence, PAE, PDE, pLDDT) |
| `boltz2_utils/process_boltz_results_nolgd.py` | Simplified version without ligand processing (legacy) |
| `boltz2_utils/generate_boltz2_yaml.py` | YAML input builder: stem contacts, additional constraints, MSA integration, molecule-type support |
| `boltz2_utils/generate_boltz2_slurm.py` | SLURM script generator (one per YAML file) |
| `boltz_R_utils/install_packages.R` | One-time R + conda setup |
| `boltz_R_utils/processor.R` | Confidence/PAE/PDE/pLDDT reading, melting, plotting |

### Project-specific templates (copy and edit per project)

| File | Edits required |
|------|----------------|
| `templates/input_file_generator.py` | `project`, `molecule_type`, `SEQUENCES`, `LIGANDS`, `additional_pairs`, `custom_msa` |
| `templates/output_file_processing.py` | `project`, `models`, `LGD_DICT`, `NAME_MAP`, `suffixes` |
| `templates/process.R` | `project`, `ligand_number` |
| `templates/orca_steps_wsl.sh` | ORCA paths, molecule name, charge, multiplicity, solvent |
| `templates/ligand_prep.sh` | `MOL`, `NC`, `S` (molecule name, net charge, verbosity) |
| `templates/leap.sh` | Ligand name, model range, ion counts |
| `templates/multiwfn_steps_windows.ps1` | Multiwfn path, Molden directory, output directory |

---

## 10. Appendix: Contact Constraints Format

The `additional_pairs` parameter in `input_file_generator.py` accepts
three formats, from simplest to most powerful.

### Format A: Flat list (backward-compatible)

A simple list of `(res_i, res_j)` 1-based residue indices.  Both
residues are assigned to the sequence chain (default `"A"`).  The same
set of pairs is applied to **every** (sequence × ligand) entry.

```python
additional_pairs = [(10, 30)]                       # single pair
additional_pairs = [(1, 44), (2, 43), (3, 42)]     # manual stem pairs
additional_pairs = None                             # no extra pairs
```

### Format B: Per-sequence dictionary

A dict mapping each sequence name to its own list of pairs.  Use the
same keys as `SEQUENCES`.  Useful when different sequences need
different constraints, regardless of the ligand.

Each pair can be:

- A 2-tuple `(res_i, res_j)` — both on the sequence chain
- A 4-tuple `(chain_i, res_i, chain_j, res_j)` — explicit chains
  (for cross-chain contacts, e.g. nucleic-acid to ligand)

```python
additional_pairs = {
    "Seq1": [("A", 5, "A", 36)],          # chain–chain, explicit chains
    "Seq2": [("A", 3, "B", 10)],          # cross-chain: DNA → ligand
}
```

### Format C: Per-sequence + per-ligand dictionary (with wildcard)

A two-level dict: outer key = sequence name, inner key = ligand name.
This lets you specify different constraints for different conditions
(e.g. "free" vs. a specific bound ligand).

Use `"*"` as the ligand key to provide a fallback for any ligands not
explicitly listed.  This avoids repeating the same set of pairs for
multiple ligands.

```python
additional_pairs = {
    "Seq1": {
        "free": [("A", 1, "A", 40)],       # only when no ligand
        "L01":  [("A", 5, "B", 10)],       # only with L01 bound
        "*":    [("A", 2, "A", 39)],        # fallback for all other ligands
    },
}
```

### What happens when constraints are defined

When `additional_pairs` resolves to a **non-empty** list for a given
(sequence × ligand) entry:

- The YAML and SLURM filenames include a `_constrained` suffix
  (e.g. `Seq1_L01_constrained.yaml`).
- Both an **unconstrained** and a **constrained** YAML file are
  generated for that entry, allowing direct comparison of the two
  conditions.