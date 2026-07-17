# aptaBoltz — Nucleic-Acid / Aptamer Analysis Pipeline for Boltz-2

A modular pipeline for:

1. **Boltz-2 structure prediction** — generate YAML inputs, submit SLURM
   jobs, process prediction outputs (PyMOL + RDKit), and visualise results
   (R).
2. **MD force-field parameterisation** — RESP2 charges, GAFF2 atom types,
   solvated Amber systems ready for molecular-dynamics production runs
   ([§8](#8-md-force-field-parameterisation-pipeline)).
3. **MD preparation and minimisation** — 11-step Amber production protocol
   (minimisation → equilibration → 100 ns NPT → final minimisation) as a
   SLURM array job ([§9](#9-amber-md-preparation-and-minimisation-pipeline)).

Designed primarily for **DNA and RNA aptamers** with small-molecule ligands.
Basic protein support is available — see [§7](#7-molecule-type-support).

### Scope

| This pipeline DOES: | This pipeline does NOT: |
|---------------------|-------------------------|
| Generate Boltz-2 YAML input files | Run Boltz-2 predictions themselves (submit to a cluster) |
| Generate SLURM submission scripts | Run MD production simulations (you provide the input) |
| Process prediction CIF outputs (alignment, colouring, ligand extraction, chirality) | Perform QM geometry optimisation (provide the optimised XYZ) |
| Aggregate confidence, PAE, PDE, pLDDT metrics | Host the Boltz-2 MSA server |
| Visualise PAE/PDE heatmaps and pLDDT per-residue plots | Install ORCA, AmberTools, or Multiwfn for you |
| Prepare AMBER-compatible MD systems from Boltz-2 predictions | |

### Repository Structure

```
boltz2_utils/               Python package (imported by all projects)
  generate_boltz2_yaml.py     YAML input builder
  generate_boltz2_slurm.py    SLURM script generator
  process_boltz_results.py    Output processing (PyMOL + RDKit)
  base.py                     PyMOL base styling

boltz_R_utils/              R scripts (sourced by all projects)
  processor.R                 Confidence, PAE, PDE, pLDDT reading + plotting
  install_packages.R          One-time R + conda environment setup

templates/                  Project-specific files with PLACEHOLDER_ values
  (root)                      Boltz-2 prediction scripts (§3–§6)
  MD/                         MD parameterisation + production scripts (§8–§9)
  input_file_generator.py     Define sequences, ligands, constraints (§3)
  output_file_processing.py   Python post-processing (§5)
  process.R                   R visualisation (§6)
  yaml/                       (placeholder — generated YAML files)
  msa/                        (placeholder — custom MSAs)

  MD/
    orca_steps_wsl.sh           ORCA single-point calculations (§8.4)
    multiwfn_steps_windows.ps1  RESP2 charge fitting, Windows (§8.5.1)
    multiwfn_steps_linux.sh     RESP2 charge fitting, Linux (§8.5.2)
    ligand_prep.sh              AmberTools GAFF2 parameterisation (§8.6)
    leap.sh                     Amber tLEAP solvation + ionisation (§8.8)
    python/
      generate_md_slurm.py        SLURM generator for PrepAndMin.slurm (§9.2)
      generate_pmemd_inputs.py    Amber pmemd input file generator (§9.5)
      model_prep.py               Boltz-2 → MD PDB preparation (§8.7)
      pymol_utils.py              PyMOL visualisation utilities (§9.11)
    R/
      ions.R                      Ion concentration calculator (§8.9)
      md_utils.R                  MD trajectory analysis: RMSD, atomic fluctuation (§9.10)
    pmemd/
      in/
        archive/                    Previous versions of generated .in files
        step1.in — step10.in       Amber MD input control files (§9)
        final_min.in                Post-production minimisation (§9)
    slurm/
      PrepAndMin.slurm            Amber MD prep + minimisation SLURM array template (§9.1)

CSS/                        Real usage example (corticosteron-specific aptamers)
  input_file_generator.py     CSS-specific sequences: CSS1, CSS2, CSS3
  output_file_processing.py   CSS-specific ligands: C0R, HCY, 11DC
  process.R
  MD/                         Full MD pipeline reference
    python/
      pymol_utils.py            PyMOL utilities (copy of template, §9.11)
      final_min_pymol.py        Final-min session generator (§9.11)
```

> **Note on molecule types:** This pipeline works natively with **DNA** and
> **RNA** (nucleic-acid aptamers).  Basic **protein** support is available:
> the YAML generator accepts `molecule_type="protein"`, and PyMOL styling
> colours protein chains in a uniform muted colour.  However, the stem-contact
> logic (Watson-Crick base pairing) is skipped for proteins, and the R pLDDT
> visualiser labels residues using single-letter codes read directly from the
> YAML sequence.

---

## Boltz-2 Prediction Pipeline Overview

```
┌──────────────────────────────────────────────────────────────────┐
│  input_file_generator.py                                         │
│  sequences + ligands + constraints → YAML + SLURM                │
│  Outputs:  project/yaml/*.yaml, project/*.slurm                  │
└──────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼  (transfer to cluster)
┌──────────────────────────────────────────────────────────────────┐
│  Cluster: sbatch → Boltz-2                                       │
│  Outputs:  {job}/boltz_results/input/predictions/*.cif           │
│            *_confidence.json, *_pae.npz, *_pde.npz, *_plddt.npz  │
└──────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼  (copy results back)
┌──────────────────────────────────────────────────────────────────┐
│  output_file_processing.py  (PyMOL + RDKit)                      │
│  Alignment, colouring, ligand extraction, chirality, CSV export  │
│  Outputs:  *.pse, *.sdf, *_chirality.png                         │
│            *_confidence.csv, *_pae.csv, *_pde.csv, *_plddt.csv   │
└──────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼
┌────────────────────────────────────────────────────────────────────┐
│  process.R  (R visualisation)                                      │
│  Confidence table, PAE/PDE heatmaps, pLDDT plots                   │
│  Outputs:  confidence_table.html, *_PAE.png, *_PDE.png, *_plddt.png│
└────────────────────────────────────────────────────────────────────┘
```

For the MD parameterisation pipeline, see [§8](#8-md-force-field-parameterisation-pipeline).

---

## Table of Contents

- [Scope](#scope)
- [Repository Structure](#repository-structure)
- [Boltz-2 Prediction Pipeline Overview](#boltz-2-prediction-pipeline-overview)
- [1. Prerequisites and Installation](#1-prerequisites-and-installation)
  - [1.1 Platform Requirements](#11-platform-requirements)
  - [1.2 Installation Options](#12-installation-options)
  - [1.3 R Environment Setup](#13-r-environment-setup)
  - [1.4 Verify Installation](#14-verify-installation)
- [2. Setting Up a New Project](#2-setting-up-a-new-project)
  - [2.1 Boltz-2 Prediction Project](#21-boltz-2-prediction-project)
  - [2.2 MD Parameterisation Project](#22-md-parameterisation-project)
- [3. Generating Input Files (YAML + SLURM)](#3-generating-input-files-yaml--slurm)
  - [3.1 Step-by-Step Editing Guide](#31-step-by-step-editing-guide)
  - [3.2 Generated Files](#32-generated-files)
- [4. Running Boltz-2 on a Cluster](#4-running-boltz-2-on-a-cluster)
  - [4.1 Files to Transfer](#41-files-to-transfer)
  - [4.2 Expected Outputs (per job)](#42-expected-outputs-per-job)
  - [4.3 Common Submission Patterns](#43-common-submission-patterns)
- [5. Post-Processing: Python (PyMOL + RDKit)](#5-post-processing-python-pymol--rdkit)
  - [5.1 Configuration](#51-configuration)
  - [5.2 Usage](#52-usage)
  - [5.3 What the Script Does](#53-what-the-script-does)
  - [5.4 Generated Files](#54-generated-files)
  - [5.5 Controlling Which Job Directories Are Processed](#55-controlling-which-job-directories-are-processed)
- [6. Post-Processing: R Visualisation](#6-post-processing-r-visualisation)
  - [6.1 Minimal Setup](#61-minimal-setup)
  - [6.2 Generated Files](#62-generated-files)
  - [6.3 Function Reference](#63-function-reference)
- [7. Molecule Type Support](#7-molecule-type-support)
  - [7.1 Troubleshooting — Boltz-2 Prediction Pipeline](#71-troubleshooting--boltz-2-prediction-pipeline)
- [8. MD Force-Field Parameterisation Pipeline](#8-md-force-field-parameterisation-pipeline)
  - [8.1 Pipeline Overview](#81-pipeline-overview)
  - [8.2 Platform Notes](#82-platform-notes)
  - [8.3 Step 0 — QM Geometry Optimisation (external)](#83-step-0--qm-geometry-optimisation-external)
  - [8.4 Step 1 — ORCA Single-Point Calculations](#84-step-1--orca-single-point-calculations)
  - [8.5 Step 2 — RESP2 Charge Fitting](#85-step-2--resp2-charge-fitting)
    - [8.5.1 Windows (recommended): `multiwfn_steps_windows.ps1`](#851-windows-recommended-multiwfn_steps_windowsps1)
    - [8.5.2 Linux alternative: `multiwfn_steps_linux.sh`](#852-linux-alternative-multiwfn_steps_linuxsh)
  - [8.6 Step 3 — Ligand Parameterisation](#86-step-3--ligand-parameterisation)
  - [8.7 Step 4 — Model Preparation for MD](#87-step-4--model-preparation-for-md)
  - [8.8 Step 5 — System Solvation \& Ionisation](#88-step-5--system-solvation--ionisation)
  - [8.9 Ion Concentration Calculator](#89-ion-concentration-calculator)
  - [8.10 Required Files Summary](#810-required-files-summary)
  - [8.11 Dependencies Table](#811-dependencies-table)
  - [8.12 Troubleshooting](#812-troubleshooting)
- [9. Amber MD Preparation and Minimisation Pipeline](#9-amber-md-preparation-and-minimisation-pipeline)
  - [9.1 Pipeline Overview](#91-pipeline-overview)
  - [9.2 SLURM Script Generation](#92-slurm-script-generation)
  - [9.3 The 11-Step Protocol](#93-the-11-step-protocol)
  - [9.4 Input File Descriptions](#94-input-file-descriptions)
  - [9.5 Input File Generation](#95-input-file-generation)
  - [9.6 Output Files (per task)](#96-output-files-per-task)
  - [9.7 SLURM Configuration](#97-slurm-configuration)
  - [9.8 Customisation Guide](#98-customisation-guide)
  - [9.9 Troubleshooting](#99-troubleshooting)
  - [9.10 MD Trajectory Analysis](#910-md-trajectory-analysis)
  - [9.11 PyMOL Visualisation Utilities](#911-pymol-visualisation-utilities)
- [10. File Reference](#10-file-reference)
  - [Shared utilities (imported, do not edit for each project)](#shared-utilities-imported-do-not-edit-for-each-project)
  - [Project-specific templates (copy and edit per project)](#project-specific-templates-copy-and-edit-per-project)
- [11. Appendix: Contact Constraints Format](#11-appendix-contact-constraints-format)
  - [Format A: Flat list (backward-compatible)](#format-a-flat-list-backward-compatible)
  - [Format B: Per-sequence dictionary](#format-b-per-sequence-dictionary)
  - [Format C: Per-sequence + per-ligand dictionary (with wildcard)](#format-c-per-sequence--per-ligand-dictionary-with-wildcard)
  - [What happens when constraints are defined](#what-happens-when-constraints-are-defined)

---

## 1. Prerequisites and Installation

### 1.1 Platform Requirements

| Component | Platform | Installation method |
|-----------|----------|-------------------|
| Python (≥ 3.10) | Any | conda (via `environment.yml`) or system Python |
| boltz2_utils | Any (Python) | `pip install -e ./boltz2_utils` or included in conda env |
| PyMOL | Any | `conda install -c conda-forge pymol-open-source` |
| RDKit | Any | Included in `environment.yml` via conda-forge |
| R (≥ 4.x) | Any | System R + `install_packages.R` |
| R packages | Any | `data.table`, `ggplot2`, `DT`, `khroma`, `jsonlite`, `yaml`, `bslib`, `reticulate` |
| ORCA | Linux/WSL | Download from orcaforum.kofo.mpg.de (for §8) |
| Multiwfn | Windows/Linux | Download from sobereva.com/multiwfn (for §8) |
| AmberTools | Linux/WSL | conda install or Amber website (for §8) |

### 1.2 Installation Options

**Option A — conda environment (recommended):**

```bash
# Create the 'boltz' environment with all Python dependencies
conda env create -f environment.yml
conda activate boltz
```

This installs Python 3.10+, NumPy, pandas, natsort, RDKit, PyMOL
(`pymol-open-source`), and the `boltz2_utils` package (via `pip -e`).

**Option B — standalone pip install:**

If you already have a Python environment with PyMOL and RDKit:

```bash
pip install -e ./boltz2_utils
```

This makes the following modules importable:
- `boltz2_utils.process_boltz_results` — prediction output processing
- `boltz2_utils.generate_boltz2_yaml` — YAML input file generation
- `boltz2_utils.generate_boltz2_slurm` — SLURM job script generation

### 1.3 R Environment Setup

Run the one-time setup script to install required R packages and create
the `boltz` conda environment (idempotent — safe to re-run):

```r
source("boltz_R_utils/install_packages.R")
```

This script:
1. Installs any missing R packages from CRAN (`data.table`, `ggplot2`,
   `DT`, `khroma`, `jsonlite`, `yaml`, `bslib`, `reticulate`).
2. Creates the `boltz` conda environment from `environment.yml` (if it
   does not yet exist).
3. Configures `reticulate` to use the `boltz` conda environment.

### 1.4 Verify Installation

```bash
# Python / boltz2_utils
python -c "from boltz2_utils.generate_boltz2_yaml import generate_boltz2_yamls; print('boltz2_utils OK')"

# R + reticulate + conda env
Rscript -e "source('boltz_R_utils/install_packages.R')"
```

---

## 2. Setting Up a New Project

### 2.1 Boltz-2 Prediction Project

For each new structure-prediction project (e.g. `my_project`):

1. Create the project directory structure:
   ```bash
   mkdir my_project
   mkdir my_project/yaml
   mkdir my_project/msa        # optional — only if using custom MSA
   ```

2. Copy the three Boltz-2 template files:
   ```bash
   cp templates/input_file_generator.py     my_project/
   cp templates/output_file_processing.py   my_project/
   cp templates/process.R                   my_project/
   ```

   Each file contains `PLACEHOLDER_` values that you replace with your
   project's names and sequences (see [§3](#3-generating-input-files-yaml--slurm)
   through [§6](#6-post-processing-r-visualisation)).

3. The shared utilities (`boltz2_utils/`, `boltz_R_utils/`) remain in
   the repository root — do **not** copy them.  They are imported /
   sourced automatically.

### 2.2 MD Parameterisation Project

For the MD force-field parameterisation workflow (see
[§8](#8-md-force-field-parameterisation-pipeline) and
[§9](#9-amber-md-preparation-and-minimisation-pipeline) for full details):

```bash
mkdir -p my_project/MD/QM
mkdir -p my_project/MD/ff
mkdir -p my_project/MD/pmemd/in
mkdir -p my_project/MD/slurm
# Copy ALL MD template files (parameterisation + production)
cp -r templates/MD/*       my_project/MD/
cp templates/MD/R/ions.R    my_project/MD/R/
```

This copies all MD pipeline scripts into the correct subdirectories:

| Copy command | Destination | Files |
|-------------|------------|-------|
| `cp -r templates/MD/*` | `my_project/MD/` | `leap.sh`, `ligand_prep.sh`, `python/generate_pmemd_inputs.py`, `python/generate_md_slurm.py`, `python/model_prep.py`, `orca_steps_wsl.sh`, `multiwfn_steps_*.ps1/.sh`, `pmemd/in/archive/`, `slurm/PrepAndMin.slurm` |
| `cp templates/MD/R/ions.R` | `my_project/MD/R/` | `ions.R` |

---

## 3. Generating Input Files (YAML + SLURM)

**Script:** `templates/input_file_generator.py` → copy to `{project}/` and edit.

**Edits required:** `project`, `molecule_type`, `SEQUENCES`, `LIGANDS`,
`additional_pairs`, `custom_msa`.

### 3.1 Step-by-Step Editing Guide

**Step 1 — Set project name, molecule type, and sequences:**

```python
project = "my_project"
molecule_type = "dna"       # "dna", "rna", or "protein"
SEQUENCES = {
    "Seq1": "GGGACGACGCCCGCATGTTCCATGGATAGTCTTGACTAGTCGTCCC",
    "Seq2": "GGGACGACTAGCGTATGCGCCAGAAGTATACGAGGATAGTCGTCCC",
}
stem_length = 8
max_distance = 4.0           # stem contact distance (Å), DNA/RNA only
```

**Step 2 — (Optional) Define contact constraints:**

Set `additional_pairs` to add extra residue-pair contacts.  Three formats
are accepted (flat list, per-sequence dict, per-sequence+per-ligand dict
with `"*"` wildcard).  See [Appendix 11](#11-appendix-contact-constraints-format).
Set to `None` to omit.

**Step 3 — Define ligands:**

Each ligand is a dict with `"ccd"` (PDB CCD code) and/or `"smiles"` key:

```python
LIGANDS = {
    "free": {"ccd": None, "smiles": None},               # no ligand
    "L01":  {"ccd": "C0R", "smiles": None},              # CCD-based
    "L02":  {"ccd": None, "smiles": "C1=CC=C2C(=C1)C(=O)C3=C(C2=O)C4=CC=CC=C4C3=O"},
}
```

**Step 4 — Build entry list:** (automatic — no action required)

**Step 5 — Toggle custom MSA:**

```python
custom_msa = False     # True → place {seq_name}.a3m in {project}/msa/
```

**Steps 6–7 — Generate files:**

```bash
python my_project/input_file_generator.py
```

### 3.2 Generated Files

| Path | Description |
|------|-------------|
| `{project}/yaml/*.yaml` | Boltz-2 input files (one per seq × ligand × variant) |
| `{project}/*.slurm` | SLURM submission scripts (one per YAML file) |

**Naming conventions:**

| Variant | Example filename |
|---------|-----------------|
| Standard | `Seq1_L01.yaml` / `.slurm` |
| With custom MSA | `Seq1_L01_custom_msa.yaml` / `.slurm` |
| With contact constraints | `Seq1_L01_constrained.yaml` / `.slurm` |
| MSA + constraints | `Seq1_L01_custom_msa_constrained.yaml` / `.slurm` |

---

## 4. Running Boltz-2 on a Cluster

The pipeline generates SLURM scripts; you transfer them to the cluster and
submit via `sbatch`.  Adjust paths below to match your cluster setup.

### 4.1 Files to Transfer

| Local path | Cluster path (example) |
|------------|----------------------|
| `{project}/` | `scratch/boltz/projects/` |
| `{project}/yaml/` | `scratch/boltz/projects/{project}/yaml/` |
| `{project}/msa/` | `scratch/boltz/projects/{project}/msa/` (if used) |
| `{project}/*.slurm` | `scratch/boltz/projects/{project}/` |

### 4.2 Expected Outputs (per job)

On successful completion, each job directory contains:

| File(s) | Contents |
|---------|----------|
| `boltz_results_input/predictions/input/*.cif` | Predicted structures (all models) |
| `*_confidence.json` | per-model confidence metrics |
| `*_pae.npz` | Predicted aligned error matrix |
| `*_pde.npz` | Predicted distance error matrix |
| `*_plddt.npz` | Per-residue pLDDT scores |

### 4.3 Common Submission Patterns

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

**Script:** `templates/output_file_processing.py` → copy to `{project}/` and edit.

**Edits required:** `project`, `models`, `LGD_DICT`, `NAME_MAP`, `suffixes`.

**Inputs required:** Boltz-2 prediction output directories (see §4.2).

### 5.1 Configuration

```python
# Display name → CCD code (used for reference fetching + chirality labels)
LGD_DICT = {"Aptamer1": "L01", "Aptamer2": "L02", "free": "free"}

# CCD code → display name (reverse mapping for file naming)
NAME_MAP = {"L01": "Aptamer1", "L02": "Aptamer2"}

project = "my_project"
models = ["Seq1", "Seq2"]       # sequence names matching SEQUENCES keys
suffixes = ["_constrained", ""] # typically both constrained and unconstrained
```

### 5.2 Usage

```bash
python my_project/output_file_processing.py
```

### 5.3 What the Script Does

For each (model × ligand × suffix × job-directory) combination:

1. **Loads CIF predictions** into PyMOL (all models from the job).
2. **Aligns and saves** a multi-object CIF and a `.pse` session file.
3. **Colours** nucleic-acid chains by base (DNA/RNA), ligand by element,
   protein chains in a muted uniform colour, constrained residues by
   component (pastel palette), and all residues by pLDDT (red→blue).
4. **Extracts ligands** as SDF files from predicted structures and fetches
   the corresponding reference ligand from the PDB.
5. **Analyses chirality** — loads SDFs into RDKit, detects chiral centres,
   and draws an annotated grid PNG.
6. **Aggregates metrics** — reads confidence JSON, PAE/PDE NPZ, and pLDDT
   NPZ files; saves as CSV files.

### 5.4 Generated Files

| File(s) | Description |
|---------|-------------|
| `{experiment}_{job}.cif` | Aligned multi-model CIF |
| `{experiment}_{job}.pse` | PyMOL session (base-coloured, aligned) |
| `{experiment}_{job}_plddt.pse` | PyMOL session (pLDDT-coloured) |
| `{ligand_name}.sdf` | Extracted ligand SDF |
| `{ligand_name}_ref.sdf` | Reference ligand from PDB |
| `{experiment}_chirality.png` | Chirality annotation grid |
| `{experiment}_{job}_confidence.csv` | Confidence metrics (one row per model) |
| `{experiment}_{job}_pae.csv` | PAE matrix (one row per residue) |
| `{experiment}_{job}_pde.csv` | PDE matrix (one row per residue) |
| `{experiment}_{job}_plddt.csv` | pLDDT per residue |

### 5.5 Controlling Which Job Directories Are Processed

The `job_dirs` parameter of `process_single_experiment()` accepts:

| Value | Behaviour |
|-------|-----------|
| `None` | Process **all** J-directories |
| `"last"` | Process only the highest-numbered directory |
| `"J1121077"` | Process only that specific directory |
| `["J1121077", "J1123593"]` | Process only the listed directories |

---

## 6. Post-Processing: R Visualisation

**Script:** `templates/process.R` → copy to `{project}/` and edit.

**Edits required:** `project`, `ligand_number`.

**Inputs required:** CSV files generated by §5 (scanned automatically).

### 6.1 Minimal Setup

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

Set `ligand_number` to the residue index at which the ligand region begins
(i.e. the number of nucleic-acid / protein residues).  Example: 46-mer DNA
aptamer → `ligand_number = 47`.

### 6.2 Generated Files

| File | Description |
|------|-------------|
| `{project}/confidence_table.html` | Interactive confidence summary (DT table) |
| `{project}/*_{region}_pae.png` | PAE heatmap, one per region (DNA, Ligand) |
| `{project}/*_{region}_pde.png` | PDE heatmap, one per region (DNA, Ligand) |
| `{project}/*_plddt.png` | pLDDT per-residue faceted line plot |

### 6.3 Function Reference

| Function | Input | Output |
|----------|-------|--------|
| `process_confidence()` | All `*_confidence.csv` files | `data.table`, one row per model |
| `table_confidence()` | The returned `data.table` | HTML file `confidence_table.html` |
| `process_pxe(type=...)` | All `*_pae.csv` or `*_pde.csv` files | Named list of `data.table`s (one per experiment) |
| `plot_pxe()` | One element from `process_pxe()` | Faceted PNG per region (DNA / Ligand) |
| `plot_plddt()` | All `*_plddt.csv` files in project | Faceted line-plot PNG per experiment |

---

## 7. Molecule Type Support

The `molecule_type` parameter in `input_file_generator.py` controls how
the pipeline handles different biopolymers:

| Type | YAML key | Stem contacts | PyMOL styling | pLDDT labels |
|------|----------|---------------|---------------|--------------|
| `"dna"` | `dna` | ✓ Watson-Crick base-pair contacts | Base-coloured (A, G, T, C) | Single-letter DNA |
| `"rna"` | `rna` | ✓ Watson-Crick base-pair contacts | Base-coloured (A, G, U, C) | Single-letter RNA |
| `"protein"` | `protein` | ✗ (no base-pair logic) | Uniform muted (`warmpink`) | 1-letter amino-acid codes |

- **DNA and RNA** are fully supported: stem contacts, base-identity
  colouring, and per-residue labelling all work natively.
- **Protein** support is basic: YAML generation and PyMOL styling work,
  but stem-contact constraints are skipped.  Protein-specific features
  (secondary-structure colouring, domain parsing) are not implemented.

The ligand-extraction and chirality-analysis components are
molecule-type agnostic and work for any system.

### 7.1 Troubleshooting — Boltz-2 Prediction Pipeline

| Issue | Likely Cause | Solution |
|-------|-------------|----------|
| No YAML files generated | Python error in the script | Run `python {project}/input_file_generator.py` and read the traceback |
| No SLURM files generated | YAML step did not complete | Fix YAML generation first; SLURM scripts depend on YAML output |
| PyMOL: `command not found` | `pymol-open-source` not installed | `conda install -c conda-forge pymol-open-source` |
| RDKit: chirality grid empty | Ligand SDF extraction failed | Verify `LGD_DICT` / `NAME_MAP` mapping matches the YAML ligand names |
| R: `object 'process_confidence' not found` | `processor.R` not sourced | Check path in `source('boltz_R_utils/processor.R')` |
| R: conda environment not active | `install_packages.R` not run first | Run `source('boltz_R_utils/install_packages.R')` once |
| R: CSV files not found | Python post-processing not run | Run `output_file_processing.py` first |
| `_constrained` jobs missing | `additional_pairs` resolves to empty | Check contact-constraint format in Appendix 11 |
| Custom MSA not used | File not found or named wrong | Place `{seq_name}.a3m` in `{project}/msa/` and set `custom_msa = True` |
| PyMOL session colours wrong | Unexpected chain/residue naming in CIF | Check the CIF output; may need adjustments in `base.py` |

---


## 8. MD Force-Field Parameterisation Pipeline

In addition to the Boltz-2 prediction pipeline, the `templates/` directory
contains scripts for generating AMBER-compatible force-field parameters
from quantum-mechanical (QM) calculations.  This is a separate workflow used
to prepare **ligand parameters** for classical molecular-dynamics (MD)
simulations that complement the Boltz-2 structure predictions.

### 8.1 Pipeline Overview

```
                            ┌──────────────────────────────────────────┐
                            │  QM Geometry Optimisation (external*)    │
                            │  →  QM/<MOL>.xyz  (or .mol2)             │
                            └──────────────────────────────────────────┘
                                            │
                                            ▼
                            ┌──────────────────────────────────────────┐
                            │  orca_steps_wsl.sh           (WSL/Linux) │
                            │  2 × SP (gas + solvent CPCM-SMD)         │
                            │  →  SP_gas.molden                        │
                            │  →  SP_solv.molden                       │
                            └──────────────────────────────────────────┘
                                            │
                          ┌─────────────────┴─────────────────┐
                          ▼                                   ▼
         ┌──────────────────────────────┐    ┌──────────────────────────────┐
         │  multiwfn_steps_windows.ps1  │    │  multiwfn_steps_linux.sh     │
         │  (Windows, recommended)      │    │  (Linux/Unix, alternative)   │
         │  →  QM/<mol>.chg (RESP2)     │    │  →  QM/<mol>.chg (RESP2)     │
         └──────────────────────────────┘    └──────────────────────────────┘
                          │                                   │
                          └─────────────────┬─────────────────┘
                                            │
                                            ▼
                            ┌──────────────────────────────────────────┐
                            │  ligand_prep.sh             (WSL/Linux)  │
                            │  AmberTools GAFF2 + RESP2 charges        │
                            │  →  ff/<lig>_resp.prepin                 │
                            │  →  ff/<lig>_resp.frcmod                 │
                            │  →  ff/<lig>_resp.pdb                    │
                            └──────────────────────────────────────────┘
                                            │
                                            ▼
                            ┌──────────────────────────────────────────┐
                            │  model_prep.py              (PyMOL)      │
                            │  Boltz-2 prediction CIFs → MD-ready PDB  │
                            │  →  pdb_for_md/<PREFIX>_<LIGAND>_constrained│
                            │      /input_model_{i}.pdb                │
                            └──────────────────────────────────────────┘
                                            │
                                            ▼
                            ┌──────────────────────────────────────────┐
                            │  leap.sh                   (WSL/Linux)   │
                            │  tLEAP solvation + ionisation            │
                            │  →  leap/<PREFIX>_<LIGAND>_constrained/   │
                            │      cplx_{i}.prmtop / .rst7 / .pdb      │
                            └──────────────────────────────────────────┘
                                            │
                                            ▼
                                  MD Production (user-defined)

*The initial optimised geometry can come from any QM package (ORCA,
Gaussian, etc.) or from public databases.  Non-optimised structures
may be used but RESP2 charges will be less reliable.
```

### 8.2 Platform Notes

The pipeline involves tools that span multiple operating systems.  Three
workflows are possible:

| Workflow | ORCA | Multiwfn | AmberTools | Key consideration |
|----------|------|----------|------------|-------------------|
| **WSL + Windows** (recommended) | WSL/Linux | Windows (.ps1) | WSL/Linux | Use `multiwfn_steps_windows.ps1` for the latest Multiwfn developments |
| **Linux-only** | Linux | Linux (.sh) | Linux | Single OS; use `multiwfn_steps_linux.sh` |
| **Windows-only** | Windows | Windows | WSL | ORCA on Windows via cmd; AmberTools via WSL |

The repository defaults to the **WSL + Windows** split (most common for
users running ORCA and AmberTools on a Linux/WSL GPU cluster and Multiwfn
on a local Windows machine).  If you work on a single Linux system, use
`multiwfn_steps_linux.sh` instead of the PowerShell script.

### 8.3 Step 0 — QM Geometry Optimisation (external)

Before running the RESP2 charge pipeline, you need an optimised 3D geometry
of your ligand.  This step is **not** covered by a template because the
choice of QM package and workflow is project-specific.

**Expected file:** `QM/<mol>_opt.xyz` — optimised geometry in XYZ format.
(Both `.xyz` and `.mol2` are accepted — see Step 3 below for the fallback.)

If you already have a CCD ligand from the PDB, you may use its geometry
directly.  For novel ligands, run a geometry optimisation in your preferred
QM package (ORCA, Gaussian, psi4, etc.) and place the output XYZ in `QM/`.

### 8.4 Step 1 — ORCA Single-Point Calculations

**Script:** `templates/orca_steps_wsl.sh` → copy to project and edit.

Runs two consecutive single-point energy calculations with ORCA on an
optimised ligand geometry (XYZ format):

- **Gas phase** (B3LYP-D3/def2-TZVP with RIJCOSX)
- **Solvent phase** (same level, CPCM-SMD continuum solvation)

Both calculations produce Molden-format files with an Nval valence-electron
header (required for Multiwfn RESP fitting with def2 ECP basis sets).

**Usage:**
```bash
./orca_steps_wsl.sh QM/<mol>_opt.xyz [charge] [multiplicity] [solvent]
```

**Inputs required:** `QM/<mol>_opt.xyz`

**Outputs:**
- `SP_gas.molden` — Gas-phase Molden file
- `SP_solv.molden` — Solvent-phase Molden file

**Edits required in template:** ORCA paths (`$ORCA`, `$orca_2mkl`), number
of processors (`$nprocs`), memory (`$maxcore`), and the `$keyword_SP` if a
different level of theory is desired.

### 8.5 Step 2 — RESP2 Charge Fitting

Two equivalent scripts are provided:

#### 8.5.1 Windows (recommended): `multiwfn_steps_windows.ps1`

Reads the gas- and solvent-phase Molden files (copied from WSL/Linux
to the Windows `<MoldenDir>`) and invokes Multiwfn via a batch wrapper
to compute RESP charges for each phase.  It then combines the two charge
sets into RESP2 charges using the mixing parameter delta (default: 0.5,
i.e. RESP2-0.5).

**Usage:**
```powershell
.\multiwfn_steps_windows.ps1 <input_name> [-delta <value>]
```

**Inputs required:**
- `SP_gas.molden` — copied from WSL to Windows
- `SP_solv.molden` — copied from WSL to Windows

**Deliverable:** `QM/<input_name>.chg` — RESP2 charges for the ligand
(also produces intermediate `gas.chg` and `solv.chg`)

**Edits required in template:** Multiwfn executable path, Molden directory,
output directory.

#### 8.5.2 Linux alternative: `multiwfn_steps_linux.sh`

A Linux/Bash equivalent of the PowerShell script.  Uses stdin piping to
interact with Multiwfn (no batch file wrapper needed).

**Usage:**
```bash
./multiwfn_steps_linux.sh <input_name> [delta]
```

**Inputs / Outputs:** Same as the Windows version.

**Note:** Using Multiwfn on Windows (`multiwfn_steps_windows.ps1`) is
recommended for the latest developments.  This Linux script is provided
as a convenience for single-OS setups.

### 8.6 Step 3 — Ligand Parameterisation

**Script:** `templates/ligand_prep.sh` → copy to project and edit.

Takes the optimised `.mol2` (or `.xyz` as fallback) and RESP2 `.chg` files
from the `QM/` directory and runs the AmberTools pipeline:

1. **Antechamber** — assigns GAFF2 atom types and BCC charges.
2. **Charge substitution** — replaces BCC charges with RESP2 charges.
3. **Prepgen** — generates an Amber prepin library file.
4. **Parmchk2** — checks for missing force-field parameters (frcmod).

If `QM/<MOL>.mol2` is not found, the script automatically looks for
`QM/<MOL>.xyz` and attempts conversion via antechamber or OpenBabel.

**Usage:**
```bash
./ligand_prep.sh
```
(Edit `MOL`, `NC`, and `S` at the top of the script for your ligand.)

**Inputs required:**
- `QM/<MOL>.mol2` (preferred) or `QM/<MOL>.xyz` (fallback) — optimised
  ligand structure
- `QM/<MOL>_opt.chg` — RESP2 charges from Step 2

**Outputs:**
- `ff/<MOL>_resp.prepin` — Amber prepin library
- `ff/<MOL>_resp.frcmod` — Force-field modifications
- `ff/<MOL>_resp.pdb` — PDB-like file for tLeap complex preparation

### 8.7 Step 4 — Model Preparation for MD

**Script:** `templates/MD/python/model_prep.py` → copy to project and edit.

Converts Boltz-2 prediction CIF files into PDB models suitable for
Amber tLEAP.  This involves loading each prediction into PyMOL,
performing system-specific modifications, optionally aligning the
reference ligand, and saving ready-for-MD PDB files.

**System-specific modifications (DNA aptamer example):**
- Removal of the 5' terminal phosphate
- Capping the terminal hydroxyl
- Aligning the reference ligand PDB (from `ligand_prep.sh`) onto the
  predicted ligand
- Reassigning chains (DNA → chain A, ligand → chain B)

When `--lgd free` is used, all ligand-related operations (ligand
alignment, reference PDB loading, chain B reassignment) are skipped.
The script processes the DNA-only structure directly.

**CAUTION:** The PyMOL operations in `prep_model()` are specific to
**DNA aptamers** (5' phosphate removal, terminal hydroxyl capping).
If your system is RNA or protein, edit the marked block in the script
before running.

**Usage (command line):**
```bash
# With ligand
cd my_project/MD
python python/model_prep.py --seq <PREFIX> --lgd <LIGAND> --lgd-file <LGD_FILE> [options]

# Without ligand (apo predictions)
cd my_project/MD
python python/model_prep.py --seq <PREFIX> --lgd free [options]
```

> **Important:** Run this script from the `MD/` directory (e.g. `CSS/MD/`), not
> from `MD/python/`.  The script uses relative paths (`../` for Boltz-2
> prediction output and `ff/` for ligand parameter files) that resolve
> correctly only from `MD/`.

Auto-detects the latest Boltz-2 job directory if `--job` is omitted.

**Inputs required:**
- Boltz-2 prediction output in `../<PREFIX>_<LIGAND><suffix>/<JOB>/`
- `ff/<lgd_file>.pdb` — reference ligand PDB from `ligand_prep.sh`
  (not required when `--lgd free`)

**Outputs:**
- `pdb_for_md/<PREFIX>_<LIGAND>_constrained/input_model_{i}.pdb`
- `pdb_for_md/<PREFIX>_<LIGAND>_constrained/README.txt`

**Edits required in template:** Sequence name (`--seq`), ligand name
(`--lgd`), reference PDB name (`--lgd-file`, omit when using `--lgd free`),
model range, and the PyMOL modification block (marked with `EDIT THIS BLOCK`).

### 8.8 Step 5 — System Solvation & Ionisation

**Script:** `templates/leap.sh` → copy to project and edit.

Creates a single tLEAP input script that:
- Loads force fields (GAFF2, DNA.OL21, OPC water)
- Loads the ligand prepin and frcmod (skipped with `-l free`)
- Loads each PDB model from `pdb_for_md/<PREFIX>_<LIGAND>_constrained/`
- Solvates in an OPC water box (14.0 Å buffer)
- Adds Na⁺/Mg²⁺/Cl⁻ ions with automatic Cl⁻ distribution (MgCl₂ + NaCl)

Runs tLEAP on the combined input, then converts each rst7 to PDB with ambpdb.

**Usage:**
```bash
./leap.sh -p <prefix> -l <ligand> -s <start> -e <end> -n <na> -m <mg> -c <cl>
```

**Examples:**
```bash
# CSS1, ligand hcy, models 0-24, 65 Na⁺, 1 Mg²⁺, 22 Cl⁻
./leap.sh -p CSS1 -l hcy -s 0 -e 24 -n 65 -m 1 -c 22

# No-ligand mode (no prepin/frcmod needed)
./leap.sh -p CSS1 -l free -s 0 -e 24 -n 65 -m 1 -c 22
```

**Inputs required:**
- `ff/<ligand>_resp.prepin` — from Step 3 (not required with `-l free`)
- `ff/<ligand>_resp.frcmod` — from Step 3 (not required with `-l free`)
- `pdb_for_md/<PREFIX>_<LIGAND>_constrained/input_model_{i}.pdb` — from Step 4

**Deliverables** (in `leap/<PREFIX>_<LIGAND>_constrained/`):
- `tleap.in` — Combined tLEAP input script
- `cplx_{i}.prmtop` — Amber topology (one per model)
- `cplx_{i}.rst7` — Amber restart/coords (one per model)
- `cplx_{i}.pdb` — PDB conversion (one per model)

**Edits required in template:** Default `-p` prefix, default `-l` ligand
name, and model/ion values can be changed via flags or by editing the
default values at the top of the script.

### 8.9 Ion Concentration Calculator

**Script:** `templates/ions.R` → copy to project and edit.  Or use
`CSS/MD/R/ions.R` as a reference implementation.

A simple R script (base R only, no external packages) that calculates
the number of Na⁺, Mg²⁺, and Cl⁻ ions required to reach target ionic
concentrations for a given water-box size and system charge.

Two methods are computed:
- **SPLIT** — classic approximation `N0 = Nw × conc / 55.5`
- **SLTCAP** — self-consistent solution (Machado + SLTCAP, see
  https://doi.org/10.1021/acs.jctc.7b01254).  This is the preferred
  method; it does not require `N0 >> |Q|`.

**Usage:**
```r
source("ions.R")
```

**Inputs to edit:** `C_NaCl`, `C_MgCl2`, `Nw` (water molecules from
leap.log), `Q` (net system charge).

**Output:** Prints the recommended `leap.sh -n <na> -m <mg> -c <cl>`
command directly to the console.

### 8.10 Required Files Summary

| Step | Script | Input files (path pattern) | Output files |
|------|--------|---------------------------|--------------|
| 0 | (external) | — | `QM/<mol>.xyz` (or `.mol2`) |
| 1 | `orca_steps_wsl.sh` | `QM/<mol>.xyz` | `SP_gas.molden`, `SP_solv.molden` |
| 2 | `multiwfn_steps_*.ps1/.sh` | `SP_gas.molden`, `SP_solv.molden` | `QM/<mol>_opt.chg` |
| 3 | `ligand_prep.sh` | `QM/<mol>.mol2` (or `.xyz`), `QM/<mol>_opt.chg` | `ff/<mol>_resp.*` (prepin, frcmod, pdb) |
| 4 | `model_prep.py` | `../<SEQ>_<LGD><suffix>/<JOB>/boltz_results_input/predictions/input/*.cif`, `ff/<mol>_resp.pdb` (not required when `--lgd free`) | `pdb_for_md/<PREFIX>_<LIGAND>_constrained/input_model_{i}.pdb` |
| 5 | `leap.sh` | `ff/<lig>_resp.prepin` (not required with `-l free`), `ff/<lig>_resp.frcmod` (not required with `-l free`), PDB models from Step 4 | `leap/<PREFIX>_<LIGAND>_constrained/cplx_{i}.prmtop/.rst7/.pdb` |

### 8.11 Dependencies Table

| Script | Platform | Dependencies |
|--------|----------|-------------|
| `orca_steps_wsl.sh` | WSL/Linux | ORCA (≥ 6.x), orca_2mkl |
| `multiwfn_steps_windows.ps1` | Windows | Multiwfn, PowerShell |
| `multiwfn_steps_linux.sh` | Linux | Multiwfn (Linux binary), Bash |
| `ligand_prep.sh` | WSL/Linux | AmberTools (antechamber, prepgen, parmchk2); optional: OpenBabel (for .xyz → .mol2 fallback) |
| `model_prep.py` | Any (PyMOL) | PyMOL (pymol-open-source preferred) |
| `leap.sh` | WSL/Linux | AmberTools (tleap, ambpdb) |
| `ions.R` | Any (R) | Base R only (no external packages) |

### 8.12 Troubleshooting

| Issue | Likely Cause | Solution |
|-------|-------------|----------|
| `SP_gas.molden` not found | WSL→Windows copy not performed | Copy Molden files from WSL to `<MoldenDir>` before running the PS1 script |
| Multiwfn fails silently (Windows) | Batch wrapper issue | Check host output for errors; verify Multiwfn path in the PS1 script |
| Multiwfn fails (Linux) | Missing input or wrong path | Ensure `$Multiwfn` path is correct and the molden files exist |
| `mol2` file not found in `QM/` | Wrong filename or directory | Name the file `<MOL>.mol2` and place in `QM/`; or provide `<MOL>.xyz` as fallback |
| `.xyz` → `.mol2` conversion fails | Missing converter | Install `antechamber` (AmberTools) or `obabel` (OpenBabel) |
| RESP2 charges look wrong | Atom ordering mismatch | Verify identical atom ordering between `.mol2` and `.chg` files |
| PyMOL: `prep_model()` fails | Missing job directory | Check that the Boltz-2 prediction CIFs are in the expected path; use `--job` explicitly |
| PyMOL: modifications don't match system | Default DNA-specific operations | Edit the "EDIT THIS BLOCK" section in `model_prep.py` for your system (RNA/protein/other) |
| `leap.sh`: `CSS1_` directory not found | Wrong prefix | Use `-p` to specify your project prefix (e.g. `-p MYPROJ`) |
| tLEAP reports missing parameters | frcmod has ATTN warnings | Review `ff/<ligand>_resp.frcmod` and adjust parmchk2 output |
| OPC water model unrecognised | AmberTools version | Requires AmberTools ≥ 18 for `leaprc.water.opc` |
| Ion counts are wrong | Incorrect Nw or Q | Run `leap.sh` once with approximate values, extract Nw from leap.log, then use `ions.R` to recalculate |

---

## 9. Amber MD Preparation and Minimisation Pipeline

In addition to the Boltz-2 prediction and ligand-parameterisation workflows,
the repository provides a SLURM-based **Amber MD preparation and minimisation
pipeline** (`PrepAndMin.slurm`).  Two companion generator scripts configure
the pipeline for each project:
`generate_pmemd_inputs.py` (§9.5) generates the Amber input control files,
and `generate_md_slurm.py` (§9.2) generates the SLURM array job script.
The pipeline takes the solvated Amber systems from `leap.sh` (§8.8) and runs
an 11-step protocol (minimisation + equilibration + 100 ns production + final
minimisation).

### 9.1 Pipeline Overview

```
generate_pmemd_inputs.py                          (§9.5)
  │
  │  Generates all 11 pmemd/in/*.in files from
  │  --macromol-type, --macromol-start/end, --ligand-idx, --production-ns
  │  Output:  pmemd/in/step{1-10}.in, final_min.in
  ▼
generate_md_slurm.py                               (§9.2)
  │
  │  Reads PrepAndMin.slurm template, replaces placeholders
  │  Output:  slurm/<EXPERIMENT>_PrepAndMin.slurm
  ▼
leap.sh                                          (§8.8)
  │
  │  Deliverables:
  │    leap/<EXPERIMENT>/cplx_{N}.prmtop
  │    leap/<EXPERIMENT>/cplx_{N}.rst7
  ▼
PrepAndMin.slurm  (SLURM array job — submit via sbatch)
  │
  │  11-step protocol on compute-node temporary storage
  │  Step  1:  pmemd     — solvent minimisation, restrained
  │  Step  2:  pmemd.cuda— NVT heating, restrained              GPU
  │  Step  3:  pmemd     — solute minimisation, medium restraints
  │  Step  4:  pmemd     — solute minimisation, weak restraints
  │  Step  5:  pmemd     — unrestrained minimisation
  │  Step  6:  pmemd.cuda— NPT equilibration, mild restraints   GPU
  │  Step  7:  pmemd.cuda— NPT equilibration, reduced restraints GPU
  │  Step  8:  pmemd.cuda— NPT, backbone-only restraints        GPU
  │  Step  9:  pmemd.cuda— short unrestrained NPT, 2 fs dt      GPU
  │  Step 10:  pmemd.cuda— 100 ns production NPT                GPU
  │  Step 11:  pmemd.cuda— final minimisation                   GPU
  │
  ▼
Output directories:
  pmemd/out/<EXPERIMENT>/J<jobid>/task_<N>/
  ├── prep/
  │   ├── step1{N}.out, .ncrst, .nc, .info      — per-step files
  │   ├── step2{N}.out, .ncrst, .nc, .info
  │   ├── ...
  │   ├── pdb/
  │   │   ├── step1{N}.pdb
  │   │   └── ...
  │   ├── traj/
  │   │   ├── step10{N}_stripped.prmtop          — stripped topology
  │   │   ├── step10{N}_stripped.nc              — stripped trajectory
  │   │   ├── step10_2drmsd{N}.dat               — 2D RMSD matrix
  │   │   └── step10_atomicfluct{N}.dat          — RMSF by residue
  │   └── cpptraj/
  │       ├── step10_stripped{N}.cpptraj
  │       └── step10_analysis{N}.cpptraj
  ├── final_min/
  │   ├── final_min{N}.out, .ncrst, .info
  │   ├── pdb/
  │   │   ├── final_min{N}_cpptrajed.pdb
  │   │   └── final_min{N}_stripped.pdb
  │   └── cpptraj/
  │       ├── final_min_stripped{N}.cpptraj
  │       └── final_min{N}.cpptraj
  └── run.info                                    — job metadata
```

### 9.2 SLURM Script Generation

**Script:** `templates/MD/python/generate_md_slurm.py` → copy to project (`cp -r templates/MD/*`)
and run from the `MD/` directory.

Reads the `slurm/PrepAndMin.slurm` template, replaces `__PLACEHOLDER__` markers with
your project parameters, and writes a configured `.slurm` file.

**Usage:**

```bash
cd my_project/MD
python python/generate_md_slurm.py \
    --experiment CSS1_HCY_constrained \
    --array 1,2,3,5,7 \
    --scratchdir scratch/boltz/projects/CSS/MD
```

For a **no-ligand** run, use `"free"` in the experiment name:

```bash
python python/generate_md_slurm.py \
    --experiment CSS1_free_constrained \
    --array 0-4 \
    --scratchdir scratch/boltz/projects/CSS/MD
```

**Output:**  `slurm/<experiment>_PrepAndMin.slurm` (e.g. `slurm/CSS1_HCY_constrained_PrepAndMin.slurm`)

The output filename includes the experiment name so multiple projects can
coexist in the same `slurm/` directory without ambiguity.

**Arguments:**

| Arg | Required | Default | Description |
|-----|----------|---------|-------------|
| `--experiment` | yes | — | Experiment name; sets `-J`, `EXPERIMENT=`, and output filename |
| `--array` | no | `0-4` | SLURM array task IDs (e.g. `0-4`, `1,2,3,5,7`) |
| `--scratchdir` | yes | — | Relative path under `$SLURM_SUBMIT_DIR` to the scratch directory containing `pmemd/` and `leap/` |
| `--gpu` | no | `gpu-h100,gpu-l40` | GPU partition(s) for `#SBATCH -p` |
| `--residues` | no | `1-46` | Residue mask for cpptraj `align` (e.g. `1-46`) |
| `--solvent` | no | `WAT` | Solvent residue name for cpptraj `strip` (e.g. `WAT`) |
| `--sequential` | no | `false` | If set, limits array concurrency to 1 (`%1` throttle). Tasks run one at a time |

After generation, transfer the output `.slurm` file to the cluster together
with the `pmemd/in/` input files, then submit:

```bash
sbatch slurm/CSS1_HCY_constrained_PrepAndMin.slurm
```

### 9.3 The 11-Step Protocol

| Step | Engine | Restraints | Length | dt | Purpose |
|------|--------|-----------|--------|----|---------|
| 1 | `pmemd` (CPU) | Heavy atoms, 5.0 kcal/mol·A² | 1000 steps SD | — | Relax solvent/ions |
| 2 | `pmemd.cuda` (GPU) | Heavy atoms, 5.0 kcal/mol·A² | 15000 steps | 1 fs | NVT heating to 300 K |
| 3 | `pmemd` (CPU) | Heavy atoms, 2.0 kcal/mol·A² | 1000 steps SD | — | Minimise solute |
| 4 | `pmemd` (CPU) | Heavy atoms, 0.1 kcal/mol·A² | 1000 steps SD | — | Minimise solute (weak) |
| 5 | `pmemd` (CPU) | **None** | 1000 steps SD | — | Unrestrained minimisation |
| 6 | `pmemd.cuda` (GPU) | Heavy atoms, 1.0 kcal/mol·A² | 5000 steps | 1 fs | First NPT equilibration |
| 7 | `pmemd.cuda` (GPU) | Heavy atoms, 0.5 kcal/mol·A² | 5000 steps | 1 fs | NPT, weaker restraints |
| 8 | `pmemd.cuda` (GPU) | Backbone only, 0.5 kcal/mol·A² | 10000 steps | 1 fs | NPT, backbone-only |
| 9 | `pmemd.cuda` (GPU) | **None** | 5000 steps | 2 fs | Unrestrained NPT (warm-up) |
| 10 | `pmemd.cuda` (GPU) | **None** | 50000000 steps | 2 fs | **100 ns production** |
| 11 | `pmemd.cuda` (GPU) | **None** | 20000 steps SD+CG | — | Final minimisation |

### 9.4 Input File Descriptions

All Amber input control files live in `pmemd/in/` within the project directory.

> **These files are now generated by `generate_pmemd_inputs.py` (§9.5).** Run the
> generator with your system parameters instead of editing manually.  The
> generator handles residue numbering, backbone atom selection, and production
> length.  If you still need to hand-edit, note that each file must end with a
> blank line (`/`) to avoid Amber "End of line" errors.

| File | Step | Key parameters to adjust for a new system |
|------|------|-------------------------------------------|
| `step1.in` | 1 | `:N_MASK` (solute residue range), `restraint_wt` |
| `step2.in` | 2 | `:N_MASK`, `restraint_wt`, `nstlim` (heating length), `dt` |
| `step3.in` | 3 | `:N_MASK`, `restraint_wt` (medium) |
| `step4.in` | 4 | `:N_MASK`, `restraint_wt` (weak) |
| `step5.in` | 5 | (none — unrestrained) |
| `step6.in` | 6 | `:N_MASK`, `restraint_wt` (mild), first NPT block |
| `step7.in` | 7 | `:N_MASK`, `restraint_wt` (reduced) |
| `step8.in` | 8 | `:N_DNA_BB` (backbone atom mask), `N_BB_ATOMS` list |
| `step9.in` | 9 | (none — unrestrained, 2 fs dt transition) |
| `step10.in` | 10 | `nstlim` (production length), `cut`, output frequencies |
| `final_min.in` | 11 | `maxcyc`, `ncyc` (minimisation depth) |

### 9.5 Input File Generation

**Script:** `templates/MD/python/generate_pmemd_inputs.py` → copy to project and run
from the `MD/` directory.

Generates all 11 `pmemd/in/*.in` files from your system parameters, removing
the need for manual editing.  Existing files in `pmemd/in/` are automatically
moved to `pmemd/in/archive/` before generation.

**Usage:**

```bash
cd my_project/MD
python python/generate_pmemd_inputs.py \
    --macromol-type dna \
    --macromol-start 1 \
    --macromol-end 46 \
    --ligand-idx 47 \
    --production-ns 100
```

For a **no-ligand** system, simply omit `--ligand-idx`:

```bash
cd my_project/MD
python python/generate_pmemd_inputs.py \
    --macromol-type dna \
    --macromol-start 1 \
    --macromol-end 46 \
    --production-ns 100
```

The output is written directly to `pmemd/in/` in the current working directory.

**Arguments:**

| Arg | Required | Default | Description |
|-----|----------|---------|-------------|
| `--macromol-type` | no | `dna` | `dna`, `rna`, or `protein` — determines backbone atom set |
| `--macromol-start` | no | `1` | First residue of the macromolecule (used for backbone mask) |
| `--macromol-end` | no | `46` | Last residue of the macromolecule |
| `--ligand-idx` | no | — | Ligand residue index (omit if no ligand present) |
| `--production-ns` | no | `100` | Production simulation length in nanoseconds |

**What changes per file:**

| File | Generated from parameters |
|------|---------------------------|
| `step1.in` – `step4.in` | Solute restraint mask (`:1-<MAX_END>&!@H=`) — `MAX_END` = max(`macromol_end`, `ligand_idx`) |
| `step5.in` | Static (no restraints) |
| `step6.in`, `step7.in` | Same solute mask as steps 1-4 |
| `step8.in` | Backbone mask (`:START-END@ATOMS`) + optional ligand clause (`|(:LIG&!@H=)`). Backbone atoms depend on `macromol_type` |
| `step9.in` | Static (no restraints, 2 fs) |
| `step10.in` | `nstlim` computed as `production_ns * 1000 / 0.002` |
| `final_min.in` | Static (no restraints, 20000 steps) |

**Backbone atom sets by molecule type:**

| `--macromol-type` | Amber atom selection |
|-------------------|---------------------|
| `dna` | `@P,OP1,OP2,O5',C5',C4',C3',O3'` |
| `rna` | `@P,OP1,OP2,O5',C5',C4',C3',O3',O2'` |
| `protein` | `@N,CA,C,O` |

### 9.6 Output Files (per task)

| File(s) | Contents |
|---------|----------|
| `prep/step{N}_{TASK}_J{JOBID}.out` | Amber log file for step N |
| `prep/step{N}_{TASK}_J{JOBID}.ncrst` | Amber restart file (binary NetCDF) |
| `prep/step{N}_{TASK}_J{JOBID}.nc` | Trajectory file (NetCDF) — written for steps 2, 6-10 |
| `prep/step{N}_{TASK}_J{JOBID}.info` | Amber mdinfo progress log |
| `prep/pdb/step{N}_{TASK}_J{JOBID}.pdb` | PDB snapshot (converted with ambpdb) |
| `final_min/final_min_{TASK}_J{JOBID}.out` | Amber log for final minimisation |
| `final_min/final_min_{TASK}_J{JOBID}.ncrst` | Minimised restart file |
| `final_min/final_min_{TASK}_J{JOBID}_cpptrajed.ncrst` | Autoimaged, aligned restart (cpptraj) |
| `final_min/final_min_{TASK}_J{JOBID}_stripped.ncrst` | Autoimaged, aligned, solvent-stripped restart |
| `final_min/final_min_{TASK}_J{JOBID}_stripped.prmtop` | Solvent-stripped topology |
| `final_min/pdb/final_min_{TASK}_J{JOBID}_cpptrajed.pdb` | PDB of full minimised structure |
| `final_min/pdb/final_min_{TASK}_J{JOBID}_stripped.pdb` | PDB of solvent-stripped structure |
| `run.info` | `SLURM_JOBID TASK_ID hostname` metadata |

### 9.7 SLURM Configuration

The generated `.slurm` script uses **SLURM array jobs** so that each
complex (`cplx_N`) runs as an independent task.  The `--array` argument
to `generate_md_slurm.py` (§9.2) sets `#SBATCH --array` to match the
`cplx_N` files produced by `leap.sh`:

```bash
#SBATCH --array=0-4       # task IDs — one per complex
```

Similarly, `--gpu` sets `#SBATCH -p` and `--experiment` sets `#SBATCH -J`.
Use `--sequential` to append `%1`, limiting the array to one concurrent task
at a time (useful when GPU memory or node resources are shared).

**Restart protection:** The loop checks for existing `.ncrst` files before
each step.  If a `.ncrst` is found, that step is skipped.  This allows
safe resubmission after a time-limit or transient failure.

**Temporary storage:** All MD steps run on compute-node local storage
(`$TMPDIR`), avoiding network filesystem latency.  Output is copied to
scratch only after successful completion.

**Scratch directory layout (on cluster):**

```
scratch/boltz/projects/<PROJECT>/MD/
  pmemd/
    in/                  — copied to each job (step*.in, final_min.in)
    out/
      <EXPERIMENT>/
        J<jobid>/
          task_<N>/
            prep/
            final_min/
            run.info
  leap/
    <EXPERIMENT>/
      cplx_<N>.prmtop    — sourced by the job
      cplx_<N>.rst7      — sourced by the job
```

### 9.8 Customisation Guide

To adapt this pipeline for a **new system** (different DNA sequence,
different ligand, protein, RNA, etc.):

1. **Copy the template files** (see §2.2).
2. **Generate pmemd input files** with `generate_pmemd_inputs.py`:
   ```bash
   # With ligand
   python python/generate_pmemd_inputs.py \
       --macromol-type dna \
       --macromol-start 1 \
       --macromol-end 46 \
       --ligand-idx 47 \
       --production-ns 100

   # Without ligand (omit --ligand-idx)
   python python/generate_pmemd_inputs.py \
       --macromol-type dna \
       --macromol-start 1 \
       --macromol-end 46 \
       --production-ns 100
   ```
   See §9.5 for full argument details.
3. **Generate the SLURM script** with `generate_md_slurm.py`:
   ```bash
   python python/generate_md_slurm.py \
       --experiment <EXPERIMENT> \
       --array <TASK_IDS> \
       --scratchdir <SCRATCH_DIR> \
       --residues <N_DNA> \
       --solvent <SOLVENT>
   ```
   See §9.2 for full argument details.
4. **Transfer** the generated `.slurm` file, `pmemd/in/` input files, and
   `leap/<EXPERIMENT>/` outputs to the cluster.
5. **Submit**:
   ```bash
   sbatch slurm/<EXPERIMENT>_PrepAndMin.slurm
   ```

### 9.9 Troubleshooting

| Issue | Likely Cause | Solution |
|-------|-------------|----------|
| `pmemd: command not found` | AmberTools not loaded | `module load amber` |
| `pmemd.cuda: command not found` | Amber GPU module not available | Use `module load amber/22` or check cluster docs |
| Step fails with `DISANG` error | `nmropt=1` but no restraint file | Set `nmropt=0` in all `.in` files (restraints off) |
| Step fails with `IMAGE NOT PRESENT` | Bad initial coordinates | Check leap-generated `.rst7` with `ambpdb` |
| `cpptraj: command not found` | AmberTools cpptraj not in PATH | `module load amber/22` |
| `ambpdb: command not found` | AmberTools not loaded | Use `cpptraj` instead if ambpdb unavailable |
| `segmentation fault` in pmemd.cuda | GPU out of memory | Reduce system size, increase `ntasks-per-node` or use `pmemd` (CPU) |
| Job hits time limit mid-pipeline | `#SBATCH -t` too short | Increase wall time; resubmit — restart protection will skip completed steps |
| Restart protection skips wrong steps | Filename mismatch | Check `SUFFIX` construction; `_${TASK}_J${JOBID}` must match existing `.ncrst` files |

### 9.10 MD Trajectory Analysis

**Script:** `templates/MD/R/md_utils.R` → copy to `{project}/MD/R/` (included
in the `cp -r templates/MD/*` copy step in §2.2).

Provides R functions for reading and plotting MD trajectory outputs
generated by the SLURM pipeline (§9).  Requires `data.table`, `stringr`,
`ggplot2`, and `ggrepel`.

#### Functions

| Function | Description |
|----------|-------------|
| `read_rmsd(file, simulation_time, extract_frequency, time_step)` | Reads a 2D RMSD `.dat` file (from cpptraj `2drms`) into a long-format `data.table` with columns `t1`, `t2`, `rmsd`, `model`, `experiment`, `seq`, `ligand`. Optionally subsamples via `time_step`. |
| `read_atomicfluct(file)` | Reads an atomic fluctuation `.dat` file (from cpptraj `atomicfluct`) into a `data.table` with columns `res`, `atomfluct`, `experiment`, `model`, `seq`, `ligand`. |
| `plot_rmsd(rmsd, max_t, time_step, scale)` | Generates a 2D RMSD heatmap (`geom_raster`) faceted by experiment and model, using a viridis colour scale. |
| `plot_atomfluct(atomfluct, res_cutoff, label_cutoff, scale)` | Plots per-residue atomic fluctuation as a line/point chart with automatic labelling of high-fluctuation regions via `ggrepel`. |
| `theme_custom(scaling)` | Returns a publication-quality `ggplot2` theme with configurable font sizes. |

#### Path Convention

Both `read_rmsd()` and `read_atomicfluct()` extract experiment and model
metadata from the file path, which must follow the pattern:

```
<project>/<category>/<method>/out/<experiment>_<replicate>/<model>/<file>.dat
```

- Path component `[5]` → experiment (e.g. `CSS1_rep1` → `CSS1·rep1`)
- Path component `[7]` → model name (e.g. `wildtype`)

#### Usage

```r
source("md_utils.R")

# Read all RMSD and atomic fluctuation files
rmsd_files  <- list.files("pmemd/out", pattern = "2drmsd\\.dat$", full.names = TRUE, recursive = TRUE)
fluct_files <- list.files("pmemd/out", pattern = "atomicfluct\\.dat$", full.names = TRUE, recursive = TRUE)

rmsd_dt  <- lapply(rmsd_files, read_rmsd) |> rbindlist()
fluct_dt <- lapply(fluct_files, read_atomicfluct) |> rbindlist()

# Plot with optional time subsampling
plot_rmsd(rmsd_dt, max_t = 100, time_step = 1, scale = 0.7)
plot_atomfluct(fluct_dt, res_cutoff = 46, label_cutoff = 4)
```

**Note:** `read_rmsd()` requires the global variables `sim_time` (total
simulation time in ns) and `extract_freq` (frame extraction interval) to
be defined before calling, or passed explicitly via `simulation_time` and
`extract_frequency`.

### 9.11 PyMOL Visualisation Utilities

**Script:** `templates/MD/python/pymol_utils.py` → copy to
`{project}/MD/python/` (included in the `cp -r templates/MD/*` copy step
in §2.2).

Provides reusable PyMOL functions for loading, aligning, and colouring
structures produced by the Amber MD pipeline.  Designed to be imported
by project-specific entry-point scripts (e.g. `final_min_pymol.py`).

#### Functions

| Function | Description |
|----------|-------------|
| `process_final_min(project, experiment, job)` | Load stripped PDB files from the final-minimisation stage, align all models, apply base/element colouring, and save a `.pse` session.  Then reset nucleic acids to white and overlay constraint-coloured residues for a second `.pse` session. |
| `load_constraints(project, experiment)` | Parse contact constraints from the experiment's YAML file.  Returns deduplicated residue-pair tuples or `None` if the experiment is not constrained. |
| `color_constrained_residues(constraints)` | Build a graph from constraint pairs, find connected components via BFS, and assign each component a unique colour from the `CONSTRAINT_COLORS` palette.  Automatically falls back to residue-number-only selection when chain IDs are absent in the PDB. |

#### Usage

```python
# Import all three functions
from pymol_utils import (
    process_final_min,
    load_constraints,
    color_constrained_residues,
)

# Process a single experiment/job pair
process_final_min("my_project", "Seq1_L01_constrained", "J1234567")
```

Or as a standalone script (see `CSS/MD/python/final_min_pymol.py` for a
complete example):

```python
import os
from pymol_utils import process_final_min

output_root = "my_project/MD/pmemd/out"
for exp in sorted(os.listdir(output_root)):
    for job in os.listdir(os.path.join(output_root, exp)):
        process_final_min("my_project", exp, job)
```

#### Chain-ID Fallback

PDB files produced by Amber's `ambpdb` often have blank chain IDs, while
the YAML constraints reference chain `"A"`.  `color_constrained_residues()`
detects this automatically: if the referenced chains are not found among the
loaded objects, the chain filter is dropped and residues are selected by
residue number only.

---

## 10. File Reference

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
| **MD parameterisation templates** | |
| `templates/MD/orca_steps_wsl.sh` | ORCA paths, molecule name, charge, multiplicity, solvent |
| `templates/MD/multiwfn_steps_windows.ps1` | Multiwfn path, Molden directory, output directory |
| `templates/MD/multiwfn_steps_linux.sh` | Multiwfn path, Molden directory, output directory |
| `templates/MD/ligand_prep.sh` | `MOL`, `NC`, `S` (molecule name, net charge, verbosity) |
| `templates/MD/python/model_prep.py` | `seq`, `lgd`, `lgd_file` (omit when `lgd=free`), model range; PyMOL modification block (EDIT THIS BLOCK) |
| `templates/MD/leap.sh` | Default `-p` prefix, `-l` ligand name (use `-l free` for no-ligand mode), model range, ion counts |
| `templates/MD/R/ions.R` | `C_NaCl`, `C_MgCl2`, `Nw`, `Q` (ion concentrations, water count, system charge) |
| `templates/MD/R/md_utils.R` | `sim_time`, `extract_freq` (globals for `read_rmsd()`); all other parameters are function arguments (§9.10) |
| **MD preparation + minimisation templates** | |
| `templates/MD/python/generate_pmemd_inputs.py` | `macromol_type`, `macromol_start`, `macromol_end`, `ligand_idx`, `production_ns` — generates all 11 `pmemd/in/*.in` files from system parameters (§9.5) |
| `templates/MD/pmemd/in/step1.in` … `step10.in` | Now **generated** by `generate_pmemd_inputs.py` (§9.5); no manual editing needed. The `archive/` directory holds previous versions |
| `templates/MD/pmemd/in/final_min.in` | Now **generated** by `generate_pmemd_inputs.py` (§9.5) |
| `templates/MD/python/generate_md_slurm.py` | `experiment`, `array`, `scratchdir`, `gpu`, `residues`, `solvent` — generates a configured `.slurm` from `PrepAndMin.slurm` (§9.2) |
| `templates/MD/python/pymol_utils.py` | No edits required — reusable functions imported by project scripts (§9.11) |
| `templates/MD/slurm/PrepAndMin.slurm` | Template with `__PLACEHOLDER__` markers; used by `generate_md_slurm.py` (§9.1) |

---

## 11. Appendix: Contact Constraints Format

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
