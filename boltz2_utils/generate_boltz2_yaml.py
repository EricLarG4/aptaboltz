#import os
import yaml
from pathlib import Path
from typing import Optional, Union


# ── Type aliases ──────────────────────────────────────────────────────────────

ResiduePair = Union[
    tuple[int, int],               # (res_i, res_j) — both on dna_chain (backward-compat)
    tuple[str, int, str, int],     # (chain_i, res_i, chain_j, res_j) — explicit chains
]

AdditionalPairs = Union[
    list[ResiduePair],                              # Flat list — same pairs for all entries
    dict[str, list[ResiduePair]],                    # Per-sequence: {seq_name: pairs}
    dict[str, dict[str, list[ResiduePair]]],         # Per-seq + per-ligand: {seq: {lig: pairs}}
]


# ── Custom YAML dumper: 4-space indent, indented block sequences ───────────────


class InlineList(list):
    """Marker class so the YAML dumper renders this list inline."""
    pass


class BoltzDumper(yaml.Dumper):
    """
    Custom YAML dumper:
    - Block sequences are indented under their mapping key (4-space indent).
    - Block mappings use a 2-space indent, so sequence items render as
      '    - key:' (1 space after dash) rather than '    -   key:'.
    """

    def increase_indent(self, flow=False, indentless=False):
        if flow:
            return super().increase_indent(flow=True, indentless=indentless)
        if isinstance(self.event, yaml.SequenceStartEvent):
            # Sequences: 4-space indent, never indentless
            self.indents.append(self.indent)
            self.indent = (
                self.best_indent
                if self.indent is None
                else self.indent + self.best_indent
            )
        else:
            # Mappings: 2-space indent
            self.indents.append(self.indent)
            self.indent = 0 if self.indent is None else self.indent + 2


def inline_list_representer(dumper, data):
    return dumper.represent_sequence("tag:yaml.org,2002:seq", data, flow_style=True)


BoltzDumper.add_representer(InlineList, inline_list_representer)


# ── Core builder ──────────────────────────────────────────────────────────────


def build_stem_contacts(
    chain_id: str,
    sequence: str,
    stem_length: int = 8,
    max_distance: float = 5.0,
) -> list[dict]:
    """
    Generate Watson-Crick stem contact constraints for the first and last
    `stem_length` residues of *sequence* (1-based indexing).
    """
    seq_len = len(sequence)
    contacts = []
    for i in range(1, stem_length + 1):
        j = seq_len - i + 1
        contacts.append(
            {
                "contact": {
                    "token1": InlineList([chain_id, i]),
                    "token2": InlineList([chain_id, j]),
                    "max_distance": max_distance,
                    "force": True,
                }
            }
        )
    return contacts


def build_additional_contacts(
    chain_id: str,
    additional_pairs: list[ResiduePair],
    max_distance: float = 5.0,
) -> list[dict]:
    """
    Convert user-specified residue pairs into Boltz-2 contact constraint dicts.

    Each pair in *additional_pairs* can be one of two formats:

    **2-tuple (backward-compatible)**
        ``(res_i, res_j)``
        Both residues are assigned to *chain_id* (the dna_chain).  Example::

            (5, 10)          # residues 5 and 10 on chain "A"

    **4-tuple (explicit chains)**
        ``(chain_i, res_i, chain_j, res_j)``
        Each residue specifies its own chain identifier.  This enables
        cross-chain contacts (e.g. DNA residue to ligand residue)::

            ("A", 5, "B", 10)   # chain A residue 5  <->  chain B residue 10

    Every pair is converted to a force-applied contact constraint.
    """
    contacts = []
    for pair in additional_pairs:
        if len(pair) == 2:
            res_i, res_j = pair
            chain_i = chain_j = chain_id
        else:
            chain_i, res_i, chain_j, res_j = pair

        contacts.append(
            {
                "contact": {
                    "token1": InlineList([chain_i, res_i]),
                    "token2": InlineList([chain_j, res_j]),
                    "max_distance": max_distance,
                    "force": True,
                }
            }
        )
    return contacts


def build_boltz2_input(
    seq_name: str,
    sequence: str,
    ligand_name: str,
    ccd: Optional[str] = None,
    smiles: Optional[str] = None,
    chain_id: str = "A",
    ligand_chain: str = "B",
    cyclic: bool = False,
    stem_length: int = 8,
    additional_pairs: Optional[list[ResiduePair]] = None,
    max_distance: float = 5.0,
    msa: Optional[str] = None,
    molecule_type: str = "dna",
) -> dict:
    """
    Build the Boltz-2 input dictionary for one (sequence, ligand) pair.

    Priority: *ccd* is used when non-empty; otherwise *smiles* is used.
    Raises ValueError when neither is provided.

    molecule_type : ``"dna"``, ``"rna"``, or ``"protein"``.
                    Controls the sequence-type key in the YAML output.
                    Stem contacts are only generated for ``"dna"`` and
                    ``"rna"`` (which have Watson-Crick base pairing).
    additional_pairs : optional list of residue pairs — see
                       :func:`build_additional_contacts` for the accepted formats.
                       Each pair is added as a force-applied contact constraint
                       on the chain, in addition to any stem contacts derived
                       from *stem_length*.
    msa : optional path to an MSA file; when provided it is added as an
          ``msa`` key inside the sequence block.
    """
    if not ccd and not smiles and ligand_name.lower() != "free":
        raise ValueError(
            f"[{seq_name} / {ligand_name}] Either 'ccd' or 'smiles' must be provided."
        )

    # Determine the sequence-type key for the YAML
    valid_types = ("dna", "rna", "protein")
    if molecule_type not in valid_types:
        raise ValueError(
            f"molecule_type must be one of {valid_types}, got '{molecule_type}'"
        )

    # Sequence entry
    seq_body: dict = {
        "id": chain_id,
        "sequence": sequence,
        "cyclic": cyclic if molecule_type in ("dna", "rna") else False,
    }
    if msa is not None:
        seq_body["msa"] = msa

    seq_entry = {molecule_type: seq_body}

    # Contact constraints — stem contacts only for nucleic acids
    contacts = []
    if molecule_type in ("dna", "rna"):
        contacts = build_stem_contacts(chain_id, sequence, stem_length, max_distance)

    # Merge any user-specified additional contacts
    if additional_pairs:
        contacts.extend(
            build_additional_contacts(chain_id, additional_pairs, max_distance)
        )

    # if ligand is equal to "free" (case-insensitive), omit the ligand block entirely
    if ligand_name.lower() == "free":
        result: dict = {"sequences": [seq_entry]}
    else:
        # Ligand entry — CCD takes priority over SMILES
        ligand_body: dict = {"id": ligand_chain}
        if ccd:
            ligand_body["ccd"] = ccd
        else:
            ligand_body["smiles"] = smiles

        ligand_entry = {"ligand": ligand_body}

        result: dict = {"sequences": [seq_entry, ligand_entry]}

    if contacts:
        result["constraints"] = contacts
    return result


# ── Pair resolution helper ────────────────────────────────────────────────────


def _resolve_pairs(
    additional_pairs: AdditionalPairs | None,
    seq_name: str,
    ligand_name: str,
) -> list[ResiduePair]:
    """
    Resolve the *additional_pairs* specification for a given (seq, ligand) pair.

    Resolution logic
    -----------------
    1. ``None`` or a flat ``list[ResiduePair]`` → returned as-is (or empty list).
    2. ``dict[str, list[ResiduePair]]`` → lookup by *seq_name*.
    3. ``dict[str, dict[str, list[ResiduePair]]]`` → lookup
       ``[seq_name][ligand_name]``; if the ligand is not listed, fall back to
       ``[seq_name]["*"]`` (wildcard) or empty list.

    Examples
    ---------
    **Flat list** (backward-compatible)::

        additional_pairs = [(5, 10), (1, 20)]

    **Per-sequence** (same pairs for all ligands, including "free")::

        additional_pairs = {
            "CSS1": [(1, 10), (2, 9)],
            "CSS2": [("A", 3, "B", 8)],
        }

    **Per-sequence + per-ligand** (different pairs for different conditions)::

        additional_pairs = {
            "CSS1": {
                "free": [("A", 1, "A", 10)],        # only for free (no ligand)
                "C0R":  [("A", 5, "B", 10)],         # only for C0R
                "*":    [("A", 2, "A", 9)],           # fallback for all other ligands
            },
        }
    """
    if additional_pairs is None:
        return []

    # Case 1: flat list → return as-is
    if isinstance(additional_pairs, list):
        return additional_pairs

    # Case 2: per-sequence dict → {seq_name: list[ResiduePair]}
    if all(isinstance(v, list) for v in additional_pairs.values()):
        return additional_pairs.get(seq_name, [])

    # Case 3: per-sequence + per-ligand dict → {seq: {lig: list[ResiduePair]}}
    seq_dict = additional_pairs  # dict[str, dict[str, list[ResiduePair]]]
    if seq_name in seq_dict:
        lig_dict = seq_dict[seq_name]
        # Exact match
        if ligand_name in lig_dict:
            return lig_dict[ligand_name]
        # Wildcard fallback
        if "*" in lig_dict:
            return lig_dict["*"]
    return []


# ── File writer ───────────────────────────────────────────────────────────────


def generate_yaml_file(
    seq_name: str,
    sequence: str,
    ligand_name: str,
    ccd: Optional[str] = None,
    smiles: Optional[str] = None,
    output_dir: str = "yaml",
    stem_length: int = 8,
    max_distance: float = 5.0,
    msa: Optional[str] = None,
    additional_pairs: Optional[list[ResiduePair]] = None,
    molecule_type: str = "dna",
    **kwargs,
) -> Path:
    """
    Generate a single Boltz-2 YAML file and return its path.

    molecule_type    : ``"dna"``, ``"rna"``, or ``"protein"``.
                       Controls the sequence-type key in the YAML output.
    stem_length      : number of Watson-Crick base-pair constraints to add at each
                       end of the stem (set to 0 to omit the constraints block).
    max_distance     : upper-bound distance (Å) for each stem contact.
    msa              : optional path to an MSA file; when provided it is written as
                       an ``msa`` key inside the sequence block.
    additional_pairs : optional list of residue pairs — see
                       :func:`build_additional_contacts` for the accepted formats.

    Filename pattern:
        <seq_name>_<ligand_name>.yaml                                     (no extra features)
        <seq_name>_<ligand_name>_custom_msa.yaml                           (with MSA)
        <seq_name>_<ligand_name>_constrained.yaml                          (with additional_pairs)
        <seq_name>_<ligand_name>_custom_msa_constrained.yaml               (with both MSA and additional_pairs)
    """
    data = build_boltz2_input(
        seq_name=seq_name,
        sequence=sequence,
        ligand_name=ligand_name,
        ccd=ccd,
        smiles=smiles,
        stem_length=stem_length,
        max_distance=max_distance,
        msa=msa,
        additional_pairs=additional_pairs,
        molecule_type=molecule_type,
        **kwargs,
    )

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Build filename suffix based on optional features
    suffix_parts = []
    if msa is not None:
        suffix_parts.append("_custom_msa")
    if additional_pairs:
        suffix_parts.append("_constrained")
    suffix = "".join(suffix_parts)
    
    filename = out_dir / f"{seq_name}_{ligand_name}{suffix}.yaml"
    with open(filename, "w") as fh:
        yaml.dump(
            data,
            fh,
            Dumper=BoltzDumper,
            default_flow_style=False,
            indent=4,
            sort_keys=False,
            allow_unicode=True,
        )

    return filename


# ── Batch generator ───────────────────────────────────────────────────────────


def generate_boltz2_yamls(
    entries: list[dict],
    output_dir: str = "yaml",
    stem_length: int = 8,
    max_distance: float = 5.0,
    msa: Optional[str] = None,
    additional_pairs: AdditionalPairs | None = None,
    molecule_type: str = "dna",
    **kwargs,
) -> list[Path]:
    """
    Generate one YAML file per (sequence, ligand) entry.

    Each dict in *entries* must contain::

        seq_name    : str   — sequence name (used in filenames)
        sequence    : str   — nucleotide / amino-acid sequence
        ligand_name : str   — e.g. "C0R", "HCY", "11DC", or "free"
        ccd         : str | None
        smiles      : str | None

    Optional keys that override the function-level defaults per entry::

        chain_id, ligand_chain, cyclic, stem_length, max_distance, msa

    molecule_type    : ``"dna"``, ``"rna"``, or ``"protein"``.
                       Controls the sequence-type key in the YAML output.
    stem_length      : number of Watson-Crick base-pair constraints at each stem
                       end; set to 0 to omit stem constraints entirely (only
                       applies to ``"dna"`` and ``"rna"``).
    max_distance     : upper-bound distance (Å) for each contact.
    msa              : optional path to an MSA file; when provided it is written
                       as an ``msa`` key inside the sequence block of every entry
                       that does not supply its own ``msa`` key.  Pass ``None``
                       (default) to omit the field entirely.
    additional_pairs : flexible specification of additional contact constraints.
                       See :func:`_resolve_pairs` for the full list of accepted
                       formats — flat list, per-sequence dict, or per-sequence +
                       per-ligand dict with ``"*"`` wildcard support.

    Returns a list of Paths of the generated files.
    """

    print(f"\nGenerating Boltz-2 YAML files in {output_dir} …\n")

    generated = []
    for entry in entries:
        # Resolve additional_pairs for this (seq_name, ligand_name)
        resolved_pairs = _resolve_pairs(
            additional_pairs,
            entry["seq_name"],
            entry["ligand_name"],
        )
        # Per-entry additional_pairs override the resolved value
        entry_pairs = entry.get("additional_pairs", resolved_pairs)

        # Always generate the unconstrained (base) version
        path_base = generate_yaml_file(
            seq_name=entry["seq_name"],
            sequence=entry["sequence"],
            ligand_name=entry["ligand_name"],
            ccd=entry.get("ccd") or None,
            smiles=entry.get("smiles") or None,
            output_dir=output_dir,
            stem_length=entry.get("stem_length", stem_length),
            max_distance=entry.get("max_distance", max_distance),
            msa=entry.get("msa", msa),
            additional_pairs=None,
            molecule_type=entry.get("molecule_type", molecule_type),
            **{
                k: entry[k]
                for k in ("chain_id", "ligand_chain", "cyclic")
                if k in entry
            },
            **kwargs,
        )
        print(f"  ✓ {path_base}")
        generated.append(path_base)

        # If there are additional pairs, also generate the constrained version
        if entry_pairs:
            path_constrained = generate_yaml_file(
                seq_name=entry["seq_name"],
                sequence=entry["sequence"],
                ligand_name=entry["ligand_name"],
                ccd=entry.get("ccd") or None,
                smiles=entry.get("smiles") or None,
                output_dir=output_dir,
                stem_length=entry.get("stem_length", stem_length),
                max_distance=entry.get("max_distance", max_distance),
                msa=entry.get("msa", msa),
                additional_pairs=entry_pairs,
                molecule_type=entry.get("molecule_type", molecule_type),
                **{
                    k: entry[k]
                    for k in ("chain_id", "ligand_chain", "cyclic")
                    if k in entry
                },
                **kwargs,
            )
            print(f"  ✓ {path_constrained}")
            generated.append(path_constrained)
    
    print(f"\n{len(generated)} file(s) written to {output_dir} \n")

    return generated


# ── Example usage ─────────────────────────────────────────────────────────────

if __name__ == "__main__":
    """
    Example: generate YAML files for a small aptamer project with one
    sequence and two ligands (one CCD-based, one SMILES-based).
    """

    SEQUENCES = {
        "Apt1": "GGGACGACGCCCGCATGTTCCATGGATAGTCTTGACTAGTCGTCCC",
        "Apt2": "GGGACGACTAGCGTATGCGCCAGAAGTATACGAGGATAGTCGTCCC",
    }

    LIGANDS = {
        "free": {"ccd": None, "smiles": None},
        "L01":  {"ccd": "C0R", "smiles": None},
        "L02":  {
            "ccd": None,
            "smiles": "C1=CC=C2C(=C1)C(=O)C3=C(C2=O)C4=CC=CC=C4C3=O",
        },
    }

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

    # ── Example: flat-list additional_pairs ────────────────────────────────
    print("Generating YAML files with flat-list additional_pairs …")
    paths = generate_boltz2_yamls(
        entries,
        output_dir="example_project/yaml",
        stem_length=8,
        additional_pairs=[(1, 44), (2, 43)],
    )

    # ── Example: per-sequence additional_pairs ────────────────────────────
    print("\nGenerating YAML files with per-sequence additional_pairs …")
    per_seq_pairs = {
        "Apt1": [("A", 5, "A", 36)],
        "Apt2": [("A", 3, "B", 10)],
    }
    paths = generate_boltz2_yamls(
        entries,
        output_dir="example_project/yaml_per_seq",
        stem_length=8,
        additional_pairs=per_seq_pairs,
    )

    # ── Example: per-sequence + per-ligand with wildcard ─────────────────
    print("\nGenerating YAML files with per-seq + per-ligand additional_pairs …")
    per_lig_pairs = {
        "Apt1": {
            "free": [("A", 1, "A", 44)],
            "L01":  [("A", 5, "B", 10)],
            "*":    [("A", 2, "A", 43)],
        },
    }
    paths = generate_boltz2_yamls(
        entries,
        output_dir="example_project/yaml_per_lig",
        stem_length=8,
        additional_pairs=per_lig_pairs,
    )

    print("\nAll examples completed.")
