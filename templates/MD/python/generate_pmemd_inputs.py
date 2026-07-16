#!/usr/bin/env python
"""
Amber pmemd input file generator
==================================

Generates all 11 Amber pmemd ``.in`` input files (steps 1-10 + final_min)
tailored to your system parameters.  Output is written to ``pmemd/in/``
in the current working directory.

Usage::

    python python/generate_pmemd_inputs.py \\
        --macromol-type dna \\
        --macromol-start 1 \\
        --macromol-end 46 \\
        --ligand-idx 47 \\
        --production-ns 100

Existing files in ``pmemd/in/`` are moved to ``pmemd/in/archive/`` before
generation.  See ``--help`` for all options.
"""

import argparse
import math
import pathlib
import shutil


BACKBONE_ATOMS = {
    "dna": "P,OP1,OP2,O5',C5',C4',C3',O3'",
    "rna": "P,OP1,OP2,O5',C5',C4',C3',O3',O2'",
    "protein": "N,CA,C,O",
}

MOL_TYPE_LABEL = {
    "dna": "DNA",
    "rna": "RNA",
    "protein": "Protein",
}


def fmt_nstlim(prod_ns: int) -> int:
    return int(prod_ns * 1000 / 0.002)


def fmt_int_with_commas(n: int) -> str:
    return f"{n:,}"


def build_step1(solute_end: int, mac_start: int, mac_end: int,
                ligand_idx: int | None, mol_label: str, lig_label: str) -> str:
    return f"""Minimize solvent, heavy atom restraints on solute
 &cntrl
  imin=1,                   ! Energy minimization flag (1 = minimization)  
  ntmin=2,				    ! Only steepest descent method is used
  
  ! Minimization steps
  maxcyc=1000,              ! Total number of minimization steps (1000 steps)
  ncyc=10,                  ! Number of steps for steepest descent if ntmin=1 (here, 1000 steps of steepest descent, no conjugate gradient)
  
  ! Nonbonded interactions
  cut=8.0,                 ! Nonbonded cutoff distance in A (8 A)

  ! Boundary conditions
  ntb=1,                    ! Constant volume (no pressure coupling)
  
  ! Constraints
  ntc=1, ntf=1,             ! no SHAKE 
  
  ! Output control
  ntpr=50,                  ! Print energy and other information every 50 steps
  ntwx=500,                 ! Write coordinates (trajectory) every 500 steps
  ntwr=500,                 ! Write restart files every 1000 steps

  ! Restraints
  ntr=1,                    ! Positional restraints are applied
  restraint_wt=5.0,         ! Restraint force constant in kcal/mol-A2 (5.0 kcal/mol-A2)
  restraintmask='(:{mac_start}-{solute_end} & !@H=)', ! Restraints applied to the heavy atoms of solute ({mol_label} residues {mac_start}-{mac_end}{lig_label})

  ! NMR restraints
  nmropt=0                  ! Disable NMR restraints 
 /
"""


def build_step2(solute_end: int, mac_start: int, mac_end: int,
                ligand_idx: int | None, mol_label: str, lig_label: str) -> str:
    return f"""MD simulation with positional restraints
 &cntrl
  imin=0,                   ! MD simulation (not minimization)
  irest=0,                  ! Starting fresh, no restart
  ntx=1,                    ! Coordinates are read from the standard coordinate file
  ig=-1,		     		! Random number generator will use a seed based on the current time

  ! Simulation time control
  nstlim=15000,             ! 15,000 steps for the simulation (15 ps total with a 1 fs timestep)
  dt=0.001,                 ! Timestep of 1 fs

  ! Temperature control
  temp0=300.0,              ! Target temperature (300 K)
  tempi=300.0,              ! Initial temperature (300 K); initial velocities to the atoms assigned according to Maxwell-Boltzmann distribution
  ntt=3, gamma_ln=5.0,      ! Langevin thermostat with friction coefficient of 5.0 (1/ps)

  ! Boundary conditions
  ntb=1,                    ! Constant volume (NVT ensemble)
  ntp=0,                    ! No pressure coupling (constant volume)
  iwrap=0,					! Coordinates will not be wrapped
  nscm=0,				    ! Center of mass motion will not be removed

  ! Constraints
  ntc=2, ntf=2,             ! SHAKE on bonds involving hydrogen atoms

  ! Nonbonded interactions
  cut=8.0,                 ! Nonbonded interactions cutoff at 8 A

  ! Output control
  ntpr=50,                 ! Print energy and other information every 50 steps
  ntwx=500,                 ! Write coordinates (trajectory) every 500 steps
  ntwr=500,                ! Write restart files every 1000 steps
  ntwv=-1,					! Velocities not written to file

  ! Restraints
  ntr=1,                    ! Apply positional restraints
  restraint_wt=5.0,         ! Restraint force constant of 5.0 kcal/mol-A2
  restraintmask=":{mac_start}-{solute_end}&!@H=", ! Apply restraints to the heavy atoms of solute ({mol_label} residues {mac_start}-{mac_end}{lig_label})

  ! NMR restraints
  nmropt=0                  ! Disable NMR restraints 
 /
"""


def build_step3(solute_end: int, mac_start: int, mac_end: int,
                ligand_idx: int | None, mol_label: str, lig_label: str) -> str:
    return f"""Minimize solute with medium restraints
 &cntrl
  imin=1,                   ! Energy minimization flag (1 = minimization)
  ntmin=2,			        ! Only steepest descent method is used

  ! Minimization steps
  maxcyc=1000,              ! Total number of minimization steps (1000 steps)
  ncyc=10,                  ! Number of steps for steepest descent if ntmin=1 (here, 1000 steps of steepest descent)

  ! Nonbonded interactions
  cut=8.0,                 ! Nonbonded cutoff distance in A (8 A)

  ! Boundary conditions
  ntb=1,                    ! Constant volume (no pressure coupling)
  
  ! Constraints
  ntc=1, ntf=1,             ! no SHAKE 
  
  ! Output control
  ntpr=50,                  ! Print energy and other information every 50 steps
  ntwx=500,                 ! Write coordinates (trajectory) every 500 steps
  ntwr=500,                 ! Write restart files every 1000 steps

  ! Restraints
  ntr=1,                    ! Positional restraints are applied
  restraint_wt=2.0,         ! Restraint force constant in kcal/mol-A2 (2.0 kcal/mol-A2)
  restraintmask=":{mac_start}-{solute_end}&!@H=", ! Restraints applied to the heavy atoms of solute ({mol_label} residues {mac_start}-{mac_end}{lig_label})

  ! NMR restraints
  nmropt=0                  ! Disable NMR restraints
 /
"""


def build_step4(solute_end: int, mac_start: int, mac_end: int,
                ligand_idx: int | None, mol_label: str, lig_label: str) -> str:
    return f"""Minimize solute with weak restraints
 &cntrl
  imin=1,                   ! Energy minimization flag (1 = minimization)
  ntmin=2,			        ! Only steepest descent method is used

  ! Minimization steps
  maxcyc=1000,              ! Total number of minimization steps (1000 steps)
  ncyc=10,                  ! Number of steps for steepest descent if ntmin=1 (here, 1000 steps of steepest descent)

  ! Nonbonded interactions
  cut=8.0,                 ! Nonbonded cutoff distance in A (8 A)

  ! Boundary conditions
  ntb=1,                    ! Constant volume (no pressure coupling)
  
  ! Constraints
  ntc=1, ntf=1,             ! no SHAKE 
  
  ! Output control
  ntpr=50,                  ! Print energy and other information every 50 steps
  ntwx=500,                 ! Write coordinates (trajectory) every 500 steps
  ntwr=500,                 ! Write restart files every 1000 steps

  ! Restraints
  ntr=1,                    ! Positional restraints are applied
  restraint_wt=0.1,         ! Restraint force constant in kcal/mol-A2 (0.1 kcal/mol-A2)
  restraintmask=":{mac_start}-{solute_end}&!@H=", ! Restraints applied to the heavy atoms of solute ({mol_label} residues {mac_start}-{mac_end}{lig_label})

  ! NMR restraints
  nmropt=0                  ! Disable NMR restraints
 /
"""


def build_step5() -> str:
    return """Final minimization without restraints
 &cntrl
  imin=1,                   ! Energy minimization flag (1 = minimization)
  ntmin=2,			        ! Only steepest descent method is used

  ! Minimization steps
  maxcyc=1000,              ! Total number of minimization steps (1000 steps)
  ncyc=10,                  ! Number of steps for steepest descent if ntmin=1 (here, 1000 steps of steepest descent)

  ! Nonbonded interactions
  cut=8.0,                 ! Nonbonded cutoff distance in A (8 A)

  ! Boundary conditions
  ntb=1,                    ! Constant volume (no pressure coupling)
  
  ! Constraints
  ntc=1, ntf=1,             ! no SHAKE 
  
  ! Output control
  ntpr=50,                  ! Print energy and other information every 50 steps
  ntwx=500,                 ! Write coordinates (trajectory) every 500 steps
  ntwr=500,                 ! Write restart files every 1000 steps

  ! Restraints
  ntr=0,                    ! No positional restraints applied

  ! NMR restraints
  nmropt=0                  ! Disable NMR restraints
 /
"""


def build_step6(solute_end: int, mac_start: int, mac_end: int,
                ligand_idx: int | None, mol_label: str, lig_label: str) -> str:
    return f"""Initial NPT with mild positional restraints
 &cntrl
  imin=0,                    ! MD simulation flag (0 = no minimization, i.e., run dynamics)
  irest=0,                   ! Initial restart flag (0 = start from scratch, no previous trajectory)
  ntx=1,                     ! Input coordinate type (1 = read in coordinates from the restart file)
  ig=-1,		     		 ! Random number generator will use a seed based on the current time

  ! Simulation time control
  nstlim=5000,               ! Number of timesteps for the simulation (5000 timesteps)
  dt=0.001,                  ! Timestep size in ps (0.001 ps)

  ! Temperature control
  temp0=300.0,               ! Target temperature in Kelvin (300 K)
  tempi=300.0,               ! Initial temperature in Kelvin (300 K)
  ntt=3,                     ! Temperature coupling (3 = Langevin dynamics)
  gamma_ln=5.0,              ! Langevin collision frequency (5 ps^-1)

  ! Pressure control
  ntb=2,                     ! Periodic boundary conditions with constant pressure (NPT ensemble)
  ntp=1,                     ! Pressure coupling (1 = isotropic pressure coupling)
  pres0=1.0,                 ! Reference pressure (1 atm)
  barostat=2,                ! Monte Carlo barostat (2 = Monte Carlo barostat)
  mcbarint=100,              ! Monte Carlo barostat interval (attempt volume change every 100 steps)
  taup=1.0,                  ! Pressure coupling time constant (1.0 ps)
  iwrap=0,					 ! Coordinates will not be wrapped
  nscm=0,				     ! Center of mass motion will not be removed

  ! Constraints
  ntc=2,                     ! Constrained bonds involving hydrogen atoms (2 = SHAKE algorithm applied)
  ntf=2,                     ! Force evaluation for constrained bonds (2 = forces on constrained bonds are calculated)

  ! Nonbonded interactions
  cut=8.0,                   ! Nonbonded cutoff distance in A (8 A)

  ! Output control
  ntpr=50,                   ! Print energy and other information every 50 steps
  ntwx=500,                  ! Write coordinates (trajectory) every 500 steps
  ntwr=500,                  ! Write restart files every 1000 steps
  ntwv=-1,					 ! Velocities not written to file

  ! Restraints
  ntr=1,                     ! Apply positional restraints (1 = apply restraints)
  restraint_wt=1.0,          ! Restraint weight for the positional restraints (mild restraint strength)
  restraintmask=":{mac_start}-{solute_end}&!@H=", ! Atom selection for the restraints ({mol_label} residues {mac_start}-{mac_end}{lig_label})

  ! NMR restraints
  nmropt=0                  ! Disable NMR restraints
 /
"""


def build_step7(solute_end: int, mac_start: int, mac_end: int,
                ligand_idx: int | None, mol_label: str, lig_label: str) -> str:
    return f"""Continue NPT, reduce positional restraints
 &cntrl
  imin=0,                    ! MD simulation flag (0 = no minimization, i.e., run dynamics)
  irest=1,                   ! Restart flag (1 = continue from a previous trajectory or restart file)
  ntx=5,                     ! Input coordinate type (5 = read coordinates and velocities from a restart file)
  ig=-1,		     		 ! Random number generator will use a seed based on the current time

  ! Simulation time control
  nstlim=5000,               ! Number of timesteps for the simulation (5000 timesteps)
  dt=0.001,                  ! Timestep size in ps (0.001 ps)

  ! Temperature control
  temp0=300.0,               ! Target temperature in Kelvin (300 K)
  tempi=300.0,               ! Initial temperature in Kelvin (300 K)
  ntt=3,                     ! Temperature coupling (3 = Langevin dynamics)
  gamma_ln=5.0,              ! Langevin collision frequency (5 ps^-1)

  ! Pressure control
  ntb=2,                     ! Periodic boundary conditions with constant pressure (NPT ensemble)
  ntp=1,                     ! Pressure coupling (1 = isotropic pressure coupling)
  pres0=1.0,                 ! Reference pressure (1 atm)
  barostat=2,                ! Monte Carlo barostat (2 = Monte Carlo barostat)
  mcbarint=100,              ! Monte Carlo barostat interval (attempt volume change every 100 steps)
  taup=1.0,                  ! Pressure coupling time constant (1.0 ps)
  iwrap=0,					 ! Coordinates will not be wrapped
  nscm=0,				     ! Center of mass motion will not be removed

  ! Constraints
  ntc=2,                     ! Constrained bonds involving hydrogen atoms (2 = SHAKE algorithm applied)
  ntf=2,                     ! Force evaluation for constrained bonds (2 = forces on constrained bonds are calculated)

  ! Nonbonded interactions
  cut=8.0,                   ! Nonbonded cutoff distance in A (8 A)

  ! Output control
  ntpr=50,                   ! Print energy and other information every 50 steps
  ntwx=500,                  ! Write coordinates (trajectory) every 500 steps
  ntwr=500,                  ! Write restart files every 1000 steps
  ntwv=-1,					 ! Velocities not written to file

  ! Restraints
  ntr=1,                     ! Apply positional restraints (1 = apply restraints)
  restraint_wt=0.5,          ! Restraint weight for the positional restraints (milder restraint strength)
  restraintmask=":{mac_start}-{solute_end}&!@H=", ! Atom selection for the restraints ({mol_label} residues {mac_start}-{mac_end}{lig_label})

  ! NMR restraints
  nmropt=0                  ! Disable NMR restraints
 /
"""


def build_step8(mac_start: int, mac_end: int,
                backbone_atoms: str,
                ligand_idx: int | None,
                mol_label: str, lig_label: str) -> str:
    bb_selection = f"(:{mac_start}-{mac_end}@{backbone_atoms})"
    if ligand_idx is not None:
        lig_selection = f"|(:{ligand_idx}&!@H=)"
    else:
        lig_selection = ""
    mask = f"{bb_selection}{lig_selection}"
    return f"""Backbone-only restraints
 &cntrl
  imin=0,                    ! MD simulation flag (0 = no minimization, i.e., run dynamics)
  irest=1,                   ! Restart flag (1 = continue from a previous trajectory or restart file)
  ntx=5,                     ! Input coordinate type (5 = read coordinates and velocities from a restart file)
  ig=-1,		     		 ! Random number generator will use a seed based on the current time

  ! Simulation time control
  nstlim=10000,              ! Number of timesteps for the simulation (10000 timesteps)
  dt=0.001,                  ! Timestep size in ps (0.001 ps)

  ! Temperature control
  temp0=300.0,               ! Target temperature in Kelvin (300 K)
  tempi=300.0,               ! Initial temperature in Kelvin (300 K)
  ntt=3,                     ! Temperature coupling (3 = Langevin dynamics)
  gamma_ln=5.0,              ! Langevin collision frequency (5 ps^-1)

  ! Pressure control
  ntb=2,                     ! Periodic boundary conditions with constant pressure (NPT ensemble)
  ntp=1,                     ! Pressure coupling (1 = isotropic pressure coupling)
  pres0=1.0,                 ! Reference pressure (1 atm)
  barostat=2,                ! Monte Carlo barostat (2 = Monte Carlo barostat)
  mcbarint=100,              ! Monte Carlo barostat interval (attempt volume change every 100 steps)
  taup=1.0,                  ! Pressure coupling time constant (1.0 ps)
  iwrap=0,					 ! Coordinates will not be wrapped
  nscm=0,				     ! Center of mass motion will not be removed

  ! Constraints
  ntc=2,                     ! Constrained bonds involving hydrogen atoms (2 = SHAKE algorithm applied)
  ntf=2,                     ! Force evaluation for constrained bonds (2 = forces on constrained bonds are calculated)

  ! Nonbonded interactions
  cut=8.0,                  ! Nonbonded cutoff distance in A (8 A)

  ! Output control
  ntpr=50,                   ! Print energy and other information every 50 steps
  ntwx=500,                  ! Write coordinates (trajectory) every 500 steps
  ntwr=500,                  ! Write restart files every 1000 steps
  ntwv=-1,					 ! Velocities not written to file

  ! Restraints
  ntr=1,                     ! Apply positional restraints (1 = apply restraints)
  restraint_wt=0.5,          ! Restraint weight for the positional restraints (milder restraint strength)
  restraintmask="{mask}", ! {mol_label} backbone atoms{lig_label}

  ! NMR restraints
  nmropt=0                  ! Disable NMR restraints
 /
"""


def build_step9() -> str:
    return """Short NPT without positional restraints
 &cntrl
  imin=0,                    ! MD simulation flag (0 = no minimization, i.e., run dynamics)
  irest=1,                   ! Restart flag (1 = continue from a previous trajectory or restart file)
  ntx=5,                     ! Input coordinate type (5 = read coordinates and velocities from a restart file)
  ig=-1,		     		 ! Random number generator will use a seed based on the current time

  ! Simulation time control
  nstlim=5000,              ! Number of timesteps for the simulation (5000 timesteps)
  dt=0.002,                  ! Timestep size in ps (0.002 ps)

  ! Temperature control
  temp0=300.0,               ! Target temperature in Kelvin (300 K)
  tempi=300.0,               ! Initial temperature in Kelvin (300 K)
  ntt=3,                     ! Temperature coupling (3 = Langevin dynamics)
  gamma_ln=5.0,              ! Langevin collision frequency (5 ps^-1)

  ! Pressure control
  ntb=2,                     ! Periodic boundary conditions with constant pressure (NPT ensemble)
  ntp=1,                     ! Pressure coupling (1 = isotropic pressure coupling)
  pres0=1.0,                 ! Reference pressure (1 atm)
  barostat=2,                ! Monte Carlo barostat (2 = Monte Carlo barostat)
  mcbarint=100,              ! Monte Carlo barostat interval (attempt volume change every 100 steps)
  taup=1.0,                  ! Pressure coupling time constant (1.0 ps)
  iwrap=0,					 ! Coordinates will not be wrapped
  nscm=1000,				 ! Removal of translational and rotational center-of-mass (COM) motion every 1000 steps (default)

  ! Constraints
  ntc=2,                     ! Constrained bonds involving hydrogen atoms (2 = SHAKE algorithm applied)
  ntf=2,                     ! Force evaluation for constrained bonds (2 = forces on constrained bonds are calculated)

  ! Nonbonded interactions
  cut=8.0,                  ! Nonbonded cutoff distance in A (8 A)

  ! Output control
  ntpr=50,                   ! Print energy and other information every 50 steps
  ntwx=500,                  ! Write coordinates (trajectory) every 500 steps
  ntwr=500,                  ! Write restart files every 1000 steps
  ntwv=-1,					 ! Velocities not written to file

  ! Restraints
  ntr=0,                     ! No positional restraints applied

  ! NMR restraints
  nmropt=0                  ! Disable NMR restraints
 /
"""


def build_step10(prod_nstlim: int, prod_ns: int, nstlim_str: str) -> str:
    return f"""Long NPT production, {prod_ns} ns
 &cntrl
  imin=0,                    ! MD simulation flag (0 = no minimization, i.e., run dynamics)
  irest=1,                   ! Restart flag (1 = continue from a previous trajectory or restart file)
  ntx=5,                     ! Input coordinate type (5 = read coordinates and velocities from a restart file)
  ig=-1,		     		 ! Random number generator will use a seed based on the current time
  
  ! Simulation time control
  nstlim={prod_nstlim},            ! Number of timesteps for the simulation ({nstlim_str} timesteps = {prod_ns} ns with dt=0.002 ps)
  dt=0.002,                  ! Timestep size in ps (0.002 ps)

  ! Temperature control
  temp0=300.0,               ! Target temperature in Kelvin (300 K)
  tempi=300.0,               ! Initial temperature in Kelvin (300 K)
  ntt=3,                     ! Temperature coupling (3 = Langevin dynamics)
  gamma_ln=5.0,              ! Langevin collision frequency (5 ps^-1)

  ! Pressure control
  ntb=2,                     ! Periodic boundary conditions with constant pressure (NPT ensemble)
  ntp=1,                     ! Pressure coupling (1 = isotropic pressure coupling)
  pres0=1.0,                 ! Reference pressure (1 atm)
  barostat=2,                ! Monte Carlo barostat (2 = Monte Carlo barostat)
  mcbarint=100,              ! Monte Carlo barostat interval (attempt volume change every 100 steps)
  iwrap=0,					 ! Coordinates will not be wrapped
  nscm=1000,				 ! Removal of translational and rotational center-of-mass (COM) motion every 1000 steps (default)

  ! Constraints
  ntc=2,                     ! Constrained bonds involving hydrogen atoms (2 = SHAKE algorithm applied)
  ntf=2,                     ! Force evaluation for constrained bonds (2 = forces on constrained bonds are calculated)

  ! Nonbonded interactions
  cut=9.0,                   ! Nonbonded cutoff distance in A (9 A)

  ! Output control
  ntpr=500,                  ! Print energy information every 5000 steps
  ntwx=5000,                 ! Write coordinates to the trajectory file every 5000 steps
  ntwr=50000,                ! Write restart files every 50000 steps
  ntwv=-1,					 ! Velocities not written to file

  ! Restraints
  ntr=0,                     ! No positional restraints (0 = no restraints)

  ! NMR restraints
  nmropt=0                  ! Disable NMR restraints
 /
"""


def build_final_min() -> str:
    return """Minimize solute with medium restraints
 &cntrl
  imin=1,                   ! Energy minimization flag (1 = minimization)
  ntx=1,					! Only coordinates read (no velocities; default value)
  irest=0,					! Run a new simulation (no restart, velocities ignored; default value)
  ntmin=1,			        ! Steepest descent method is used until ncyc is reached, then conjugate gradient (default)

  ! Minimization steps
  maxcyc=20000,             ! Total number of minimization steps (20,000 steps)
  ncyc=4000,                ! Method switched from steepest descent to conjugate gradient after 4000 cycles (if ntmin=1)

  ! Nonbonded interactions
  cut=8.0,                 ! Nonbonded cutoff distance in A (8 A)

  ! Boundary conditions
  ntb=1,                    ! Constant volume (no pressure coupling; Constant pressure is not used in minimization)
  
  ! Constraints
  ntc=1, ntf=1,             ! no SHAKE (complete interaction is calculated; default value)
  
  ! Output control
  ntpr=50,                  ! Print energy and other information every 50 steps
  ntwx=500,                 ! Write coordinates (trajectory) every 500 steps
  ntwr=500,                 ! Write restart files every 1000 steps

  ! Restraints
  ntr=0,                    ! No positional restraints are applied (default value)

  ! NMR restraints
  nmropt=0                  ! Disable NMR restraints
 /
"""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate Amber pmemd input files for MD preparation + production."
    )
    parser.add_argument(
        "--macromol-type",
        choices=["dna", "rna", "protein"],
        default="dna",
        help="Macromolecule type (determines backbone atom set) [default: dna]",
    )
    parser.add_argument(
        "--macromol-start",
        type=int,
        default=1,
        help="First residue of the macromolecule [default: 1]",
    )
    parser.add_argument(
        "--macromol-end",
        type=int,
        default=46,
        help="Last residue of the macromolecule [default: 46]",
    )
    parser.add_argument(
        "--ligand-idx",
        type=int,
        default=None,
        help="Ligand residue index (omit if no ligand) [default: None]",
    )
    parser.add_argument(
        "--production-ns",
        type=int,
        default=100,
        help="Production simulation length in ns [default: 100]",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    out_dir = pathlib.Path.cwd() / "pmemd" / "in"
    archive_dir = out_dir / "archive"

    if out_dir.exists():
        archive_dir.mkdir(parents=True, exist_ok=True)
        for f in out_dir.glob("*.in"):
            shutil.move(str(f), str(archive_dir / f.name))
            print(f"Archived: {f.name} -> archive/")
    else:
        out_dir.mkdir(parents=True)

    solute_end = args.macromol_end
    ligand_idx = args.ligand_idx
    if ligand_idx is not None:
        solute_end = max(solute_end, ligand_idx)

    mol_label = MOL_TYPE_LABEL[args.macromol_type]
    backbone_atoms = BACKBONE_ATOMS[args.macromol_type]

    lig_label = ""
    if ligand_idx is not None:
        lig_label = f" + residue {ligand_idx}"

    prod_nstlim = fmt_nstlim(args.production_ns)
    nstlim_str = fmt_int_with_commas(prod_nstlim)

    files = {
        "step1.in": build_step1(
            solute_end, args.macromol_start, args.macromol_end,
            ligand_idx, mol_label, lig_label),
        "step2.in": build_step2(
            solute_end, args.macromol_start, args.macromol_end,
            ligand_idx, mol_label, lig_label),
        "step3.in": build_step3(
            solute_end, args.macromol_start, args.macromol_end,
            ligand_idx, mol_label, lig_label),
        "step4.in": build_step4(
            solute_end, args.macromol_start, args.macromol_end,
            ligand_idx, mol_label, lig_label),
        "step5.in": build_step5(),
        "step6.in": build_step6(
            solute_end, args.macromol_start, args.macromol_end,
            ligand_idx, mol_label, lig_label),
        "step7.in": build_step7(
            solute_end, args.macromol_start, args.macromol_end,
            ligand_idx, mol_label, lig_label),
        "step8.in": build_step8(
            args.macromol_start, args.macromol_end,
            backbone_atoms, ligand_idx,
            mol_label, lig_label),
        "step9.in": build_step9(),
        "step10.in": build_step10(prod_nstlim, args.production_ns, nstlim_str),
        "final_min.in": build_final_min(),
    }

    for name, content in files.items():
        path = out_dir / name
        path.write_text(content, encoding="utf-8", newline="\n")
        print(f"Generated: {path}")

    print(f"\nDone — {len(files)} files written to {out_dir}/")


if __name__ == "__main__":
    main()
