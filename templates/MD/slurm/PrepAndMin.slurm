#!/bin/bash
# =============================================================================
# PrepAndMin.slurm — Amber MD Preparation + Final Minimisation Pipeline
# =============================================================================
#
# SLURM array job that runs an 11-step Amber MD protocol:
#   Steps  1-10:  preparation (minimisation, NVT heating, NPT equilibration,
#                 NPT production)
#   Step  11:     final energy minimisation
#
# !! ADJUST FOR YOUR SYSTEM:
#   - JOB NAME  (#SBATCH -J)
#   - PARTITION (#SBATCH -p)
#   - ARRAY     (#SBATCH --array)
#   - EXPERIMENT name
#   - SCRATCHDIR path
#   - STEPS table (add/remove steps)
#   - cpptraj align and strip masks
#
# Reference: CSS1_HCY_constrained (46-nt DNA aptamer + HCY ligand)
# =============================================================================

#############################
# SLURM OPTIONS
#############################
#SBATCH -J __JOB_NAME__              !! ADJUST: job name
#SBATCH -t 23:59:59
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --gres gpu:1
#SBATCH -p __GPU_PARTITION__          !! ADJUST: cluster partition(s)
#SBATCH --array=__ARRAY__            !! ADJUST: array task IDs matching
                                     !!   your leap output files (cplx_N)
#############################

echo "#############################"
echo "User:" $USER
echo "Date:" $(date)
echo "Host:" $(hostname)
echo "PWD:" $(pwd)
echo "SLURM_JOBID:" $SLURM_JOBID
echo "SLURM_ARRAY_JOB_ID:" $SLURM_ARRAY_JOB_ID
echo "SLURM_ARRAY_TASK_ID:" $SLURM_ARRAY_TASK_ID
echo "#############################"

#############################
# IDENTIFIERS
#############################

TASK_ID=${SLURM_ARRAY_TASK_ID}
SUFFIX=_${TASK_ID}_J${SLURM_ARRAY_JOB_ID}

!! ADJUST: experiment name (must match leap/ directory name)
EXPERIMENT="__EXPERIMENT__"

###############################
# DIRECTORIES (SAFE ISOLATION)
###############################

!! ADJUST: path to scratch directory containing pmemd/ and leap/
SCRATCHDIR=__SCRATCHDIR__

RUNROOT="${SCRATCHDIR}/pmemd/out/${EXPERIMENT}/J${SLURM_ARRAY_JOB_ID}/task_${TASK_ID}"
PREPOUTDIR="${RUNROOT}/prep"
MINOUTDIR="${RUNROOT}/final_min"

mkdir -p "${PREPOUTDIR}/pdb"
mkdir -p "${MINOUTDIR}/pdb"

##############################
# TEMP WORKDIR (COMPUTE NODE)
##############################

TMPBASE=${TMPDIR:-/tmp}
WORKDIR="${TMPBASE}/md${SUFFIX}"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

echo "Running in TMPDIR: $WORKDIR"
echo "Output root: $RUNROOT"

#############################
# COPY INPUT FILES
#############################

echo "Staging input files"

cp ${SCRATCHDIR}/pmemd/in/step*.in .
cp ${SCRATCHDIR}/pmemd/in/final_min.in .

!! ADJUST: paths to leap-generated topology and coordinates
cp ${SCRATCHDIR}/leap/${EXPERIMENT}/cplx_${TASK_ID}.prmtop start.prmtop
cp ${SCRATCHDIR}/leap/${EXPERIMENT}/cplx_${TASK_ID}.rst7 start.rst7

#############################
# METADATA
#############################

echo "${SLURM_JOBID} ${TASK_ID} $(hostname)" > run.info

#############################
# MD + FINAL MIN TABLE
#############################

# Steps 1 to 10 are the preparation
# Step 10 has a 100 ns prod length (adjusted in step10.in)
# Step 11 is the final minimisation
#
# !! ADJUST: modify this table to add/remove/reorder steps
#
# Format: step  engine       input     ref     traj  infile         prefix
#   input:  "start" = start.rst7, otherwise ${input}${SUFFIX}.ncrst
#   ref:    "start" = start.rst7, "none" = no -ref, otherwise ${ref}${SUFFIX}.ncrst
#   traj:   "yes" = write .nc trajectory, "no" = skip

STEPS=(
"1  pmemd        start   start   no  step1.in      step1"
"2  pmemd.cuda   step1   start   yes step2.in      step2"
"3  pmemd        step2   start   no  step3.in      step3"
"4  pmemd        step3   start   no  step4.in      step4"
"5  pmemd        step4   none    no  step5.in      step5"
"6  pmemd.cuda   step5   step5   yes step6.in      step6"
"7  pmemd.cuda   step6   step5   yes step7.in      step7"
"8  pmemd.cuda   step7   step5   yes step8.in      step8"
"9  pmemd.cuda   step8   none    yes step9.in      step9"
"10 pmemd.cuda   step9   none    yes step10.in     step10"
"11 pmemd.cuda   step10  none    no  final_min.in  final_min"
)

echo "Starting MD pipeline"

#############################
# LOOP
#############################

for S in "${STEPS[@]}"; do
    read -r STEP ENGINE INPUT REF TRAJ INFILE PREFIX_BASE <<< "$S"

    PREFIX="${PREFIX_BASE}${SUFFIX}"

    # restart protection — skip if .ncrst already exists
    if [ -f "${PREFIX}.ncrst" ]; then
        echo "Skipping step ${STEP} (${PREFIX}) exists"
        continue
    fi

    # coordinate selection
    if [ "$INPUT" = "start" ]; then
        COORD_IN="start.rst7"
    else
        COORD_IN="${INPUT}${SUFFIX}.ncrst"
    fi

    CMD="$ENGINE -O \
        -i $INFILE \
        -p start.prmtop \
        -c $COORD_IN \
        -o ${PREFIX}.out \
        -r ${PREFIX}.ncrst \
        -inf ${PREFIX}.info"

    # restraints reference
    if [ "$REF" != "none" ]; then
        if [ "$REF" = "start" ]; then
            CMD="$CMD -ref start.rst7"
        else
            CMD="$CMD -ref ${REF}${SUFFIX}.ncrst"
        fi
    fi

    # trajectory
    if [ "$TRAJ" = "yes" ]; then
        CMD="$CMD -x ${PREFIX}.nc"
    fi

    echo "Running step ${STEP}: ${PREFIX_BASE}"
    eval $CMD
done

echo "Pipeline finished"

#############################
# COPY BACK PREP RESULTS
#############################

cp step*.out step*.ncrst step*.nc step*.info -t "${PREPOUTDIR}/"

for i in {1..10}; do
    ambpdb -p start.prmtop -c step${i}${SUFFIX}.ncrst > "${PREPOUTDIR}/pdb/step${i}${SUFFIX}.pdb"
done

#############################
# FINAL MIN OUTPUTS
#############################

cp final_min*.out final_min*.ncrst final_min*.info -t "${MINOUTDIR}/"

#############################
# CPPTRAJ
#############################

echo "Running cpptraj"

## ADJUST: alignment and strip masks for your system
##   align :N_DNA   — DNA/RNA residues for alignment  [currently :1-46]
##   strip :WAT     — solvent residue name             [usually :WAT]
##   :1-46 is for 46-nt DNA — change to match your solute

cat <<EOF > final_min_stripped${SUFFIX}.cpptraj
parm start.prmtop
trajin final_min${SUFFIX}.ncrst
autoimage
!! ADJUST: solute residues for alignment
align :__RESIDUES__
!! ADJUST: solvent residue name to strip
strip :__SOLVENT__ nobox parmout final_min${SUFFIX}_stripped.prmtop
trajout final_min${SUFFIX}_stripped.ncrst netcdf nobox
go
EOF

cat <<EOF > final_min${SUFFIX}.cpptraj
parm start.prmtop
trajin final_min${SUFFIX}.ncrst
autoimage
!! ADJUST: solute residues for alignment
align :__RESIDUES__
trajout final_min${SUFFIX}_cpptrajed.ncrst netcdf nobox
go
EOF

cpptraj -i final_min_stripped${SUFFIX}.cpptraj
cpptraj -i final_min${SUFFIX}.cpptraj

#############################
# FINAL PDB OUTPUT
#############################

ambpdb -p start.prmtop \
       -c final_min${SUFFIX}_cpptrajed.ncrst \
       > "${MINOUTDIR}/pdb/final_min${SUFFIX}_cpptrajed.pdb"

ambpdb -p final_min${SUFFIX}_stripped.prmtop \
       -c final_min${SUFFIX}_stripped.ncrst \
       > "${MINOUTDIR}/pdb/final_min${SUFFIX}_stripped.pdb"

#############################
# METADATA COPY
#############################

cp run.info "${RUNROOT}/"

echo "Congrats! All jobs completed successfully! You're amazing."
