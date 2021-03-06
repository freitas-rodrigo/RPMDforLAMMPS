#!/bin/bash
# This bash script executes sequentially a series of LAMMPS-RPMD simulations at different temperatures.
# Usage: bash job.sh

lammps="../../lmp_mpi" # Path to LAMMPS executable.

# Create directory structure for data output.
mkdir -p data/logs
mkdir -p data/thermo
mkdir -p data/processed

# Setup list of parameter to loop over.
T=(50 150 250 350 450 550 650) # Temperature list.
BEADS=(29 10 6 5 4 4 4)        # Number of beads for each temperature.

# Loop over different temperatures.
for n in $(seq 0 6)
do
    printf "Running T = ${T[n]}K simulation with ${BEADS[n]} beads...\n"
    mpirun -n ${BEADS[n]} ${lammps}                  \
                   -in in.lmp                        \
                   -partition ${BEADS[n]}x1          \
                   -log       data/logs/${T[n]}K.log \
                   -screen    none                   \
                   -var       RANDOM ${RANDOM}       \
                   -var       nbeads ${BEADS[n]}     \
                   -var       T ${T[n]}
done
