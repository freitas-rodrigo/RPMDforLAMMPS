"""
  This script collects all data from the different beads (ie, partitions) and computes the average energy of the harmonic oscillator for each temperature.

  Usage:
    python compute.py
"""
from numpy import *

################################################################################
# Input parameters.                                                            #
################################################################################

workdir='../data'               # Work directory.
T_list = arange(50,750,100)     # Temperature [K].
nbeads_list = [29,10,6,5,4,4,4] # Total number of beads.
kB = 8.6173324e-05              # Boltzman constant [eV/K].
t0 = 20                         # Equilibration time [time stride].
m = 2.5                         # Atomic mass [g/mol].
k = 1.0                         # Spring constant [eV/A^2].

################################################################################
# Reduce data from rings and compute energy.                                   #
################################################################################

# Open output file.
f = open(workdir+'/processed/energy_vs_T.dat', 'w')
f.write('# Temperature [harmonic oscillator units] | Energy [harmonic oscillator units]\n')

# Loop over each temperature (and respective total number of beads).
for (T,nbeads) in zip(T_list, nbeads_list):
    # Loop over beads and sums up partitions energy.
    for ibead in range(1,nbeads+1): 
        data = loadtxt(workdir+'/thermo/thermo_%dK_%d.dat' % (T, ibead))
        if ibead == 1:
            harm_pe = data[:,1] # Harmonic oscillators potential energy.
            ring_pe = data[:,2] # Ring-polymer springs energy.
        else:
            harm_pe += data[:,1]
            ring_pe += data[:,2]
    # Compute total energy and its average.
    E = 1.5*kB*T*nbeads - ring_pe + harm_pe/nbeads
    E = E[t0:].mean()
    f.write('%3d %.3f\n' % (T, E))
f.close()

################################################################################
