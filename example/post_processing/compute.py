"""
  This script collects all data from the different beads (or replicas) and computes the average energy of the harmonic oscillator for each temperature.
"""
from numpy import *

################################################################################
# Input variables.                                                             #
################################################################################

# Parameters.
workdir='../data' # Work directory.
T_list = arange(50,750,100)                # Temperature [K].
nbeads_list = array([29,10,6,5,4,4,4])     # Total number of beads.
L = 100                                   # Number of output events.
t0 = 20                                   # Equilibration time [time stride].
m = 2.5                                    # Atomic mass [g/mol].
k = 1.0                                    # Spring constant [eV/A^2].
kB = 8.6173324e-05                         # Boltzman constant [eV/K].

################################################################################
# Reduce data from rings and compute energy.                                   #
################################################################################

# Open output file.
f = open(workdir+'/thermo/energy_vs_T.dat', 'w')
f.write('# T E [harmonic oscillator units]\n')

# Loop over each temperature (and respective total number of beads).
for (T,nbeads) in zip(T_list, nbeads_list):

    # Loop over beads and add total energy of each bead.
    harm_pe = zeros(L) # Harmonic potential energy.
    ring_pe = zeros(L) # Ring-polymer springs energy.
    for ibead in range(1,nbeads+1): 
        data = loadtxt(workdir+'/thermo/thermo_%dK_%d.dat' % (T, ibead))
        harm_pe += data[:L,1]
        ring_pe += data[:L,2]

    # Compute total energy and convert to harmonic oscillator units.
    E = 1.5*kB*T*nbeads - ring_pe + harm_pe/nbeads
    E = E[t0:].mean()
    f.write('%.2f %.3f\n' % (T, E))

f.close()

################################################################################
