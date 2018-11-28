"""
  This script plots the average energy computed from RPMD as a function of the temperature. It also computes and plots the analytical results for the energy of a quantum and a classical harmonic oscillator at finite temperature.

  Usage:
    python plot.py
"""
from numpy import *
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rc('xtick', labelsize=14)
mpl.rc('ytick', labelsize=14)
mpl.rc('axes', labelsize=18)
mpl.rc('legend', fontsize=18)

################################################################################
# Load and prepare data.                                                       #
################################################################################

# Input parameters.
workdir='../data'
m = 2.5 # Atomic mass [g/mol].
k = 1.0 # Harmonic oscillator spring constant [eV/A^2].

# Constants & units.
kB    = 8.6173324e-05     # Boltzman constant [eV/K].
hbar  = 6.58211928e-04    # Plank constant [eV.ps].
JtoeV = 6.2415093433e+18  # Conversion of J to eV.
mol   = 6.02214129e+23    # One mol.
mvv2e = 10*JtoeV/mol      # Multiply mv^2 by this to obtain eV.
omega = sqrt(k/(mvv2e*m)) # Angular frequency [1/ps].
ho_e  = hbar*omega        # Harmonic oscillator energy units [eV].

# Analytical results.
temp = arange(1,670)*kB/ho_e    # Temperature.
n = 1.0/(2*temp)                # Quantumness measure.
E_quantum = 1.5*cosh(n)/sinh(n) # Exact quantum mechanics result.
E_classic = 3.0*temp            # Exact classical mechanics result.

# RPMD results.
T, E = loadtxt(workdir+'/processed/energy_vs_T.dat', unpack=True)
T *= kB/ho_e
E /= ho_e

################################################################################
# Plot.                                                                        #
################################################################################

# Start Matplotlib.
fig = plt.figure(figsize=(8,6))
ax  = fig.add_axes([0.15, 0.15, 0.80, 0.80])

# Plot.
ax.plot(temp, E_quantum, ls='-', c='k', lw=3, label=r'$\left<E\right>_\mathrm{quantum}=\frac{3}{2} \hbar \omega \, \mathrm{coth} \left( \frac{\hbar\omega}{2k_\mathrm{B} T} \right)$')
ax.plot(temp, E_classic, ls='-', c='royalblue', lw=3, label=r'$\left<E\right>_\mathrm{classical} = 3k_\mathrm{B} T$') 
ax.plot(T, E, 'o', c='red', ms=10, mew=0, label="RPMD") 

# Details.
ax.set_xlabel(r'$k_\mathrm{B} T$ $[\hbar \omega]$')
ax.set_ylabel('Total Energy [$\hbar\omega$]')
ax.set_ylim(0,4.5)
ax.set_xlim(0, 1.4)
ax.legend(loc='best', frameon=False)

# Save.
fig.savefig("fig_energy_vs_temperature.png", dpi=300)
plt.close()

################################################################################
