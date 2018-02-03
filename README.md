# Ring-Polymer Molecular Dynamics package for LAMMPS
This repository hosts the source code of a [LAMMPS](http://lammps.sandia.gov/) code extension in the form of a [fix command](http://lammps.sandia.gov/doc/fix.html) that allows Ring-Polymer Molecular Dynamics (RPMD) calculations to be performed. Details about the algorithm implementation and capabilities can be found in:

["Quantum effects on dislocation motion from Ring-Polymer Molecular Dynamics"  
Rodrigo Freitas, Mark Asta, and Vasily Bulatov  
arXiv preprint:1712.0462 (submitted)](https://arxiv.org/abs/1712.04629)

## Installation
Clone this repository to a subdirectory named USER-RPMD inside the `src/` directory of your LAMMPS installation:
```
cd <lammps-path>/src/
git clone https://github.com/freitas-rodrigo/QuantumEffectsDislocationMotion.git USER-RPMD
```
While still inside the `src/` of LAMMPS, activate the code extension package with the command
```
make yes-user-RPMD
```
Finally, [recompile](http://lammps.sandia.gov/doc/Section_start.html#making-lammps) LAMMPS by running `make <machine>` (for example: `make mpi`).

This package requires LAMMPS to be compiled using MPI libraries, the code will not run if you use a serial LAMMPS executable (the inter-partition communication is done using MPI library routines). But notice that you can use the parallel LAMMPS executable to run this code in serial (or in more partitions than the number of physical processors you have) since the MPI routines will simulate virtual processors.

## Usage
The only command that needs to be added to a LAMMPS input script in order to run a RPMD simulation is the
`fix rpmd`, the command's syntax is described below. Notice that `fix rpmd` should be declared _before_ any other fix that alters the forces on the particles (such as thermostats or constraints). This is necessary because one of the actions of `fix rpmd` is to scale down all forces acting on the particles, thus any forces added by fixes declared before `fix rpmd` will be altered.

The number of beads per ring polymer is determined by the use of the command-line flag `-partition` when running LAMMPS. For example, in order to run a RPMD simulation with 30 beads per ring polymer where 10 processors are used to run each of the 30 replicas of the system (resulting in a total of 300 processors) specified in the `in.lmp` LAMMPS input script we should use:
```
mpirun -n 300 lmp_mpi -in in.lmp -partition 30x10
```

### Fix syntax
```
fix ID group-ID rpmd T
```
- `ID`, `group-ID` are documented in [fix](http://lammps.sandia.gov/doc/fix.html) command. `group-ID` must be all.  
- `rpmd` = style name of this fix command  
- `T` = simulation temperature  

#### Example: 
```
fix f1 all rpmd 300
```
This fix computes a global scalar and a global vector quantities that can be accessed by various output commands. The scalar is the sum of the ring-polymers' spring energy for each partition, where the spring energy is 0.5 * k * r^2. The vector is the partition's ring-polymers contribution to the stress tensor. The tensor has 6 components and it is stored in the following order: xx, yy, zz, xy, xz, yz. The scalar and vector values calculated by this fix are extensive.

## Example: quantum harmonic oscillator.
Inside the example directory of this repository you should find scripts that compute the temperature dependence of the energy of a quantum harmonic oscillator. This is an instructional example because it shows how to use and post-process the data from fix rpmd.

`in.lmp`: LAMMPS script to simulate 125 harmonic oscillators using RPMD.  
`job.sh`: bash script to run the LAMMPS-RPMD simulations.

By Running `./job.sh` the simulations will run for 7 different temperatures with the number of beads ranging from 29 to 4. In my 2012 personal laptop it takes ~2min to run this set of simulations. Inside post_processing there are python scripts to gather the data, compute the total energy, and plot the result. You will need python's modules numpy and matplotlib installed to run these scripts.

`compute.py`: Reduce the data from all replicas and compute thte total energy.  
`plot.py`: plot the temperature dependence of the total energy.

If the scripts run sucessfully you would be able to recover the plot below. An extended analysis of these simulations is given in the Supplementary Information section of the [paper](https://arxiv.org/abs/1712.04629).

## Computing the quantum energy and stress

Computing the quantum energy and stress from the output of LAMMPS and fix rpmd can be nontrivial due how scattered the data can be between many files. Here is an algorithm and some code excerpts to help you out:

### Energy:

Let us assume you declared your fix rpmd as
```
fix f1 all rpmd 300
```
To extract the potential energy of your interatomic potential and the ring-polymer springs energy you should collect the data from these two variables:
```
variable pair_pe equal pe/atom
variable ring_pe equal f_f2/atoms
```
where ``pair_pe`` is a variable for the interatomic potential energy and ``ring_pe`` is for the ring-polymer springs energy. To compute the total energy for the quantum system according to equation 3 from the Supplementary Information you would use in python, assuming you saved the interatomic potential energy and the ring-polymer potential energy as columns in the file "potential_energy_P.dat" where P is the replica number:
```python
for ibead in range(1,nbeads+1):
    data = loadtxt('potential_energy_%d.dat' % ibead)
    pair_pe += data[:,0]
    ring_pe += data[:,1]
E = 1.5*kB*T*nbeads - ring_pe + harm_pe/nbeads
```
where kB is the Boltzmann constant and T is the system temperature declared in the fix.

### Stress tensor:

Let us assume the same setup as for the energy calculation. For the stress tensor we will be using equation 4 of the Supplementary Information. Suppose we are interested in the pressure only.
xx component of the stress tensor. In the LAMMPS script we would use
```
compute  c1 all pressure NULL pair # Compute virial part of the pressure.
variable pair_P equal c_c1
variable ring_P equal (f_f2[1]+f_f2[2]+f_f2[3])/(3*vol)
```  
```python
for ibead in range(1,nbeads+1):
    data = loadtxt('pressure_%d.dat' % ibead)
    pair_P += data[:,0]
    ring_P += data[:,1]
P = nbeads*kB*T/v * eVA3tobar - ring_P + pair_P/nbeads
```
Where v is the volume per atom and eVA3tobar is a conversion factor from eV/A^2 to bars (this fix has only been tested in metal units).
    
# Author & Contact

Rodrigo Freitas | rodrigof@berkeley.edu

