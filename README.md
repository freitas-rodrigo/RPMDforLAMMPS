# Ring-Polymer Molecular Dynamics package for LAMMPS
This repository hosts the source code of a [LAMMPS](http://lammps.sandia.gov/) code extension in the form of a [fix command](http://lammps.sandia.gov/doc/fix.html) that allows Ring-Polymer Molecular Dynamics (RPMD) simulations to be performed. Details about the algorithm implementation and capabilities can be found in:

["Quantum effects on dislocation motion from Ring-Polymer Molecular Dynamics"  
Rodrigo Freitas, Mark Asta, and Vasily Bulatov  
nature partner journal: Computational Materials (Open Access)  
DOI:10.1038/s41524-018-0112-9](https://doi.org/10.1038/s41524-018-0112-9)

## Installation
Clone this repository to a subdirectory named `USER-RPMD` inside the `src/` directory of your LAMMPS installation:
```
cd <lammps-path>/src/
git clone https://github.com/freitas-rodrigo/RPMDforLAMMPS.git USER-RPMD
```
While still inside the `src/`, activate the code extension package with the command:
```
make yes-user-RPMD
```
Finally, [recompile](http://lammps.sandia.gov/doc/Section_start.html#making-lammps) LAMMPS by running `make <machine>` (for example: `make mpi`).

This package requires LAMMPS to be compiled using MPI libraries, the code will not run if you use a serial LAMMPS executable because the RPMD inter-partition communication is done using MPI library routines. But notice that you can use the parallel LAMMPS executable to run this code in serial (or in more partitions than the number of physical processors you have) because MPI routines are able to simulate virtual processors.

## Usage
The command that needs to be added to a LAMMPS input script in order to run a RPMD simulation is the
`fix rpmd`, the command's syntax is described below. It is also necessary to add the two following commands to your LAMMPS script in order for `fix rpmd` to work properly:
```
atom_modify map array
neigh_modify delay 0 every XYZ check no
```
The first command instructs LAMMPS to create a map between the local index of atoms and their global ID's (as explained [here](http://lammps.sandia.gov/doc/atom_modify.html)), `fix rpmd` uses this map during the simulation to perform inter-partition communication. The second command instructs LAMMPS to rebuild the neighbor list every XYZ steps, as explained [here](http://lammps.sandia.gov/doc/neigh_modify.html). This command is necessary in order to synchronize the neighbor list creation for the different partitions, `fix rpmd` requires this synchronization in order to rebuild and communicate its own list of neighbor beads. It is the user's responsibility to pick a reasonable value for XYZ based on the physical properties of the system being simulated. A good strategy is to run a classical MD simulation of the same system (using the same number of processors as the number of processors per partition to be used on the RPMD simulation) and then to choose XYZ (conservatively) based on the average time between neighbor list rebuilds.

Notice that `fix rpmd` should be declared _before_ any other fix that modifies the forces on the particles (such as thermostats or constraints). This is necessary because one of the actions of `fix rpmd` is to scale all forces acting on the particles, hence any forces added by `fix` commands declared before `fix rpmd` will be modified.

The number of beads per ring polymer is determined by the use of the command-line flag `-partition` when running LAMMPS. For example, in order to run a RPMD simulation with 30 beads per ring polymer where 10 processors are used to run each of the 30 replicas of the system (resulting in a total of 300 processors) specified in the `in.lmp` LAMMPS input script, we should use:
```
mpirun -n 300 lmp_mpi -in in.lmp -partition 30x10
```

#### Fix syntax
```
fix ID group-ID rpmd T
```
- `ID` and `group-ID` are documented in [fix](http://lammps.sandia.gov/doc/fix.html) command. `group-ID` must be all.  
- `rpmd` = style name of this fix command  
- `T` = simulation temperature  

##### Example: 
```
fix f1 all rpmd 300
```
This fix computes a global scalar and a global vector quantities that can be accessed by various output commands. The scalar is the sum of the ring-polymers' spring energy for each partition, where the spring energy is 0.5 * k * r^2. The vector is the partition's ring-polymers contribution to the stress tensor. The tensor has 6 components and it is stored in the following order: xx, yy, zz, xy, xz, yz. The scalar and vector values calculated by this fix are extensive.

## Example: quantum harmonic oscillator.
Inside the [`example`](example/) directory of this repository you should find scripts for simulations to compute the temperature dependence of the energy of a quantum harmonic oscillator. The scripts included also illustrate how to post-process data from `fix rpmd`. You will need python's modules [numpy](http://www.numpy.org/) and [matplotlib](https://matplotlib.org/) installed to run the post-processing scripts.

#### Scripts included
[`in.lmp`](example/in.lmp): LAMMPS script to simulate 125 harmonic oscillators using RPMD.  
[`job.sh`](example/job.sh): bash script to run the LAMMPS-RPMD simulations.  
[`post_processing/compute.py`](example/post_processing/compute.py): reduce the data from all replicas and compute the total energy.  
[`post_processing/plot.py`](example/post_processing/plot.py): plot the temperature dependence of the total energy against analytical results.  

#### Running the example scripts
From the [`example`](example/) directory use the command `bash job.sh` to run the RPMD-LAMMPS simulations for a total of 7 different temperatures where the number of beads ranges from 29 to 4, it should take 1-3 minutes to run the simulations. Notice that you might need to edit the path to your LAMMPS executable in [`job.sh`](example/job.sh/). From inside [`example/post_processing/`](example/post_processing/) use `python compute.py` to gather the data generated by `fix rpmd` and compute the total energy per harmonic oscillator. Finally, run `python plot.py` to plot the results. If the scripts ran successfully you will obtain a `pdf` with the following plot:

<p align="center">
  <img src="https://image.ibb.co/mnxXdx/fig_energy_vs_temperature.png" width="400"/>
</p>

An extended analysis of these simulations, including convergence with respect to number of beads, is given in the Supplementary Information section of the [paper](https://doi.org/10.1038/s41524-018-0112-9).

## Computing the quantum energy and stress
Computing the quantum energy and stress from the output of LAMMPS and `fix rpmd` can be nontrivial due how scattered the data can become. Find below an algorithm and some code excerpts that can help with these tasks. 

#### Energy:
Let us assume `fix rpmd` was declared in the LAMMPS script as
```
fix f1 all rpmd 300
```
The potential energy of the interatomic potential and ring-polymers' springs can be collected into the following two LAMMPS variables:
```
variable pair_pe equal pe/atom
variable ring_pe equal f_f2/atoms
```
Where ``pair_pe`` is for the interatomic potential energy and ``ring_pe`` is for the ring-polymers springs energy. The total energy of the quantum systems is computed according to equation (2) from the [Supplementary Information](https://doi.org/10.1038/s41524-018-0112-9). Assuming the data from `pair_pe` and `ring_pe` were stored as two columns in files named `potential_energy_N.dat`, where `N` is the partition number, the following python script can be used to compute the total energy of the quantum system:
```python
for ibead in range(1,nbeads+1):
    data = loadtxt('potential_energy_%d.dat' % ibead)
    pair_pe += data[:,0]
    ring_pe += data[:,1]
E = 1.5*kB*T*nbeads - ring_pe + pair_pe/nbeads
```
where `kB` is the Boltzmann constant, `T` is the system temperature declared in `fix rpmd`, `nbeads` is the total number of beads per ring polymer, and `E` is the total energy of the quantum system.

#### Stress tensor:
Let us assume the same setup presented for the energy calculation above. The stress tensor of the quantum system is obtained through equation (3) from [Supplementary Information](https://doi.org/10.1038/s41524-018-0112-9). For simplicity we will be computing the system's pressure, given by equation (4) from [Supplementary Information](https://doi.org/10.1038/s41524-018-0112-9). In the LAMMPS script we use
```
compute c1 all pressure NULL pair
variable pair_P equal c_c1
variable ring_P equal (f_f2[1]+f_f2[2]+f_f2[3])/(3*vol)
```  
Then, assuming the data from these variables were stored as two columns in files named `pressure_N.dat` (where `N` is the partition number), the following python script can be used to compute quantum system pressure:
```python
for ibead in range(1,nbeads+1):
    data = loadtxt('pressure_%d.dat' % ibead)
    pair_P += data[:,0]
    ring_P += data[:,1]
P = nbeads*kB*T/v * eVA3tobar - ring_P + pair_P/nbeads
```
Where `v` is the volume per atom and `eVA3tobar` is the units conversion factor from electron-Volt/Angstrom^2 to bars.
    
## Author & Contact

Rodrigo Freitas | freitas@stanford.edu
