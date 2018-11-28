#include "fix_rpmd.h"   // Fix rpmd class header.
#include "atom.h"       // Atomic position and forces.
#include "universe.h"   // Number of beads.
#include "comm.h"       // Number of processors per partition.
#include "memory.h"     // Memory allocation routines.
#include "error.h"      // Error messages.
#include "string.h"     // memcpy function.
#include "force.h"      // Spring constant calculation.
#include "math_const.h" // Math constants (MY_2PI).
#include "stdlib.h"     // atof function.
#include "domain.h"     // pbc function.
#include "citeme.h"     // Create citation file.

using namespace LAMMPS_NS; // LAMMPS namespace.
using namespace FixConst;  // Fix constants namespace (mask ID).
using namespace MathConst; // Math constants namespace.

static const char cite_fix_rpmd[] =
  "fix rpmd command:\n\n"
  "@article{freitas2018,\n"
  "  author={Freitas, Rodrigo and Asta, Mark and Bulatov, Vasily V},\n"
  "  title={Quantum effects on dislocation motion from Ring-Polymer Molecular Dynamics},\n"
  "  journal={npj Computational Materials},\n"
  "  year={2018},\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

// Class constructor: arg = { fix_ID, group_ID, fix_name, arg_1, ... }
FixRPMD::FixRPMD(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (lmp->citeme) lmp->citeme->add(cite_fix_rpmd);

  // Verify number of arguments.
  if (narg != 4)
    error->universe_one(FLERR,"Illegal fix rpmd command");

  // Setup fix flags.
  vector_flag = 1;
  size_vector = 6;
  scalar_flag = 1;
  global_freq = 1;
  extvector = 1;
  extscalar = 1;

  double T = atof(arg[3]); // Read fix keywords (temperature [K]).

  // Compute spring constant.
  P = universe->nworlds;                  // Number of beads.
  double hbar = force->hplanck / MY_2PI;  // Planck constant [eV.ps].
  double beta = 1.0 / (force->boltz * T); // Thermodynamic beta [1/eV].
  double omega = sqrt(P) / (beta * hbar); // Spring frequency [1/ps^2].
  double m = atom->mass[atom->type[0]];   // Atomic mass [g/mol].
  k = m * omega * omega * force->mvv2e;   // Spring constant [eV/A^2].

  energy = 0.0; // Spring energy.

  x_neigh = new double * [2];
  x_neigh[0] = NULL;
  x_neigh[1] = NULL;

  comm_init(); // Initialize communication variables.
}

/* ---------------------------------------------------------------------- */

// Class destructor.
FixRPMD::~FixRPMD()
{
  memory->destroy(x_neigh); // Ring polymer variables.

  // Communications variables.
  int c_nprocs = comm->nprocs;
  memory->destroy(lc_tracker);
  memory->destroy(x_send);
  memory->destroy(x_recv);
  for (int ibead = 0; ibead < 2; ibead++)
  {
      memory->destroy(max_nbt[ibead]);
      memory->destroy(nb_nlocal[ibead]);
      memory->destroy(nb_ranks[ibead]);
      memory->destroy(nb_tracker[ibead]);
  }
}

/* ---------------------------------------------------------------------- */

// Determine at which point in the Verlet algorithm this fix is called.
int FixRPMD::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;   // Force rescaling + spring contribution.
  mask |= PRE_NEIGHBOR; // Rebuild tracker after atoms migrate.
  return mask;
}

/* ---------------------------------------------------------------------- */

// Called at the beginning of each run.
void FixRPMD::init()
{
  build_tracker(); // Build atom tracker before run starts.
}

/* ---------------------------------------------------------------------- */

// Called after atoms migrate but before neighbor list is built.
void FixRPMD::pre_neighbor()
{
  build_tracker(); // Rebuild atom tracker after atoms migrate.
}

/* ---------------------------------------------------------------------- */

// Called after pairwise forces are computed.
void FixRPMD::post_force(int vflag)
{
  int nlocal = atom->nlocal; // Number of atom in this processor.
  double **x = atom->x;      // Coordinates.
  double **f = atom->f;      // Forces.
  double dx, dy, dz;         // Distance between connected atoms.

  // Scale pair-wise forces by 1/P.
  for(int i = 0; i < nlocal; i++)
  {
    f[i][0] /= P;
    f[i][1] /= P;
    f[i][2] /= P;
  }

  comm_exchange(x); // Update the coordinate of atoms on neighbor beads.

  // Clean variables to accumulate spring stress and energy.
  energy = 0.0;
  for (int i = 0; i < 6; i++)
    stress[i] = 0.0;

  // Add spring force to atoms (loop over neighbor beads).
  for (int ibead = 0; ibead < 2; ibead++)
  {
    for (int i = 0; i < nlocal; i++)
    {
      // Compute distance between atoms connected by springs (apply pbc).
      dx = x[i][0] - x_neigh[ibead][3*i+0];
      dy = x[i][1] - x_neigh[ibead][3*i+1];
      dz = x[i][2] - x_neigh[ibead][3*i+2];
      domain->minimum_image(dx, dy, dz);

      // Add spring force and accumulate spring energy and stress.
      f[i][0] -= k*dx;
      f[i][1] -= k*dy;
      f[i][2] -= k*dz;
      energy += dx*dx + dy*dy + dz*dz;
      if (ibead == 0) {
        stress[0] += dx*dx;
        stress[1] += dy*dy;
        stress[2] += dz*dz;
        stress[3] += dx*dy;
        stress[4] += dx*dz;
        stress[5] += dy*dz;
      }
    }
  }
  energy *= 0.5;   // Eliminate double counting.
  energy *= k/2.0; // k/2 factor from spring energy.
  // k factor from spring and conversion to stress units.
  for (int i = 0; i < 6; i++)
    stress[i] *= k*force->nktv2p;
}

/* ---------------------------------------------------------------------- */

// Collect stress accumulator from each processor.
double FixRPMD::compute_vector(int n)
{
  double all;
  MPI_Allreduce(&stress[n],&all,1,MPI_DOUBLE,MPI_SUM,world);
  return all;
}

/* ---------------------------------------------------------------------- */

double FixRPMD::compute_scalar()
{
  double all;
  MPI_Allreduce(&energy,&all,1,MPI_DOUBLE,MPI_SUM,world);
  return all;
}

/* ---------------------------------------------------------------------- */
/* Interpartition communication methods.                                  */
/* ---------------------------------------------------------------------- */

// Initialize communication variables and flags.
void FixRPMD::comm_init()
{
  int me = universe->me;
  int c_nprocs = comm->nprocs;
  int u_nprocs = universe->nprocs;
  int iworld = universe->iworld;

  // Initializes buffer for neighbor beads' atoms.
  max_nlocal = 0;
  x_recv = NULL;
  x_send = NULL;

  // Initialize trackers.
  lc_tracker = NULL;
  nb_tracker[0] = new int * [c_nprocs];
  nb_tracker[1] = new int * [c_nprocs];
  for (int ibead = 0; ibead < 2; ibead++)
  {
    for (int nrank = 0; nrank < c_nprocs; nrank++)
      nb_tracker[ibead][nrank] = NULL;
  }

  // Initialize maximum buffer size.
  max_lct = 0;
  max_nbt[0] = new int [c_nprocs];
  max_nbt[1] = new int [c_nprocs];
  for (int ibead = 0; ibead < 2; ibead++)
  {
    for (int nrank = 0; nrank < c_nprocs; nrank++)
      max_nbt[ibead][nrank] = 0;
  }

  // Initialize vector with number of atoms on neighbor bead's processors.
  nb_nlocal[0] = new int [c_nprocs];
  nb_nlocal[1] = new int [c_nprocs];

  // Define rank of neighbor beads' processors.
  // This processor will communicates with the processor of its neighbor beads
  // in a specific order. First it talks to its equivalent processor in the next
  // bead. Then it talks to the next one, and so on. "Periodic boundary
  // conditions" need to be applied to make sure we talk to the correct 
  // beads (partitions).
  nb_ranks[0] = new int [c_nprocs];
  nb_ranks[1] = new int [c_nprocs];

  // Set the rank of equivalent processor in the next bead (partition).
  nb_ranks[0][0] = me - c_nprocs;
  nb_ranks[1][0] = me + c_nprocs;
  if(nb_ranks[0][0] < 0) nb_ranks[0][0] += u_nprocs; // "PBC".
  if(nb_ranks[1][0] >= u_nprocs) nb_ranks[1][0] -= u_nprocs; // "PBC".

  // Define limit of the processors of next bead (partition).
  int upper_rank[2];
  upper_rank[0] = iworld * c_nprocs;
  upper_rank[1] = (iworld+2) * c_nprocs;
  if (upper_rank[0] == 0) upper_rank[0] += u_nprocs; // "PBC".
  if (upper_rank[1] > u_nprocs) upper_rank[1] -= u_nprocs; // "PBC".

  // Set up rank of neighbor beads.
  for (int nrank = 0; nrank < c_nprocs; nrank++)
  {
    nb_ranks[0][nrank] = nb_ranks[0][0] + c_nprocs - nrank;
    nb_ranks[1][nrank] = nb_ranks[1][0] + nrank;
    // Impose "PBC".
    if (nb_ranks[0][nrank] >= upper_rank[0]) 
      nb_ranks[0][nrank] -= c_nprocs;
    if (nb_ranks[1][nrank] >= upper_rank[1]) 
      nb_ranks[1][nrank] -= c_nprocs;
  }
}

/* ---------------------------------------------------------------------- */

// Construct local atom tracker and communicates it to neighbors.
void FixRPMD::build_tracker()
{
  int natoms = atom->natoms;
  int nlocal = atom->nlocal;
  int c_nprocs = comm->nprocs;

  // Grow local tracker size if necessary.
  if (nlocal > max_lct)
  {
    max_lct = nlocal + 100;
    memory->grow(lc_tracker, max_lct,"fix_rpmd:lc_tracker");
  }

  // Build local tracker.
  int t_size = 0;
  for (int gID = 1; gID < natoms+1; gID++)
  {
    if((atom->map(gID) >= 0) && (atom->map(gID) < nlocal))
      lc_tracker[t_size++] = gID;
  }

  // Exchange tracker with neighbor beads.
  // First the nlocal from the neighbor beads are obtained. Then the local 
  // buffer for the tracker is grown if necessary. Finally the trackers are
  // exchanged.

  // Loop over neighbor beads.
  for (int ibead = 0; ibead < 2; ibead++)
  {
    int prev_bead = ibead%2;
    int next_bead = (ibead+1)%2;
    // Loop over processors of neighbor beads.
    for (int nrank = 0; nrank < c_nprocs; nrank++)
    {
      // Send and receive local number of atoms.
      MPI_Sendrecv(
        &nlocal,                      1, MPI_INT, nb_ranks[prev_bead][nrank], 0,
        &nb_nlocal[next_bead][nrank], 1, MPI_INT, nb_ranks[next_bead][nrank], 0,
                   universe->uworld, MPI_STATUS_IGNORE);

      // Grow buffer to receive tracker.
      if (nb_nlocal[next_bead][nrank] > max_nbt[next_bead][nrank])
      {
        max_nbt[next_bead][nrank] = nb_nlocal[next_bead][nrank] + 100;
        memory->grow(nb_tracker[next_bead][nrank], max_nbt[next_bead][nrank]*sizeof(int), "FixRPMD:nb_tracker");
      }

      // Send and receive atoms tracker.
      MPI_Sendrecv(lc_tracker, nlocal, MPI_INT, nb_ranks[prev_bead][nrank], 0,
                      nb_tracker[next_bead][nrank], nb_nlocal[next_bead][nrank],
                                       MPI_INT, nb_ranks[next_bead][nrank],0,
                      universe->uworld, MPI_STATUS_IGNORE);
    }
  }
}

/* ---------------------------------------------------------------------- */

// Exchange atomic coordinates with neighbor beads.
void FixRPMD::comm_exchange(double **x)
{
  int nlocal = atom->nlocal;
  int c_nprocs = comm->nprocs;

  // Grow buffer to store coordinates (if necessary).
  if (nlocal > max_nlocal)
  {
    max_nlocal = nlocal + 100;
    int size = 3*max_nlocal * sizeof(double);
    memory->grow(x_send, size, "FixRPMD:x_send");
    memory->grow(x_recv, size, "FixRPMD:x_recv");
    memory->grow(x_neigh[0], size, "FixRPMD:x_neigh[0]");
    memory->grow(x_neigh[1], size, "FixRPMD:x_neigh[1]");
  }

  // Exchange atoms with neighbor beads.

  // Loop over neighbor beads.
  for (int ibead = 0; ibead < 2; ibead++)
  {
    int prev_bead = ibead%2;
    int next_bead = (ibead+1)%2;
    // Loop over processors of neighbor beads.
    for (int nrank = 0; nrank < c_nprocs; nrank++)
    {

      // Wrap local coordinates.
      double *x_ptr = x_send;
      int sshared = 0;
      for (int i = 0; i < nb_nlocal[next_bead][nrank]; i++) 
      {
        // Finds local index for each atom in neigh_bead.
        int index = atom->map(nb_tracker[next_bead][nrank][i]);

        // Copy coordinates to x_send in sequence. (x_ptr > x_send)
        if ((index >= 0) && (index < nlocal)) {
          sshared++;
          memcpy(x_ptr, x[index], 3*sizeof(double));
          x_ptr += 3;
        }
      }

      // Send and receive coordinates.
      int rshared;
      MPI_Sendrecv(&sshared, 1, MPI_INT, nb_ranks[next_bead][nrank], 0,
                   &rshared, 1, MPI_INT, nb_ranks[prev_bead][nrank], 0,
                   universe->uworld, MPI_STATUS_IGNORE);
      MPI_Sendrecv(x_send,3*sshared,MPI_DOUBLE, nb_ranks[next_bead][nrank], 0,
                   x_recv,3*rshared,MPI_DOUBLE, nb_ranks[prev_bead][nrank], 0,
                   universe->uworld, MPI_STATUS_IGNORE);

      // Unpack received coordinates.
      x_ptr = x_recv;
      for (int i = 0; i < nb_nlocal[prev_bead][nrank]; i++)
      {
        // Find local index for each atom in neigh_bead.
        int index = atom->map(nb_tracker[prev_bead][nrank][i]);

        // Copy coordinates from x_recv to x_neigh (x_ptr > x_recv)
        if ((index >= 0) && (index < nlocal)) {
          x_neigh[ibead][3*index+0] = x_ptr[0];
          x_neigh[ibead][3*index+1] = x_ptr[1];
          x_neigh[ibead][3*index+2] = x_ptr[2];
          x_ptr += 3;
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */
