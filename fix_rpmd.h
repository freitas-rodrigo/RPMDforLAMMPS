#ifdef FIX_CLASS

FixStyle(rpmd,FixRPMD)

#else

#ifndef LMP_FIX_RPMD_H
#define LMP_FIX_RPMD_H

#include "fix.h"

// LAMMPS namespace.
namespace LAMMPS_NS {

class FixRPMD : public Fix {
 public:
  FixRPMD(class LAMMPS *, int, char **);
  ~FixRPMD();
  int setmask();
  virtual void init();
  virtual void pre_neighbor();
  virtual void post_force(int);
  double compute_scalar();
  double compute_vector(int);

 protected:
  // Ring polymer variables.
  int P;            // Number of beads.
  double k;         // Spring constant.
  double energy;    // Spring energy accumulator.
  double stress[6]; // Spring stress accumulator.
  double **x_neigh; // Coordinates of atoms on neighbor beads.

  // Communication methods.
  void comm_init();              // Initializes comm. variables and flags.
  void build_tracker();          // Builds tracker of atoms in this partition.
  void comm_exchange(double **); // Communicates coordinates to neighbor beads.

  // Communication variables.
  int max_nlocal;      // Size of buffer to store neighbor beads coordinates
  double *x_send;      // Buffer to send coordinates to neighbor beads.
  double *x_recv;      // Buffer to receive neighbor beads coordinates.
  int *lc_tracker;     // Tracker of global ID for local atoms.
  int **nb_tracker[2]; // Tracker of global ID for neighbor beads' atoms.
  int max_lct;         // Size of buffer to store local tracker.
  int *max_nbt[2];     // Size of buffer to store neighbor beads' tracker.
  int *nb_nlocal[2];   // Number of atoms on neighbor beads' processors.
  int *nb_ranks[2];    // Rank of neighbor beads processors.
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
