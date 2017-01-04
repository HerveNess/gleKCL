/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(gle2baths,FixGLE2B)

#else

#ifndef LMP_FIX_GLE2B_H
#define LMP_FIX_GLE2B_H

#include "fix.h"

namespace LAMMPS_NS {

class FixGLE2B : public Fix {
 public:
  FixGLE2B(class LAMMPS *, int, char **);
  virtual ~FixGLE2B();
  int setmask();
  virtual void init();
  virtual void initial_integrate(int);
  virtual void final_integrate();
  void g12_system_bath(double, double **, double **);
  void calc_pair_force(int, int, double **);
  void calc_force_on_bath(double **);
  void calc_force_on_bath_gen(double **, double **);
  void calc_polarisation_force(int, double *, double **, double **);
  void calc_polarisation_force_diag(int, double *, double **);
  virtual void initial_integrate_respa(int, int, int);
  virtual void final_integrate_respa(int, int);
  virtual void reset_dt();

/*  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  int pack_restart(int, double *);
  void unpack_restart(int, int);
  int size_restart(int);
  int maxsize_restart();
  void init_s_gld();
*/  
 protected:
  FILE *fp1, *fp2;
  char *infile; // infile: file to read in virtual DOF tau & omega
  char *infile2; // infile2: file to read in bath atom information (incl. c coeffs, mass, mapping to atom id     
  // bathgroup, sysgroup: integer ids of system and bath groups
  // bathgroupbit, sysgroupbit: bitwise interpretation of groups
  int bathgroupL, sysgroup, bathgroupbitL, sysgroupbit;
  int bathgroupR, bathgroupbitR;
  // numinbath, numinbath2: count of atoms in thermal bath group
  int numinbathL,numinbath2; 
  int numinbathR; 
  // *gle_map: map of read in bath atom index to global atom id
  int *gle_mapL;
  int *gle_mapR;
  // nvdof: number of virtual DOF in system
  int nvdofL;
  int nvdofR;
  int me;
  // gle_temp_bath: temperature of thermal bath 
  // gle_delta: increment of change in position of atoms used to determine gamma12
  // gle_mubar: 'mass' of virtual degrees of freedom
  double gle_temp_bathL, gle_delta, gle_mubarL;
  double gle_temp_bathR, gle_mubarR;
  // gle_rseed: random number generator seed
  int gle_rseed;
  // *gle_tau: tau values of virtual DOF 
  // *gle_omega: omega values of virtual DOF
  double *gle_tauL, *gle_omegaL;
  double *gle_tauR, *gle_omegaR;
  // *gle_mbath: mass of bath atoms
  double *gle_mbathL;
  double *gle_mbathR;
  // ***gle_c: c coefficients of bath atoms
  double ***gle_cL;
  double ***gle_cR;
  // **gle_sk: sk values for virtual DOF
  double **gle_skL;
  double **gle_skR;
  // ****gle_g12: gamma12 values 
  double ****gle_g1L;
  double ****gle_g1R;
  // *gle_ak, *gle_bk:
  double *gle_akL, *gle_bkL;
  double *gle_akR, *gle_bkR;
  double **gle_farray1,**gle_fp1L,**gle_fm1L,**gle_flagrangeL;
  double **gle_farray2,**gle_fp1R,**gle_fm1R,**gle_flagrangeR;
  double dtv,dtf,dtv2;
  double *step_respa;
  double **cutsq;
  class RanMars *random;
  class Pair *pair;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix gle requires a valid group id for bath atoms

Self-explanatory.

E: Fix gle requires a valid group id for system atoms

Self-explanatory.

E: Fix gle dof file must exist, have the correct number of lines, > 0 lines and be formatted correctly.

Self-explanatory.

E: Fix gle bath file must exist, have the correct number of lines, > 0 lines and be formatted correctly.

Self-explanatory.

E: Fix gle needs to use a pair potential that allows for a single pair to be calculated.

Self-explanatory.

*/
