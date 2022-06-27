/* -*- c++ -*- ----------------------------------------------------------
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

FixStyle(spring/self/move,FixSpringSelfMove)

#else

#ifndef LMP_FIX_SPRING_SELF_MOVE_H
#define LMP_FIX_SPRING_SELF_MOVE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSpringSelfMove : public Fix {
 public:
  FixSpringSelfMove(class LAMMPS *, int, char **);
  ~FixSpringSelfMove();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void initial_integrate(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double compute_scalar();

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  int pack_restart(int, double *);
  void unpack_restart(int, int);
  int size_restart(int);
  int maxsize_restart();

 private:
  char *xvarstr,*yvarstr,*zvarstr;
  double k,espring;
  double **xoriginal;         // original coords of atoms
  int xvar,yvar,zvar;
  int xvarstyle,yvarstyle,zvarstyle;
  int xflag, yflag, zflag;
  int ilevel_respa;
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
