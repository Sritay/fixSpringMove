/* ----------------------------------------------------------------------
  This fix works with the LAMMPS software package (post Sept 2021 version)
  This fix enables atoms to be constrained by a 1-D spring to a moveble 
  virtual point, which is initialized to the atom position at fix
  initiation.
  To compile: place fix_spring_move.cpp and fix_spring_move.h in the 
  LAMMPS src folder aand recompile.
  Usage: fix fixID groupID spring/self/move k xyz v_1 v_2 v_3
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Sritay Mistry (University of Edinburgh)
------------------------------------------------------------------------- */

#include <cstdlib>
#include <cstring>
#include "fix_spring_move.h"
#include "atom.h"
#include "modify.h"
#include "update.h"
#include "domain.h"
#include "respa.h"
#include "memory.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "variable.h"

enum{EQUAL,ATOM};
using namespace LAMMPS_NS;
using namespace FixConst;

/* ----------------------------------------------------------------------
   Constructor to read variables from input file and initialise positions
------------------------------------------------------------------------- */

FixSpringSelfMove::FixSpringSelfMove(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  xoriginal(NULL), xvarstr(NULL), yvarstr(NULL), zvarstr(NULL)
{
  if ((narg < 8) || (narg > 8))
    error->all(FLERR,"Illegal fix spring/self/move command");

  restart_peratom = 1;
  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;
  respa_level_support = 1;

  k = utils::numeric(FLERR,arg[3],false,lmp);
  if (k <= 0.0) error->all(FLERR,"Illegal fix spring/self/move command");

/* read equal style variables from input file to add to atom velocity */

    if (strcmp(arg[5],"NULL") == 0) xvarstr = NULL;
    else if (strstr(arg[5],"v_") == arg[5]) {
      int n = strlen(&arg[5][2]) + 1;
      xvarstr = new char[n];
      strcpy(xvarstr,&arg[5][2]);
    } else error->all(FLERR,"Illegal fix spring move command");
    if (strcmp(arg[6],"NULL") == 0) yvarstr = NULL;
    else if (strstr(arg[6],"v_") == arg[6]) {
      int n = strlen(&arg[6][2]) + 1;
      yvarstr = new char[n];
      strcpy(yvarstr,&arg[6][2]);
    } else error->all(FLERR,"Illegal fix spring move command");
    if (strcmp(arg[7],"NULL") == 0) zvarstr = NULL;
    else if (strstr(arg[7],"v_") == arg[7]) {
      int n = strlen(&arg[7][2]) + 1;
      zvarstr = new char[n];
      strcpy(zvarstr,&arg[7][2]);
    } else error->all(FLERR,"Illegal fix spring move command");

/* variable reading section ends*/

/* set direction of motion specified in input file */
 
  xflag = yflag = zflag = 1;

  if (narg == 5) {
    if (strcmp(arg[4],"xyz") == 0) {
      xflag = yflag = zflag = 1;
    } else if (strcmp(arg[4],"xy") == 0) {
      zflag = 0;
    } else if (strcmp(arg[4],"xz") == 0) {
      yflag = 0;
    } else if (strcmp(arg[4],"yz") == 0) {
      xflag = 0;
    } else if (strcmp(arg[4],"x") == 0) {
      yflag = zflag = 0;
    } else if (strcmp(arg[4],"y") == 0) {
      xflag = zflag = 0;
    } else if (strcmp(arg[4],"z") == 0) {
      xflag = yflag = 0;
    } else error->all(FLERR,"Illegal fix spring/self/move command");
  }

  // perform initial allocation of atom-based array
  // register with Atom class

  xoriginal = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);
  atom->add_callback(1);

  /* save initial unwrapped atom position in the variable xoriginal */

  double **x = atom->x;
  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) domain->unmap(x[i],image[i],xoriginal[i]);
    else xoriginal[i][0] = xoriginal[i][1] = xoriginal[i][2] = 0.0;
  }

  espring = 0.0;

}

/* ---------------------------------------------------------------------- */

FixSpringSelfMove::~FixSpringSelfMove()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,0);
  atom->delete_callback(id,1);

  // delete locally stored array

  memory->destroy(xoriginal);
}

/* ---------------------------------------------------------------------- */

int FixSpringSelfMove::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ----------------------------------------------------------------------
   Verfiy equal style variables read from input file
------------------------------------------------------------------------- */

void FixSpringSelfMove::init()
{
  if (strstr(update->integrate_style,"respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }

/* checking variable existence and ensuring only equal style variables 
   are provided as input arguments */

    if (xvarstr) {
      xvar = input->variable->find(xvarstr);
      if (xvar < 0) error->all(FLERR,
                               "Variable name for fix move does not exist");
      if (input->variable->equalstyle(xvar)) xvarstyle = EQUAL;
      else error->all(FLERR,"Variable for fix move is invalid style");
    }
    if (yvarstr) {
      yvar = input->variable->find(yvarstr);
      if (yvar < 0) error->all(FLERR,
                               "Variable name for fix move does not exist");
      if (input->variable->equalstyle(yvar)) yvarstyle = EQUAL;
      else error->all(FLERR,"Variable for fix move is invalid style");
    }
    if (zvarstr) {
      zvar = input->variable->find(zvarstr);
      if (zvar < 0) error->all(FLERR,
                               "Variable name for fix move does not exist");
      if (input->variable->equalstyle(zvar)) zvarstyle = EQUAL;
      else error->all(FLERR,"Variable for fix move is invalid style");
    }
}

/* ---------------------------------------------------------------------- */

void FixSpringSelfMove::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixSpringSelfMove::min_setup(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
  Evaluate equal style variable and move xoriginal and atom postions
  by the evaluated value
------------------------------------------------------------------------- */

void FixSpringSelfMove::initial_integrate(int /*vflag*/)
{

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double mov1x = 0.0;
  double mov1y = 0.0;
  double mov1z = 0.0;
  
    modify->clearstep_compute();

    if (xvarstr){ 
      mov1x = input->variable->compute_equal(xvar);}
    
    if (yvarstr)
      mov1y = input->variable->compute_equal(yvar);

    if (zvarstr)
      mov1z = input->variable->compute_equal(zvar);

    modify->addstep_compute(update->ntimestep + 1);

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {

      	xoriginal[i][0] += mov1x;
      	xoriginal[i][1] += mov1y;
      	xoriginal[i][2] += mov1z;
        x[i][0] += mov1x;
        x[i][1] += mov1y;
        x[i][2] += mov1z;
      }

}

/* ----------------------------------------------------------------------
   Calculate forces and energies due to this fix 
------------------------------------------------------------------------- */

void FixSpringSelfMove::post_force(int vflag)
{

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  double dx,dy,dz;
  double unwrap[3];

  espring = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - (xoriginal[i][0]);
      dy = unwrap[1] - (xoriginal[i][1]);
      dz = unwrap[2] - (xoriginal[i][2]);
      if (!xflag) dx = 0.0;
      if (!yflag) dy = 0.0;
      if (!zflag) dz = 0.0;
      f[i][0] -= k*dx;
      f[i][1] -= k*dy;
      f[i][2] -= k*dz;
      espring += k * (dx*dx + dy*dy + dz*dz);
    }

  espring *= 0.5;

}

/* ---------------------------------------------------------------------- */

void FixSpringSelfMove::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSpringSelfMove::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   Compute total energy from all atoms due to this fix
------------------------------------------------------------------------- */

double FixSpringSelfMove::compute_scalar()
{
  double all;
  MPI_Allreduce(&espring,&all,1,MPI_DOUBLE,MPI_SUM,world);
  return all;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixSpringSelfMove::memory_usage()
{
  double bytes = atom->nmax*3 * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixSpringSelfMove::grow_arrays(int nmax)
{
  memory->grow(xoriginal,nmax,3,"fix_spring/self/move:xoriginal");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixSpringSelfMove::copy_arrays(int i, int j, int delflag)
{
  xoriginal[j][0] = xoriginal[i][0];
  xoriginal[j][1] = xoriginal[i][1];
  xoriginal[j][2] = xoriginal[i][2];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixSpringSelfMove::pack_exchange(int i, double *buf)
{
  buf[0] = xoriginal[i][0];
  buf[1] = xoriginal[i][1];
  buf[2] = xoriginal[i][2];
  return 3;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixSpringSelfMove::unpack_exchange(int nlocal, double *buf)
{
  xoriginal[nlocal][0] = buf[0];
  xoriginal[nlocal][1] = buf[1];
  xoriginal[nlocal][2] = buf[2];
  return 3;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixSpringSelfMove::pack_restart(int i, double *buf)
{
  buf[0] = 4;
  buf[1] = xoriginal[i][0];
  buf[2] = xoriginal[i][1];
  buf[3] = xoriginal[i][2];
  return 4;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixSpringSelfMove::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  xoriginal[nlocal][0] = extra[nlocal][m++];
  xoriginal[nlocal][1] = extra[nlocal][m++];
  xoriginal[nlocal][2] = extra[nlocal][m++];
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixSpringSelfMove::maxsize_restart()
{
  return 4;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixSpringSelfMove::size_restart(int nlocal)
{
  return 4;
}
