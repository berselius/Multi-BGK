// Finite Volume code for LHS of 1D-3V Boltzmann

// Data structure for species information
//#include "species.h"

// Sets up everything used in transport side inputs are:
//  numV - Number of velocity points in each coordinate direction
//  numX - Number of physical points in x-direction. Note - if doing this with
//  MPI, this is the number of points that are handled by the current process.
//  Periodic + MPI not currently implemented
void initialize_transport(int numV, int numX, int nspec, double *xnodes,
                          double *dxnodes, double Lx, double **vel, int ord,
                          double timestep, int bc);

/*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
// these are what you actually call from the outside
// f is the input distribution - first coordinate is the x location, second is
// species, third is the v location f_conv is the output
//
// f array sizes: f[numX][nspec][numV*numV*numV]

// id is the process id, so it knows if it's at a boundary or interior of domain
// for MPI communication

void advectOne(double ***f, double *PoisPot, double **qm, double m, int sp);

void advectTwo(double ***f, double *PoisPot, double **qm, double m, int sp);

void advectTwo_x(double ***f, double ***f_conv, int sp);

void advectTwo_v(double ***f, double ***f_conv, double *PoisPot, double **qm,
                 double m, int sp);

// void upwindTwo(double ***f, double ***f_conv, double *PoisPot, double **qm,
// double m, int sp);

/*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/

// deallocates stuff allocated in the initializer
void dealloc_trans();
