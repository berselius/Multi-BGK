// Collection of 2D periodic poisson solvers
// For now, this assumes that all data is on a single node (e.g. no MPI)

#ifndef POISSON2DPERIODIC_H
#define POISSON2DPERIODIC_H

// Sets some parameters that the package will use
void init_poisson2D(unsigned Nx, unsigned Ny, double dx, double dy, double Lx,
                    double Ly);

// Solves classical Poisson with a linearized RHS
// Effectively this is a Helmholtz equation
void PoissonLinearperiodic2d(double **alpha, double **source, double **phi) {

  // Solves classical Poisson with a nonlinear RHS
  // Placeholder for now
  void PoissonNonlinearPeriodic2D(double **source, double **Te, double **phi);

  // Solves classical Poisson with a nonlinear RHS
  // Placeholder for now
  void PoissonLinearPeriodic2D_TF(double **source, double **Te, double **phi);

  // Solves classical Poisson with a nonlinear RHS
  // Placeholder for now
  void PoissonNonlinearPeriodic2D_TF(double **source, double **Te,
                                     double **phi);

#endif
