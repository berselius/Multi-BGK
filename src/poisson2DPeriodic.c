// First stab at a 2D poisson solve. This assumes a single node (no mpi) and
// periodic BCs This uses a Gauss-Seidel solve, leveraging gsl direct solve
// functions

#include "poisson2DPeriodic.h"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

static unsigned Nx, Ny;
static double dx, dy;
static double Lx, Ly;
static gsl_vector b;
static const unsigned maxiter = 100;
static const double abstol = 1e-8;

void init_poisson2D(unsigned Nx_, unsigned Ny_, double dx_, double dy_,
                    double Lx_, double Ly_) {

  Nx = Nx_;
  Ny = Ny_;
  dx = dx_;
  dy = dy_;
  Lx = Lx_;
  Ly = Ly_;
}

// Quick n dirty Gauss seidel iterative solve
// Note - I probably need an initial guess better than 0??
// The nonlinear one before used 0 on first time step, then just used a previous
// phi
void GaussSeidel(gsl_matrix *A, gsl_vector *b, gsl_vector *sol) {
  unsigned i, j;

  gsl_matrix *L = gsl_matrix_calloc(Nx * Ny, Nx * Ny);
  gsl_matrix *U = gsl_matrix_calloc(Nx * Ny, Nx * Ny);
  gsl_vector *xold = gsl_vector_calloc(Nx * Ny);
  gsl_vector *xnew = gsl_vector_calloc(Nx * Ny);
  gsl_vector *RHS = gsl_vector_calloc(Nx * Ny);

  gsl_permutation *P =
      gsl_permutation_calloc(Nx * Ny); // This should set it to identity

  unsigned iter = 0;
  double residual;

  // initialize phi
  gsl_vector_dcopy(sol, xold);

  // Fill the matrices
  for (i = 0; i < Nx * Ny; i++) {
    for (j = 0; j < i; j++)
      gsl_matrix_set(L, i, j) = gsl_matrix_get(A, i, j);
    for (j = i; j < Nx * Ny; j++)
      gsl_matrix_set(U, i, j) = gsl_matrix_get(A, i, j);
  }

  for (iter = 0; iter < maxiters; ++iter) {

    // compute new RHS
    // RHS = b - Lx_old

    // Incomprehensible blas does RHS = -Lx_old
    gsl_blas_dgemv(CblasNoTrans, -1.0, L, xold, 0.0, RHS);
    // Incomprehensible blas does RHS = RHS + b
    gsl_blas_daxpy(1.0, b, RHS);

    // Now solve Lxnew = RHS
    gsl_linalg_LU_solve(U, P, RHS, xnew);

    // Check residual (re-use RHS data block, since no longer needed)
    // Compute A*xnew
    gsl_blas_dgemv(CblasNoTrans, 1.0, A, xnew, 0.0, RHS);
    // Compute A*xnew - b
    gsl_blas_daxpy(-1.0, b, RHS);

    residual = gsl_blas_dnrm2(RHS);
    if (residual < abstol)
      break;

    gsl_vector_dswap(xold, xnew);
  }

  printf("Number of iterations %d\n", iter);
  printf("Residual %le\n", residual);

  gsl_vector_dcopy(xnew, sol);

  gsl_matrix_free(L);
  gsl_matrix_free(U);
  gsl_vector_free(xold);
  gsl_vector_free(xnew);
  gsl_vector_free(RHS);
  gsl_vector_free(P);
}

// Generic Helmholtz solve - computes
//
// -\phi_xx - \phi_yy + alpha * phi = source
//
// on a periodic domain
//
// INPUTS
// alpha(x,y): the helmholtz function (linearization of electron dist func)
// source: the RHS of the equation (charge distribution)
//
// OUTPUTS
// phi(x,y): the solution
void PoissonLinearPeriodic2D(double **alpha, double **source, double **phi) {

  unsigned i, j;

  gsl_matrix *A = gsl_matrix_calloc(Nx * Ny, Nx * Ny);

  // Construct the A matrix using second order centered difference
  unsigned SELF, RIGHT, LEFT, UP, DOWN;

  for (i = 0; i < Nx; i++) {
    for (j = 0; j < Ny; j++) {
      SELF = j + i * Ny;
      UP = (j + 1) != Ny ? j : 0;
      DOWN = j == 0 ? Ny - 1 : j;
      RIGHT = (j + 1) != Ny ? j : 0;
      LEFT = j == 0 ? Ny - 1 : j;

      gsl_matrix_set(SELF, j + RIGHT * Ny, -1.0 / dx / dx);
      gsl_matrix_set(SELF, j + LEFT * Ny, -1.0 / dx / dx);

      gsl_matrix_set(SELF, SELF, 2.0 / dx / dx + 2.0 / dy / dy + alpha[i][j]);

      gsl_matrix_set(SELF, UP + i * Ny, -1.0 / dy / dy);
      gsl_matrix_set(SELF, DOWN + i * Ny, -1.0 / dy / dy);
    }
  }

  // Because this is linearized (i.e. helmholtz) we don't need a fudge for the
  // constant solution

  // Set up RHS
  gsl_vector *b = gsl_vector_calloc(Nx * Ny);
  for (i = 0; i < Nx; i++) {
    for (j = 0; j < Ny; j++) {
      gsl_vector_set(b, j + i * Ny, source[i][j]);
    }
  }

  gsl_vector *phi_v = gsl_vector_calloc(Nx * Ny);
  for (i = 0; i < Nx; i++) {
    for (j = 0; j < Ny; j++) {
      gsl_vector_set(phi_v, j + i * Ny, phi[i][j]);
    }
  }

  // Solve Ax = b by LU
  /*
  gsl_permutation P = gsl_permutation_calloc(Nx*Ny);
  int signum;
  gsl_linalg_LU_decomp(A, P, signum);
  gsl_linalg_LU_solve(A, P, b, psi_v)
  gsl_permutation_free(P);
  */

  // Now solve Ax = b using gauss-seidel
  GaussSeidel(A, b, phi_v);

  // Turn psi back to a bare array
  for (i = 0; i < Nx; i++)
    for (j = 0; j < Ny; j++)
      phi[i][j] = gsl_vector_get(phi_v, j + i * Ny);

  gsl_matrix_free(A);
  gsl_vector_free(b);
  gsl_vector_free(phi_v);
}

void PoissonNonlinearPeriodic2D(double **source, double **Te, double **phi) {
  printf(
      "Error: Nonlinear, classical Poisson solve not yet implemented in 2D\n");
  exit(1);
}
void PoissonLinearPeriodic2D_TF(double **source, double **Te, double **phi);
{
  printf(
      "Error: Linear, Thomas-Fermi Poisson solve not yet implemented in 2D\n");
  exit(1);
}
void PoissonNonlinearPeriodic2D_TF(double **source, double **Te, double **phi);
{
  printf("Error: Nonlinear, Thomas-Fermi Poisson solve not yet implemented in "
         "2D\n");
  exit(1);
}
