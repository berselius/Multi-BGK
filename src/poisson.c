// Utility code for poisson solver(s)

#include "poisson.h"
#include "poissonNonlinPeriodic.h"
#include <mpi.h>

static int rank;
static int numRanks;
static int nspec;
static int Nx_rank;
static int Nx, dx, Lx;
static int poissonFlavor;
static int order;
static int *Nx_ranks;

static MPI_Status status;

static double *source_buf;

// All ranks
void initialize_local_poisson(int nspec_, int Nx_rank_, int order_,
                              double **PoisPot_rank, double **source,
                              double **Te) {

  nspec = nspec_;
  Nx_rank = Nx_rank_;
  order = order_;

  *PoisPot_rank = malloc((Nx_rank + 2 * order) * sizeof(double));
  *source = malloc(Nx_rank * sizeof(double));
  *Te = malloc(Nx_rank * sizeof(double));
  source_buf = malloc(2 * (Nx_rank + 1) * sizeof(double));

  // Get MPI info
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
}

// Rank 0 function
void initialize_global_poisson(int Nx_, int dx_, int Lx_, int *Nx_ranks_,
                               int poissonFlavor_, double **PoisPot_allranks,
                               double **source_allranks,
                               double **Te_arr_allranks) {

  Nx = Nx_;
  dx = dx_;
  Lx = Lx_;

  poissonFlavor = poissonFlavor_;

  Nx_ranks = Nx_ranks_;

  *PoisPot_allranks = malloc(Nx * sizeof(double));
  *source_allranks = malloc(Nx * sizeof(double));
  *Te_arr_allranks = malloc(Nx * sizeof(double));
}

// All ranks
void calculate_local_charge_source(double **n, double **Z, double *source) {
  int i, l;

  for (l = 0; l < Nx_rank; l++) {
    source[l] = 0.0;
    for (i = 0; i < nspec; i++) {
      source[l] += Z[l][i] * n[l][i];
    }
  }
}

// Rank 0 function
void gather_charge_sources(double *source, double *Te, double *source_allranks,
                           double *Te_allranks) {

  int l;
  int rankOffset;
  int rankCounter;

  // First do your own...
  for (l = 0; l < Nx_rank; l++) {
    source_allranks[l] = source[l];
    Te_allranks[l] = Te[l];
  }

  if (numRanks > 1) {
    rankOffset = Nx_rank;

    // Blarg, I need Nx per rank to do this shit
    for (rankCounter = 1; rankCounter < numRanks; rankCounter++) {

      MPI_Recv(source_buf, 2 * Nx_ranks[rankCounter], MPI_DOUBLE, rankCounter,
               2000 + rankCounter, MPI_COMM_WORLD, &status);

      for (l = 0; l < Nx_ranks[rankCounter]; l++) {
        source_allranks[l + rankOffset] = source_buf[2 * l + 0];
        Te_allranks[l + rankOffset] = source_buf[2 * l + 1];
      }

      rankOffset += Nx_ranks[rankCounter];
    }
  }
}

// Rank > 0 function
void send_charge_sources(double *source, double *Te) {
  int l;

  for (l = 0; l < Nx_rank; l++) {
    source_buf[2 * l + 0] = source[l];
    source_buf[2 * l + 1] = Te[l];
  }

  MPI_Send(source_buf, 2 * Nx_rank, MPI_DOUBLE, 0, 2000 + rank, MPI_COMM_WORLD);
}

void poisson_solve(double *source_allranks, double *Te_allranks,
                   double *PoisPot_allranks) {

  int l;

  if (poissonFlavor == 0) { // no E-field
    for (l = 0; l < Nx; l++)
      PoisPot_allranks[l] = 0.0;
  } else if (poissonFlavor == 11) // Linear Yukawa
    PoissLinPeriodic1D(Nx, source_allranks, dx, Lx, PoisPot_allranks,
                       Te_allranks);
  else if (poissonFlavor == 12) // Nonlinear Yukawa
    PoissNonlinPeriodic1D(Nx, source_allranks, dx, Lx, PoisPot_allranks,
                          Te_allranks);
  else if (poissonFlavor == 21) // Linear Thomas-Fermi
    PoissLinPeriodic1D_TF(Nx, source_allranks, dx, Lx, PoisPot_allranks,
                          Te_allranks);
  else if (poissonFlavor == 22) // Nonlinear Thomas-Fermi
    PoissNonlinPeriodic1D_TF(Nx, source_allranks, dx, Lx, PoisPot_allranks,
                             Te_allranks);
  else {
    printf("ERROR: Please set your poisson solver option to 0, 11, 12, 21, or "
           "22\n");
    exit(1);
  }
}

// Rank 0 function
void send_poisson_potential(double *PoisPot_allranks, double *PoisPot) {

  int rankOffset = 0;
  int rankCounter;
  int l;

  // Set the local ghost cells
  if (order == 1) {
    PoisPot[0] = PoisPot_allranks[Nx - 1];
    if (numRanks == 1)
      PoisPot[Nx_rank + 1] = PoisPot_allranks[0];
    else
      PoisPot[Nx_rank + 1] = PoisPot_allranks[Nx_rank];
  } else if (order == 2) {
    PoisPot[0] = PoisPot_allranks[Nx - 2];
    PoisPot[1] = PoisPot_allranks[Nx - 1];
    if (numRanks == 1) {
      PoisPot[Nx_rank + 2] = PoisPot_allranks[0];
      PoisPot[Nx_rank + 3] = PoisPot_allranks[1];
    } else {
      PoisPot[Nx_rank + 2] = PoisPot_allranks[Nx_rank];
      PoisPot[Nx_rank + 3] = PoisPot_allranks[Nx_rank + 1];
    }
  }

  // Set main body of PoisPot
  for (l = 0; l < Nx_rank; l++) {
    PoisPot[l + order] = PoisPot_allranks[l];
  }

  // Now deal with other ranks, if needed
  if (numRanks > 1) {

    rankOffset = Nx_rank;

    for (rankCounter = 1; rankCounter < numRanks - 1; rankCounter++) {

      // Send local ghost cell info
      if (order == 1) {
        source_buf[0] = PoisPot_allranks[rankOffset - 1];
        source_buf[Nx_ranks[rankCounter] + 1] =
            PoisPot_allranks[rankOffset + Nx_ranks[rankCounter]];
      } else if (order == 2) {
        source_buf[0] = PoisPot_allranks[rankOffset - 2];
        source_buf[1] = PoisPot_allranks[rankOffset - 1];
        source_buf[Nx_ranks[rankCounter] + 2] =
            PoisPot_allranks[rankOffset + Nx_ranks[rankCounter]];
        source_buf[Nx_ranks[rankCounter] + 3] =
            PoisPot_allranks[rankOffset + Nx_ranks[rankCounter] + 1];
      }

      // Set main body of PoisPot on rankCounter
      for (l = 0; l < Nx_ranks[rankCounter]; l++) {
        source_buf[l + order] = PoisPot_allranks[rankOffset + l];
      }

      MPI_Send(source_buf, Nx_ranks[rankCounter] + 2 * order, MPI_DOUBLE,
               rankCounter, rankCounter, MPI_COMM_WORLD);

      rankOffset += Nx_ranks[rankCounter];
    }

    // Deal with periodic BC for rightmost rank

    // SET LOCAL GHOST CELLS
    if (order == 1) {
      source_buf[0] = PoisPot_allranks[rankOffset - 1];
      source_buf[Nx_ranks[numRanks - 1] + 1] = PoisPot_allranks[0];
    } else if (order == 2) {
      source_buf[0] = PoisPot_allranks[rankOffset - 2];
      source_buf[1] = PoisPot_allranks[rankOffset - 1];
      source_buf[Nx_ranks[numRanks - 1] + 2] = PoisPot_allranks[0];
      source_buf[Nx_ranks[numRanks - 1] + 3] = PoisPot_allranks[1];
    }

    for (l = 0; l < Nx_ranks[numRanks - 1]; l++) {
      source_buf[l + order] = PoisPot_allranks[rankOffset + l];
    }

    MPI_Send(source_buf, Nx_ranks[numRanks - 1] + 2 * order, MPI_DOUBLE,
             numRanks - 1, numRanks - 1, MPI_COMM_WORLD);
  }
}

// Ranks > 0
void recieve_poisson_potentials(double *PoisPot) {
  int l;

  MPI_Recv(source_buf, Nx_rank + 2 * order, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD,
           &status);

  // Set Poispot
  for (l = 0; l < Nx_rank + 2 * order; l++) {
    PoisPot[l] = source_buf[l];
  }
}
