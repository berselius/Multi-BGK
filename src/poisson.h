#ifndef POISSON_H
#define POISSON_H

// Basically, local ctor for poisson routines, deals with processor local data
// Call from each rank
void initialize_local_poisson(int nspec_, int Nx_rank_, int order_,
                              double **PoisPot_rank, double **source,
                              double **Te);

// The global ctor for poisson routines
// Call from Rank 0
void initialize_global_poisson(int Nx_, double dx_, double Lx_, int *Nx_ranks_,
                               int poissonFlavor_, double **PoisPot_allranks,
                               double **source_allranks,
                               double **Te_arr_allranks);

// Calculates n_e = \sum_i Z_i n_i in each cell
void calculate_local_charge_source(double **n, double **Z, double *source);

// Gathers source and electron temperature from all other ranks
// Call from rank 0
void gather_charge_sources(double *source, double *Te, double *source_allranks,
                           double *Te_allranks);

// Rank > 0 : send relevant data for global poisson solve
void send_charge_sources(double *source, double *Te);

// Computes the global poisson solve
// Call on rank 0
void poisson_solve(double *source_allranks, double *Te_allranks,
                   double *PoisPot_allranks);

// Sends results from global poisson solve back to other ranks
// Call from rank 0
void send_poisson_potential(double *PoisPot_allranks, double *PoisPot);

// Receives results from global poisson solve
// Call from ramk > 0
void recieve_poisson_potentials(double *PoisPot);

#endif
