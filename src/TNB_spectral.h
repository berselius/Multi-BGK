#ifndef _TNB_SPECTRAL_H
#define _TNB_SPECTRAL_H

// Initializer

void initializeTNB(int Nv_in, double **c_in, double **wts_in);

//-----------------------------------------//

// Calculates loss function due to DT TNB reactions for all velocity points
// Q_TNB(c) = \int g sigma_{DT}(g) f(c) f(v_\ast) d \v_ast
void GetTNB_dt(double mu, double *in_1, double *in2, double *Q_DT);

//-----------------------------------------//

// Calculates loss function due to dd->he TNB reactions for all velocity points
// Q_TNB(c) = \int g sigma_{DDHe}(g) f(c) f(v_\ast) d \v_ast
// Note - will have to do this twice to get the D and T depletions
void GetTNB_dd_He(double mu, double *in, double *Q_DDHE);

//-----------------------------------------//

// Calculates loss function due to dd->T TNB reactions for all velocity points
// Q_TNB(c) = \int g sigma_{DDT}(g) f(c) f(v_\ast) d \v_ast
void GetTNB_dd_T(double mu, double *in, double *Q_DDT);

//-----------------------------------------//

// Calculates loss function due to TT TNB reactions for all velocity points
// Q_TNB(c) = \int g sigma_{TT}(g) f(c) f(v_\ast) d \v_ast
// void GetTNB_tt(double mu, double *in, double *Q_DDT);

//-----------------------------------------//

// Calculates the total reaction rate in a cell given the Q_TNB function
// R = \int Q_TNB(c) dc
double GetReactivity(double *Q_in);

// This is the main interface from other code
// If TNB_FLAG is 1, it simply computes the reactivity
// If TNB_FLAG is 2, it computes reactivity and does tail depeletion
// In both cases it writes reactivity info to a file in Data/
void TNB_DD(double *f_D, double *fout_D, int rank, int TNB_FLAG, double dt,
            double mu, double n, double T);

void TNB_DT(double *f_D, double *f_T, double *fout_D, double *fout_T, int rank,
            int TNB_FLAG, double dt, double mu, double n1, double n2, double T1,
            double T2);

#endif
