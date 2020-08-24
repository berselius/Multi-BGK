#ifndef _TNB_SPECTRAL_H
#define _TNB_SPECTRAL_H

// Initializer

void initializeTNB(int Nv_in, double **c_in, double **wts_in);

//-----------------------------------------//

// Calculates loss function due to DT TNB reactions for all velocity points
// Q_TNB(c) = \int g sigma_{DT}(g) f(c) f(v_\ast) d \v_ast
void GetTNB_dt(double mu, double *in_1, double *in2, double *Q_DT);

// Calculates the total dt reaction rate in a cell
// R_DT = \int Q_TNB(c) dc
double GetReactivity_dt(double mu, double *in, double *in2);

//-----------------------------------------//

// Calculates loss function due to dd->he TNB reactions for all velocity points
// Q_TNB(c) = \int g sigma_{DDHe}(g) f(c) f(v_\ast) d \v_ast
// Note - will have to do this twice to get the D and T depletions
void GetTNB_dd_He(double mu, double *in, double *Q_DDHE);

// Calculates the total dd->he reaction rate in a cell
// R_DDHe = \int Q_TNB(c) dc
double GetReactivity_dd_He(double mu, double *in);

//-----------------------------------------//

// Calculates loss function due to dd->T TNB reactions for all velocity points
// Q_TNB(c) = \int g sigma_{DDT}(g) f(c) f(v_\ast) d \v_ast
void GetTNB_dd_T(double mu, double *in, double *Q_DDT);

// Calculates the total dd->T reaction rate in a cell
// R_DDT = \int Q_TNB(c) dc
double GetReactivity_dd_T(double mu, double *in);

//-----------------------------------------//

// Calculates loss function due to TT TNB reactions for all velocity points
// Q_TNB(c) = \int g sigma_{TT}(g) f(c) f(v_\ast) d \v_ast
// void GetTNB_tt(double mu, double *in, double *Q_DDT);

// Calculates the total TT reaction rate in a cell
// R_TT = \int Q_TNB(c) dc
// double GetReactivity_tt(double mu, double *in, double *in2);

//-----------------------------------------//

// I think this does the actual tail depletion
void TNB_DD(double **f, double **f_out, int sp, int rank, int TNB_FLAG,
            double dt, double mu, double *n, double **v, double *T);

void TNB_DT(double **f, double **f_out, int sp, int sp2, int rank, int TNB_FLAG,
            double dt, double mu, double *n, double **v, double *T);

#endif
