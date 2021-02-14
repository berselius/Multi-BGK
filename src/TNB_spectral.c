#include "fourier.h"
#include "gauss_legendre.h"
#include "units/unit_data.c"
#include <fftw3.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef EPS_COLL
#define EPS_COLL 1e6
#endif

static int Nv;
static double *c;
static double *wts;
static double Lv;

static int first_DD = 1;
static int first_DT = 1;

static fftw_complex *temp_fftIn, *temp_fftOut;
static fftw_complex *fftIn_g, *fftOut_g, *qHat;

struct TNB_data {
  double a1;
  double a2;
  double a3;
  double a4;
  double a5;

  double b1;
  double b2;
  double b3;
  double b4;

  double B_G;

  double mu_reaction; // This is the relative mass for the reaction - used to
                      // check if we have the right species coming in

  double mabs; // extra data member tacked on for use in integrator

  char name[16];
};

static struct TNB_data DT;
static struct TNB_data DDHE;
static struct TNB_data DDT;

void initializeTNB(int Nv_in, double *c_in, double *wts_in) {
  printf("Initializing TNB\n");
  Nv = Nv_in;
  c = c_in;
  wts = wts_in;

  initialize_fourier(Nv, c_in);

  // Set up data structures for the reactions

  // DT
  DT.a1 = 6.927e4;
  DT.a2 = 7.454e8;
  DT.a3 = 2.050e6;
  DT.a4 = 5.2002e4;
  DT.a5 = 0;
  DT.b1 = 6.38e1;
  DT.b2 = -9.95e-1;
  DT.b3 = 6.981e-5;
  DT.b4 = 1.728e-4;
  DT.B_G = 34.33827;

  DT.mu_reaction = 2.0053e-24;

  strcpy(DT.name, "DT");

  // DD-He3
  DDHE.a1 = 5.3701e4;
  DDHE.a2 = 3.3027e2;
  DDHE.a3 = -1.2706e-1;
  DDHE.a4 = 2.9327e-5;
  DDHE.a5 = -2.5151e-9;
  DDHE.b1 = 0.0;
  DDHE.b2 = 0.0;
  DDHE.b3 = 0.0;
  DDHE.b4 = 0.0;
  DDHE.B_G = 31.3970;

  DDHE.mu_reaction = 1.672e-14;

  strcpy(DDHE.name, "DD_He3");

  // DD-T
  DDT.a1 = 5.5576e4;
  DDT.a2 = 2.1054e2;
  DDT.a3 = -3.2638e-2;
  DDT.a4 = 1.4987e-6;
  DDT.a5 = 1.8181e-10;
  DDT.b1 = 0.0;
  DDT.b2 = 0.0;
  DDT.b3 = 0.0;
  DDT.b4 = 0.0;
  DDT.B_G = 31.3970;

  DDT.mu_reaction = 1.672e-14;

  strcpy(DDT.name, "DD_T");

  // allocate bins for ffts
  fftIn_g = fftw_malloc(Nv * Nv * Nv * sizeof(fftw_complex));
  fftOut_g = fftw_malloc(Nv * Nv * Nv * sizeof(fftw_complex));
  qHat = fftw_malloc(Nv * Nv * Nv * sizeof(fftw_complex));
  temp_fftIn = fftw_malloc(Nv * Nv * Nv * sizeof(fftw_complex));
  temp_fftOut = fftw_malloc(Nv * Nv * Nv * sizeof(fftw_complex));
}

// forward declaration
void generate_conv_weights(double *conv_weights,
                           struct TNB_data *reaction_info);

// Generic TNB calculator
// This takes input data for the various flavors of TNB reactions
// Not this assumes that you have already checked the species to make sure that
// it is TNB Computes Q_TNB(c) = \int |g| \sigma(|g|) f_1(c) f_2(c_\ast) d\c_ast
void TNB_generic(struct TNB_data *reaction_info, double *f1, double *f2,
                 double *Q_TNB) {

  // Initialize Q_TNB
  for (int i = 0; i < Nv * Nv * Nv; i++) {
    Q_TNB[i] = 0.0;
  }

  // initialize weight array
  double *conv_weights = malloc(Nv * Nv * Nv * sizeof(double));

  // Check to see if the weight is generated
  char weight_filename[256] = "Data/";
  char buffer[256];
  sprintf(buffer, "TNB_weight_");
  strcat(buffer, reaction_info->name);
  strcat(weight_filename, buffer);

  FILE *fidWeights;
  if ((fidWeights = fopen(weight_filename, "r"))) {
    printf("Loading weights from %s\n", weight_filename);
    fflush(stdout);
    fread(conv_weights, sizeof(double), Nv * Nv * Nv, fidWeights);
  } else {
    printf("Weights for %s not found, generating and storing in %s\n",
           reaction_info->name, weight_filename);
    fflush(stdout);
    generate_conv_weights(conv_weights, reaction_info);
    // dump the weights we've computed into a file

    fidWeights = fopen(weight_filename, "w");
    fwrite(conv_weights, sizeof(double), Nv * Nv * Nv, fidWeights);
    printf("Weights stored for %s\n", reaction_info->name);
    if (fflush(fidWeights) != 0) {
      printf("Something is wrong with storing the weights");
      exit(0);
    }
  }

  // Weights are loaded, now compute the ffts of the functions

  // Fill fft data structions
  for (int index = 0; index < Nv * Nv * Nv; ++index) {
    fftIn_g[index][0] = f2[index];
    fftIn_g[index][1] = 0.0;
  }

  // move to fourier space - N^3 log N
  fft3D(fftIn_g, fftOut_g, 0);

  // Find inverse of W(m) ghat(m)
  // note W(m) is a real function
  // N^3
  for (int index = 0; index < Nv * Nv * Nv; ++index) {
    temp_fftIn[index][0] = fftOut_g[index][0] * conv_weights[index];
    temp_fftIn[index][1] = fftOut_g[index][1] * conv_weights[index];
  }

  fft3D(temp_fftIn, temp_fftOut, 1); // N^3 log N

  // Take product in real space
  // Just return the real part, but check on the imag
  // N^3
  double imagmax = 0;
  double imag;
  for (int index = 0; index < Nv * Nv * Nv; ++index) {
    Q_TNB[index] = f1[index] * temp_fftOut[index][0];
    imag = fabs(temp_fftOut[index][1]);
    imagmax = imag > imagmax ? imag : imagmax;
  }
}

double integrand(double r, void *args) {

  struct TNB_data *data = (struct TNB_data *)args;

  double E_COM = 0.5 * data->mu_reaction * r * r * ERG_TO_EV_CGS *
                 1e-3; // Center of mass energy in keV

  double nuclear_factor_num =
      data->a1 +
      E_COM * (data->a2 +
               E_COM * (data->a3 + E_COM * (data->a4 + E_COM * data->a5)));
  double nuclear_factor_denom =
      1.0 +
      E_COM * (data->b1 +
               E_COM * (data->b2 + E_COM * (data->b3 + E_COM * data->b4)));

  double xsec = (E_COM != 0) ? nuclear_factor_num / nuclear_factor_denom *
                                   exp(-data->B_G / sqrt(E_COM)) / E_COM
                             : 0;

  double sincval =
      (r * data->mabs < 1.0e-6) ? 1.0 : sin(r * data->mabs) / (r * data->mabs);

  return pow(r, 3.0) * xsec * sincval;
}

// Computes the convolution weights given the cross section parameters
void generate_conv_weights(double *conv_weights,
                           struct TNB_data *reaction_info) {
  for (int i = 0; i < Nv; ++i) {
    for (int j = 0; j < Nv; ++j) {
      for (int k = 0; k < Nv; ++k) {
        // add gauss legendre
        double result = gauss_legendre(64, integrand, reaction_info, 0, Lv);
        conv_weights[k + Nv * (j + Nv * i)] = result;
      }
    }
  }
}

// Calculates the total TNB reaction rate in a cell
// R_DT = \int Q_TNB(c) dc
double GetReactivity(double *Q_TNB) {
  int i, j, k;

  double react = 0;
#pragma omp parallel for private(i, j, k) reduction(+ : react)
  for (i = 0; i < Nv; i++) {
    for (j = 0; j < Nv; j++) {
      for (k = 0; k < Nv; k++) {
        react += wts[i] * wts[j] * wts[k] * Q_TNB[k + Nv * (j + Nv * i)];
      }
    }
  }
  return react;
}

// Calculates loss function due to DT TNB reactions for all velocity points
// Q_TNB(c) = \int g sigma_{DT}(g) f(c) f(v_\ast) d \v_ast
void GetTNB_dt(double mu, double *in_D, double *in_T, double *Q_DT) {

  // Check that we have the right reaction
  if (abs(mu - DT.mu_reaction) / DT.mu_reaction > 1e-3) {
    // this is not a DT reaction
    printf("Warning - this is not a DT reaction - given mu=%g, DT mu is %g\n",
           mu, DT.mu_reaction);
    return;
  }
  TNB_generic(&DT, in_D, in_T, Q_DT);
}

//-----------------------------------------//

// Calculates loss function due to dd->he TNB reactions for all velocity points
// Q_TNB(c) = \int g sigma_{DDHe}(g) f(c) f(v_\ast) d \v_ast
void GetTNB_dd_He(double mu, double *in, double *Q_DDHE) {
  // Check that we have the right reaction
  if (abs(mu - DDHE.mu_reaction) / DDHE.mu_reaction > 1e-3) {
    // this is not a DT reaction
    printf("Warning - This is not a DD reaction - given mu=%g, DD mu is %g\n",
           mu, DDHE.mu_reaction);
    return;
  }
  TNB_generic(&DDHE, in, in, Q_DDHE);
}

//-----------------------------------------//

// Calculates loss function due to dd->T TNB reactions for all velocity points
// Q_TNB(c) = \int g sigma_{DDT}(g) f(c) f(v_\ast) d \v_ast
void GetTNB_dd_T(double mu, double *in, double *Q_DDT) {
  // Check that we have the right reaction
  if (abs(mu - DDT.mu_reaction) / DDT.mu_reaction > 1e-3) {
    // this is not a DT reaction
    printf("Warning - this is not a DD reaction - given mu=%g, DD mu is %g\n",
           mu, DDT.mu_reaction);
    return;
  }
  TNB_generic(&DDT, in, in, Q_DDT);
}

//-----------------------------------------//

// Calculates loss function due to TT TNB reactions for all velocity points
// Q_TNB(c) = \int g sigma_{TT}(g) f(c) f(v_\ast) d \v_ast
// double GetTNB_tt(double mu, double *in, double *c1, int sp, int sp2) {}

// Calculates the total TT reaction rate in a cell
// R_TT = \int Q_TNB(c) dc
// double GetReactivity_tt(double mu, double *in, double *in2, int sp, int sp2)
// {}

//-----------------------------------------//

// This is the main interface from other code
// If TNB_FLAG is 1, it simply computes the reactivity
// IF TNB_FLAG is 2, it computes reactivity and does tail depeletion
void TNB_DD(double *f_D, double *fout_D, int rank, int TNB_FLAG, double dt,
            double mu, double n, double *v, double T) {

  char buffer[50];
  FILE *fpij;
  sprintf(buffer, "Data/TNB_DD_%d.dat", rank);

  double *Q_TNB_HE = malloc(Nv * Nv * Nv * sizeof(double));
  double *Q_TNB_T = malloc(Nv * Nv * Nv * sizeof(double));

  if (n > EPS_COLL) {

    // Calculate Q_TNB
    GetTNB_dd_He(mu, f_D, Q_TNB_HE);
    GetTNB_dd_T(mu, f_D, Q_TNB_T);

    // Get the reactivities
    double He_react = GetReactivity(Q_TNB_HE);
    double T_react = GetReactivity(Q_TNB_T);

    // Do tail depletion
    if (TNB_FLAG == 2) {
      for (int vx = 0; vx < Nv; vx++)
        for (int vy = 0; vy < Nv; vy++)
          for (int vz = 0; vz < Nv; vz++) {
            int index = vz + Nv * (vy + Nv * vx);
            fout_D[index] -=
                f_D[index] * (Q_TNB_HE[index] + Q_TNB_T[index]) * dt;
          }
    }

    if (first_DD) {
      fpij = fopen(buffer, "w");
      first_DD = 0;
    } else
      fpij = fopen(buffer, "a");

    fprintf(fpij, "%5.2e %5.2e %10.6e %10.6e\n", n, T, He_react, T_react);
    fclose(fpij);
  } else {
    if (first_DD) {
      fpij = fopen(buffer, "w");
      first_DD = 0;
    } else
      fpij = fopen(buffer, "a");
    fprintf(fpij, "%5.2e %5.2e %10.6e %10.6e\n", n, T, 0.0, 0.0);
    fclose(fpij);
  }

  free(Q_TNB_HE);
  free(Q_TNB_T);
}

void TNB_DT(double *f_D, double *f_T, double *fout_D, double *fout_T, int rank,
            int TNB_FLAG, double dt, double mu, double n1, double n2, double T1,
            double T2) {

  char buffer[50];
  FILE *fpij;
  sprintf(buffer, "Data/TNB_DT_%d.dat", rank);

  double *Q_TNB_D = malloc(Nv * Nv * Nv * sizeof(double));
  double *Q_TNB_T = malloc(Nv * Nv * Nv * sizeof(double));

  printf("DT TNB n:%g %g T:%g %g\n", n1, n2, T1, T2);

  if ((n1 > EPS_COLL) && (n2 > EPS_COLL)) {

    // Calculate Q_TNB
    GetTNB_dt(mu, f_D, f_T, Q_TNB_D);
    GetTNB_dt(mu, f_T, f_D, Q_TNB_T);

    // Get the reactivities
    double DT_react = GetReactivity(Q_TNB_D);

    // Do tail depletion
    if (TNB_FLAG == 2) {
      for (int vx = 0; vx < Nv; vx++)
        for (int vy = 0; vy < Nv; vy++)
          for (int vz = 0; vz < Nv; vz++) {
            int index = vz + Nv * (vy + Nv * vx);
            fout_D[index] -= f_D[index] * (Q_TNB_D[index]) * dt;
            fout_T[index] -= f_T[index] * (Q_TNB_T[index]) * dt;
          }
    }

    if (first_DT) {
      fpij = fopen(buffer, "w");
      first_DT = 0;
    } else
      fpij = fopen(buffer, "a");

    fprintf(fpij, "%5.2e %5.2e %5.2e %5.2e %10.6e \n", n1, n2, T1, T2,
            DT_react);
    fclose(fpij);
  } else {
    if (first_DT) {
      fpij = fopen(buffer, "w");
      first_DT = 0;
    } else
      fpij = fopen(buffer, "a");
    fprintf(fpij, "%5.2e %5.2e %5.2e %5.2e %10.6e \n", n1, n2, T1, T2, 0.0);
    fclose(fpij);
  }

  free(Q_TNB_D);
  free(Q_TNB_T);
}
