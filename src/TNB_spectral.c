#include "fourier.h"
#include "gauss_legendre.h"
#include "units/unit_data.c"
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef EPS_COLL
#define EPS_COLL 1e6
#endif

static int Nv;
static double **c;
static double **wts;
static double Lv;

static int first_DD = 1;
static int first_DT = 1;

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

  strcpy(DDHE.name, "DD-He3");

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

  strcpy(DDT.name, "DD-T");
}

// forward declaration
void generate_conv_weights(double *conv_weights,
                           struct TNB_data *reaction_info);

// Generic TNB calculator
// This takes input data for the various flavors of TNB reactions
// Not this assumes that you have already checked the species to make sure that
// it is TNB Computes Q_TNB(c) = \int |g| \sigma(|g|) f_1(c) f_2(c_\ast) d\c_ast
void TNB_generic(struct TNB_data *reaction_info, double *f1, double *f2, int sp,
                 int sp2, double *Q_TNB) {

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
    fread(conv_weights, sizeof(double), Nv * Nv * Nv, fidWeights);
  } else {
    generate_conv_weights(conv_weights, reaction_info);
    // dump the weights we've computed into a file

    fidWeights = fopen(weight_filename, "w");
    fwrite(conv_weights, sizeof(double), Nv * Nv * Nv, fidWeights);
    if (fflush(fidWeights) != 0) {
      printf("Something is wrong with storing the weights");
      exit(0);
    }
  }

  // More stuff to compute...
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
      }
    }
  }
}

// Calculates loss function due to DT TNB reactions for all velocity points
// Q_TNB(c) = \int g sigma_{DT}(g) f(c) f(v_\ast) d \v_ast
double GetTNB_dt(double mu, double *in, double *c1, int sp, int sp2);

// Calculates the total dt reaction rate in a cell
// R_DT = \int Q_TNB(c) dc
double GetReactivity_dt(double mu, double *in, double *in2, int sp, int sp2);

//-----------------------------------------//

// Calculates loss function due to dd->he TNB reactions for all velocity points
// Q_TNB(c) = \int g sigma_{DDHe}(g) f(c) f(v_\ast) d \v_ast
double GetTNB_dd_He(double mu, double *in, double *c1, int sp, int sp2);

// Calculates the total dd->he reaction rate in a cell
// R_DDHe = \int Q_TNB(c) dc
double GetReactivity_dd_He(double mu, double *in, double *in2, int sp, int sp2);

//-----------------------------------------//

// Calculates loss function due to dd->T TNB reactions for all velocity points
// Q_TNB(c) = \int g sigma_{DDT}(g) f(c) f(v_\ast) d \v_ast
double GetTNB_dd_T(double mu, double *in, double *c1, int sp, int sp2);

// Calculates the total dd->T reaction rate in a cell
// R_DDT = \int Q_TNB(c) dc
double GetReactivity_dd_T(double mu, double *in, double *in2, int sp, int sp2);

//-----------------------------------------//

// Calculates loss function due to TT TNB reactions for all velocity points
// Q_TNB(c) = \int g sigma_{TT}(g) f(c) f(v_\ast) d \v_ast
double GetTNB_tt(double mu, double *in, double *c1, int sp, int sp2);

// Calculates the total TT reaction rate in a cell
// R_TT = \int Q_TNB(c) dc
double GetReactivity_tt(double mu, double *in, double *in2, int sp, int sp2);

//-----------------------------------------//

// I think this does the actual tail depletion
void TNB_DD(double **f, double **f_out, int sp, int rank, int TNB_FLAG,
            double dt, double mu, double *n, double **v, double *T);

void TNB_DT(double **f, double **f_out, int sp, int sp2, int rank, int TNB_FLAG,
            double dt, double mu, double *n, double **v, double *T);
