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

void initializeTNB(int Nv_in, double **c_in, double **wts_in) {
  printf("Initializing TNB\n");
  Nv = Nv_in;
  c = c_in;
  wts = wts_in;
}

// forward declaration
void generate_conv_weights(double *conv_weights,
                           struct TNB_data *reaction_info);

// Generic TNB calculator
// This takes input data for the various flavors of TNB reactions
// Computes Q_TNB(c) = \int |g| \sigma(|g|) f_1(c) f_2(c_\ast) d\c_ast
void TNB_generic(struct TNB_data *reaction_info, double mu, double *f1,
                 double *f2, int sp, int sp2, double *Q_TNB) {

  // Initialize Q_TNB
  for (int i = 0; i < Nv * Nv * Nv; i++) {
    Q_TNB[i] = 0.0;
  }

  // initialize weight array
  double *conv_weights = malloc(Nv * Nv * Nv * sizeof(double));

  // check to see if we have the right species setup
  if ((mu - reaction_info->mu_reaction) / reaction_info->mu_reaction > 0.01) {
    return;
  }

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

  double Lv;
  double spmax;

  for (int i = 0; i < Nv; ++i) {
    for (int j = 0; j < Nv; ++j) {
      for (int k = 0; k < Nv; ++k) {

        // add gauss legendre
        double result = gauss_legendre(64, integrand, args, 0, Lv);
      }
    }
  }
}
