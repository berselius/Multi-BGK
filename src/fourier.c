// Package for computing fourier transforms of velocity grid data

#include <fftw3.h>
#include <math.h>
#include <stdlib.h>

static int N; // number of fourier modes, number of 1d grid points

static double *k_arr; // Fourier grid info
static double Lk;
static double dk;

static double *v; // velocity grid info
static double Lv;
static double dv;

static const double scale3 = pow(1.0 / (2.0 * M_PI), 1.5);
static const double scale1 = pow(1.0 / (2.0 * M_PI), 0.5);

static double *wtN;

void initialize_fourier(int modes, double *vel) {

  N = modes;

  // Gather data for making the Fourier grid
  v = vel;
  Lv = -vel[0];
  dv = vel[1] - vel[0];

  k_arr = (double *)malloc(N * sizeof(double));

  dk = (2 * M_PI / N) / dv;
  Lk = 0.5 * N * dk;
  for (int i = 0; i < N; i++) {
    k_arr[i] = -Lk + i * dk;
  }

  wtN = malloc(N * sizeof(double));
  // wtN[0] = 0.5;
  wtN[0] = 1.0;
  for (int i = 1; i < (N - 1); i++) {
    wtN[i] = 1.0;
  }
  // wtN[N - 1] = 0.5;
  wtN[N - 1] = 1.0;
}

// fftw requires some shenanigans to get grid alignment correct
void fft1D(fftw_complex *in, fftw_complex *out, int sign) {
  int i;
  double sum, prefactor, wt_factor;
  double delta, L_start, L_end, *varr;
  fftw_plan p;

  fftw_complex *temp1, *temp2;
  temp1 = fftw_malloc(N * sizeof(fftw_complex));
  temp2 = fftw_malloc(N * sizeof(fftw_complex));

  fftw_plan p_forward;
  fftw_plan p_backward;

  p_forward = fftw_plan_dft_1d(N, temp1, temp2, FFTW_FORWARD, FFTW_ESTIMATE);
  p_backward = fftw_plan_dft_1d(N, temp1, temp2, FFTW_BACKWARD, FFTW_ESTIMATE);

  if (sign == 1) {
    delta = dv;
    L_start = Lk;
    L_end = Lv;
    varr = k_arr;
    p = p_forward;
  } else if (sign == -1) {
    delta = dk;
    L_start = Lv;
    L_end = Lk;
    varr = v;
    p = p_backward;
  } else {
    printf(
        "Something weird happened in how a fourier transform was specified\n");
    exit(37);
  }

  // shift the 'v' terms in the exponential to reflect our velocity domain
  for (i = 0; i < N; i++) {
    sum = sign * (double)i * L_start * delta;

    wt_factor = wtN[i];

    // dv correspond to the velocity space scaling - ensures that the FFT is
    // properly scaled since fftw does no scaling at all
    temp1[i][0] = wt_factor * (cos(sum) * in[i][0] - sin(sum) * in[i][1]);
    temp1[i][1] = wt_factor * (cos(sum) * in[i][1] + sin(sum) * in[i][0]);
  }

  // computes fft on temp1, stores in temp2
  fftw_execute(p);

  prefactor = scale1 * delta;

  // shifts the 'eta' terms to reflect our fourier domain
  for (i = 0; i < N; i++) {
    sum = sign * L_end * varr[i];

    out[i][0] = prefactor * (cos(sum) * temp2[i][0] - sin(sum) * temp2[i][1]);
    out[i][1] = prefactor * (cos(sum) * temp2[i][1] + sin(sum) * temp2[i][0]);
  }

  fftw_free(temp1);
  fftw_free(temp2);
}

// fftw requires some shenanigans to get grid alignment correct
void fft3D(fftw_complex *in, fftw_complex *out, int sign) {
  int i, j, k, index;
  double sum, prefactor, wt_factor;
  double delta, L_start, L_end, *varr;
  fftw_plan p;

  fftw_complex *temp1, *temp2;
  temp1 = fftw_malloc(N * N * N * sizeof(fftw_complex));
  temp2 = fftw_malloc(N * N * N * sizeof(fftw_complex));

  fftw_plan p_forward;
  fftw_plan p_backward;

  p_forward =
      fftw_plan_dft_3d(N, N, N, temp1, temp2, FFTW_FORWARD, FFTW_ESTIMATE);
  p_backward =
      fftw_plan_dft_3d(N, N, N, temp1, temp2, FFTW_BACKWARD, FFTW_ESTIMATE);

  if (sign == 1) {
    delta = dv;
    L_start = Lk;
    L_end = Lv;
    varr = k_arr;
    p = p_forward;
  } else if (sign == -1) {
    delta = dk;
    L_start = Lv;
    L_end = Lk;
    varr = v;
    p = p_backward;
  } else {
    printf(
        "Something weird happened in how a fourier transform was specified\n");
    exit(37);
  }

  // shift the 'v' terms in the exponential to reflect our velocity domain
  for (index = 0; index < N * N * N; index++) {
    i = index / (N * N);
    j = (index - i * N * N) / N;
    k = index - N * (j + i * N);
    sum = sign * (double)(i + j + k) * L_start * delta;

    wt_factor = wtN[i] * wtN[j] * wtN[k];

    // dv correspond to the velocity space scaling - ensures that the FFT is
    // properly scaled since fftw does no scaling at all
    temp1[index][0] =
        wt_factor * (cos(sum) * in[index][0] - sin(sum) * in[index][1]);
    temp1[index][1] =
        wt_factor * (cos(sum) * in[index][1] + sin(sum) * in[index][0]);
  }

  // computes fft on temp1, stores in temp2
  fftw_execute(p);

  prefactor = scale3 * delta * delta * delta;

  // shifts the 'eta' terms to reflect our fourier domain
  for (index = 0; index < N * N * N; index++) {
    i = index / (N * N);
    j = (index - i * N * N) / N;
    k = index - N * (j + i * N);
    sum = sign * L_end * (varr[i] + varr[j] + varr[k]);

    out[index][0] =
        prefactor * (cos(sum) * temp2[index][0] - sin(sum) * temp2[index][1]);
    out[index][1] =
        prefactor * (cos(sum) * temp2[index][1] + sin(sum) * temp2[index][0]);
  }

  fftw_free(temp1);
  fftw_free(temp2);
}

void fft1D_noshift(fftw_complex *in, fftw_complex *out, int sign) {

  fftw_plan p;
  double rescale;

  if (sign == 1) {
    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    rescale = 1.0;
    // p = p_forward;
  } else if (sign == -1) {
    p = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    rescale = 1.0 / N;
    // p = p_backward;
  } else {
    printf("Incorrect sign set in fft call!\n");
    exit(37);
  }

  fftw_execute(p);

  for (int index = 0; index < N; index++) {
    out[index][0] *= rescale;
    out[index][1] *= rescale;
  }
}

void fft3D_noshift(fftw_complex *in, fftw_complex *out, int sign) {

  fftw_plan p;
  double N3;

  if (sign == 1) {
    p = fftw_plan_dft_3d(N, N, N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    N3 = 1.0;
    // p = p_forward;
  } else if (sign == -1) {
    p = fftw_plan_dft_3d(N, N, N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    N3 = 1.0 / (N * N * N);
    // p = p_backward;
  } else {
    printf("Incorrect sign set in fft call!\n");
    exit(37);
  }

  fftw_execute(p);

  for (int index = 0; index < N * N * N; index++) {
    out[index][0] *= N3;
    out[index][1] *= N3;
  }
}
