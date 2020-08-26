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

static const double scale3 = pow(1.0 / sqrt(2.0 * M_PI), 3.0);

static double *wtN;

static fftw_plan p_forward;
static fftw_plan p_backward;
static fftw_complex *temp;

static int inverse = 1;
static int noinverse = 0;

void initialize_fourier(int modes, double *vel) {

  N = modes;

  // Gather data for making the Fourier grid
  v = vel;
  Lv = -vel[0];
  dv = vel[1] - vel[0];

  k_arr = (double *)malloc(N * sizeof(double));

  dk = 2 * M_PI * N / dv;
  Lk = 0.5 * (N - 1) * dk;
  for (int i = 0; i < N; i++) {
    k_arr[i] = -Lk + i * dk;
  }

  wtN = malloc(N * sizeof(double));
  wtN[0] = 0.5;
  for (int i = 1; i < (N - 1); i++) {
    wtN[i] = 1.0;
  }
  wtN[N - 1] = 0.5;

  temp = fftw_malloc(N * N * N * sizeof(fftw_complex));

  // Set up fftw infrastructure
  p_forward =
      fftw_plan_dft_3d(N, N, N, temp, temp, FFTW_FORWARD, FFTW_ESTIMATE);
  p_backward =
      fftw_plan_dft_3d(N, N, N, temp, temp, FFTW_BACKWARD, FFTW_ESTIMATE);
}

// fftw requires some shenanigans to get grid alignment correct
void fft3D(fftw_complex *in, fftw_complex *out, int invert) {
  int i, j, k, index;
  double sum, prefactor, factor;
  double delta, L_start, L_end, sign, *varr;
  fftw_plan p;

  if (invert == noinverse) {
    delta = dv;
    L_start = Lk;
    L_end = Lv;
    varr = k_arr;
    sign = 1.0;
    p = p_forward;
  } else if (invert == inverse) {
    delta = dk;
    L_start = Lv;
    L_end = Lk;
    varr = v;
    sign = -1.0;
    p = p_backward;
  } else {
    printf(
        "Something weird happened in how a fourier transform was specified\n");
    exit(37);
  }
  prefactor = scale3 * delta * delta * delta;

  printf("delta %g Lstart %g Lend %g sign %g varr[0] %g\n", delta, L_start,
         L_end, sign, varr[0]);

  // shift the 'v' terms in the exponential to reflect our velocity domain
  for (index = 0; index < N * N * N; index++) {
    i = index / (N * N);
    j = (index - i * N * N) / N;
    k = index - N * (j + i * N);
    sum = sign * (double)(i + j + k) * L_start * delta;

    factor = prefactor * wtN[i] * wtN[j] * wtN[k];

    // dv correspond to the velocity space scaling - ensures that the FFT is
    // properly scaled since fftw does no scaling at all
    temp[index][0] =
        factor * (cos(sum) * in[index][0] - sin(sum) * in[index][1]);
    temp[index][1] =
        factor * (cos(sum) * in[index][1] + sin(sum) * in[index][0]);
  }

  // computes fft
  fftw_execute(p);

  // shifts the 'eta' terms to reflect our fourier domain
  for (index = 0; index < N * N * N; index++) {
    i = index / (N * N);
    j = (index - i * N * N) / N;
    k = index - N * (j + i * N);
    sum = sign * L_end * (varr[i] + varr[j] + varr[k]);

    out[index][0] = cos(sum) * temp[index][0] - sin(sum) * temp[index][1];
    out[index][1] = cos(sum) * temp[index][1] + sin(sum) * temp[index][0];
  }
}
