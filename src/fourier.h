// Utility functions for fourier transforms

#ifndef _FOURIER_H
#define _FOURIER_H

#include <fftw3.h>

void initialize_fourier(int modes, double *vel);

void fft1D(fftw_complex *in, fftw_complex *out, int invert);

void fft3D(fftw_complex *in, fftw_complex *out, int invert);

void fft1D_noshift(fftw_complex *in, fftw_complex *out, int invert);

void fft3D_noshift(fftw_complex *in, fftw_complex *out, int invert);

#endif
