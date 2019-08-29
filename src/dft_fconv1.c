/**
 * @file dft_fconv1.c
 * @brief Computation of fast convolution with N/2-point DFT (type 1)
 * @author matyunin.d
 * @date 29.08.2019
 */

#include <complex.h>
#include <math.h>

/**
 * @brief Computation of fast convolution with N/2-point DFT (type 1)
 * @param x0 Pointer to N/2-point complex-valued sequence, x0 = DFT[ s(2k)+i*s(2k+1) ], s - N-point real-valued sequence
 * @param x1 Pointer to N-point complex-valued sequence, x1 = DFT[ m(k) ], s - N-point real-valued sequence
 * @param y Pointer to output N/2-point complex-valued sequence
 * @param N Length
 */
void dft_fconv1(const double complex *x0, const double complex *x1, double complex *y, unsigned int N)
{
	double complex w0, w1, w2;
	double complex h0, h1;

	y[0] = x0[0] * creal(x1[N/2]) - conj(x0[0]) * cimag(x1[N/2]);

	for (unsigned int k = 1; k < N/2; k++) {
		w0 = 0.5 * (1 - sin( 2 * M_PI * k / N));
		w1 = 0.5 * (1 + sin( 2 * M_PI * k / N));
		w2 = 0.5 * I * cos(2 * M_PI * k / N);

		h0 = w0 * x1[k] + w1 * conj(x1[N/2 - k]);
		h1 = w2 * (x1[k] - conj(x1[N/2 - k]));

		y[k] = x0[k] * h0 + conj(x0[N/2 - k]) * h1;
	}
}

