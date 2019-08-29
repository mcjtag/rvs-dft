/**
 * @file dft_fconv2.c
 * @brief Computation of fast convolution with N/2-point DFT (type 2)
 * @author matyunin.d
 * @date 29.08.2019
 */

#include <complex.h>
#include <math.h>

/**
 * @brief Computation of fast convolution with N/2-point DFT (type 2)
 * @param x0 Pointer to N/2-point complex-valued sequence, x = DFT[ s(2k)+i*s(2k+1) ], s - N-point real-valued sequence
 * @param x1 Pointer to N/2-point complex-valued sequence, h = DFT[ m(2k)+i*m(2k+1) ], m - N-point real-valued sequence
 * @param y Pointer to output N/2-point complex-valued sequence
 * @param N Length
 */
void dft_fconv2(const double complex *x0, const double complex *x1, double complex *y, unsigned int N)
{
	double complex w0, w1;
	double complex h0, h1;

	y[0] = x0[0] * creal(x1[0]);

	for (unsigned int k = 1; k < N/2; k++) {
		w0 = 0.25 * (3.0 - cexp(-4.0 * I * M_PI * k / N));
		w1 = 0.25 * (1.0 + cexp(-4.0 * I * M_PI * k / N));

		h0 = w0 * x1[k] + w1 * conj(x1[N/2 - k]);
		h1 = w1 * (x1[k] - conj(x1[N/2 - k]));

		y[k] = x0[k] * h0 + conj(x0[N/2 - k]) * h1;
	}
}


