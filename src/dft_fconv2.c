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
 * @param x Pointer to N/2-point complex-valued sequence, x = DFT[ s(2k)+i*s(2k+1) ], s - N-point real-valued sequence
 * @param h Pointer to N/2-point complex-valued sequence, h = DFT[ m(2k)+i*m(2k+1) ], s - N-point real-valued sequence (impulse response)
 * @param y Pointer to output N/2-point complex-valued sequence
 * @param N Length
 */
void dft_fconv2(const double complex *x, const double complex *h, double complex *y, unsigned int N)
{
	double complex w0, w1;
	double complex h0, h1;

	y[0] = x[0] * creal(h[0]);

	for (unsigned int k = 1; k < N/2; k++) {
		w0 = 0.25 * (3.0 - cexp(-4.0 * I * M_PI * k / N));
		w1 = 0.25 * (1.0 + cexp(-4.0 * I * M_PI * k / N));

		h0 = w0 * h[k] + w1 * conj(h[N/2 - k]);
		h1 = w1 * (h[k] - conj(h[N/2 - k]));

		y[k] = x[k] * h0 + conj(x[N/2 - k]) * h1;
	}
}


