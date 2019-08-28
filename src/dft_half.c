/**
 * @file dft_half.c
 * @brief Computation of the N/2-point DFT of a N-point real sequence (Recovery)
 * @author matyunin.d
 * @date 28.08.2019
 */

#include <complex.h>
#include <math.h>

/**
 * @brief Recovery N-point sequence from N/2-point sequence
 * x(k) = N/2-point DFT[ s(2*k) + i*s(2*k+1) ] => y(k) = N-point DFT[ s(k) ]
 * s - N-point real-valued sequence
 * @param x Pointer to input array of N/2-point complex-valued sequence (length = N/2)
 * @param y Pointer to output recovered N-point complex-valued sequence (length = N)
 * @param N Length
 */
void dft_half(const double complex *x, double complex *y, unsigned int N)
{
	double complex w0, w1;

	w0 = (1 - I) / 2;
	w1 = (1 + I) / 2;
	y[0]= x[0] * w0 + conj(x[0]) * w1;

	for (unsigned int k = 1; k < N/2; k++) {
		w0 = (1 - I * cexp(-I * 2 * M_PI * k / N)) / 2;
		w1 = (1 + I * cexp(-I * 2 * M_PI * k / N)) / 2;
		y[k]= x[k] * w0 + conj(x[N/2-k]) * w1;
	}

	y[N/2] = creal(x[0]) - cimag(x[0]);
	for (unsigned int k = 1; k < N/2; k++)
		y[N-k] = conj(y[k]);
}


