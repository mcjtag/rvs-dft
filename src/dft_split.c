/**
 * @file dft_split.c
 * @brief Computation of the DTF of two real sequences (Splitting)
 * @author matyunin.d
 * @date 28.08.2019
 */

#include <complex.h>

/**
 * @brief Split two real sequences
 * s(k) = DFT[ xr(k) + i*xi(k)] => DFT[xr(k)] and DFT[xi(k)]
 * @param s Pointer to complex-valued sequence with two real-valued sequences
 * @param xr Pointer to output 'xr' sequence
 * @param xy Pointer to output 'xy' sequence
 * @param N Length of the sequences
 */
void dft_split(const double complex *s, complex double *xr, complex double *xi, unsigned int N)
{
	xr[0] = creal(s[0]);
	xi[0] = cimag(s[0]);
	xr[N/2] = creal(s[N/2]);
	xi[N/2] = cimag(s[N/2]);

	for (unsigned int k = 1; k < N/2; k++) {
		xr[k] = (s[k] + conj(s[N-k])) / 2;
		xi[k] = (s[k] - conj(s[N-k])) / (2*I);
	}

	for (unsigned int k = 1; k < N/2; k++) {
		xr[N/2+k] = conj(xr[N/2-k]);
		xi[N/2+k] = conj(xi[N/2-k]);
	}
}


