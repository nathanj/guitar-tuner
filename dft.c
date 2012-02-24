/**
 * Computes the discrete fourier transform (DFT).
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <stdint.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include <fftw3.h>

void dft(double complex *in, double complex *out, size_t len)
{
	int k, n, N;
	double complex c, c2;
	double complex sum;
	double arg;

	N = len;

	for (k = 0; k < N; k++) {
		sum = 0;
		for (n = 0; n < N; n++) {
			c = cexp(-I * 2 * M_PI * k * n / N);
			sum += in[n] * c;
		}
		out[k] = sum;
	}
}

void idft(double complex *in, double complex *out, size_t len)
{
	int k, n, N;
	double complex c, c2;
	double complex sum;
	double arg;

	N = len;

	for (n = 0; n < N; n++) {
		sum = 0;
		for (k = 0; k < N; k++) {
			c = cexp(I * 2 * M_PI * k * n / N);
			sum += in[k] * c;
		}
		out[n] = sum / N;
	}
}




int main()
{
	double complex b = 3 + 3 * I;
	double complex c;
	double complex values[100];
	double complex result[100];
	double complex result2[100];
	double x1[100];
	double y1[100];
	int i;
	const int n = 20;
	fftw_complex *in, *out;
	fftw_plan p;

	in = fftw_malloc(sizeof(fftw_complex) * n);
	out = fftw_malloc(sizeof(fftw_complex) * n);
	p = fftw_plan_dft_r2c_1d(n, x1, out, FFTW_ESTIMATE);

	for (i = 0; i < n; i++) {
		values[i] = 0.6 * sin(i * 3 * 2 * M_PI / n)
			+ 0.8 * cos(i * 8 * 2 * M_PI / n);
		x1[i] = creal(values[i]);
		y1[i] = 0;
		in[i] = values[i];
		printf("%f\n", creal(values[i]));
	}

	//values[0] = 5;
	//values[1] = 7;
	//values[2] = 8;
	//values[3] = 7;
	//values[4] = 3;
	//values[5] = 3;
	//values[6] = 3;
	//values[7] = 4;

	fftw_execute(p);
	dft(values, result, n);

	for (i = 0; i < n; i++) {
		printf("result[%d] = (%f %f) (%f)\n", i, creal(result[i]),
		       cimag(result[i]), cabs(result[i]));
	}

	for (i = 0; i < n; i++) {
		printf("out[%d] = (%f %f) (%f)\n", i, creal(out[i]),
		       cimag(out[i]), cabs(out[i]));
	}

	//result[2] = 0;
	result[3] = 0;
	result[4] = 0;
	result[5] = 0;
	//result[6] = 0;

	idft(result, result2, n);

	for (i = 0; i < n; i++) {
		printf("result2[%d] = (%f %f) (%f)\n", i, creal(result2[i]),
		       cimag(result2[i]), cabs(result2[i]));
	}

	//for (i = 0; i < n; i++) {
	//	printf("(%f %f)\n", x1[i], y1[i]);
	//}

	printf("%f %f %f %f\n", cabs(b), creal(b), carg(b), cimag(b));

	c = cexp(0 - I * 2 * M_PI);
	printf("%f %f\n", creal(c), cimag(c));

	fftw_destroy_plan(p);
	fftw_free(in);
	fftw_free(out);

	return 0;
}

