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
#include <pthread.h>

#include <fftw3.h>
#include <alsa/asoundlib.h>

#define SAMPLE_SIZE     (4096*2)

snd_pcm_uframes_t frames = SAMPLE_SIZE;

unsigned char *buffer;

int init_alsa(snd_pcm_t **handle)
{
	snd_pcm_hw_params_t *params;
	unsigned int rate = 44100;
	unsigned int exact_rate = rate;
	unsigned int dir = 0;
	int rc;
	int size;

	/* Open default PCM device for recording. */
	rc = snd_pcm_open(handle, "default", SND_PCM_STREAM_CAPTURE, 0);

	if (rc < 0) {
		printf("Unable to open pcm device: %s\n", snd_strerror(rc));
		exit(1);
	} else {
		printf("Using pcm device: %s\n", snd_pcm_name(*handle));
	}

	/* Allocate a hardware parameter object and fill it with default
	 * values. */
	snd_pcm_hw_params_alloca(&params);
	snd_pcm_hw_params_any(*handle, params);

	/* Interleaved. */
	rc = snd_pcm_hw_params_set_access(*handle, params, SND_PCM_ACCESS_RW_INTERLEAVED);
	if (rc < 0) {
		printf("Unable to set interleaved mode\n");
		return 1;
	}

	/* Signed 16-bit little-endian format. */
	rc = snd_pcm_hw_params_set_format(*handle, params, SND_PCM_FORMAT_S16_LE);
	if (rc < 0) {
		printf("Unable to set format\n");
		return 1;
	}

	/* Mono. */
	rc = snd_pcm_hw_params_set_channels(*handle, params, 1);
	if (rc < 0) {
		printf("Unable to set to two channel stereo\n");
		return 1;
	}

	/* 44.1 Khz */
	rc = snd_pcm_hw_params_set_rate_near(*handle, params, &exact_rate, 0);
	if (rc < 0) {
		printf("Error setting rate\n");
		return 1;
	}

	if (rate != exact_rate) {
		printf("%d Hz not supported. Using %d Hz instead", rate, exact_rate);
	}

	/* Set period size to 32 frames. */
	rc = snd_pcm_hw_params_set_period_size_near(*handle, params, &frames, &dir);
	if (rc < 0) {
		printf("Error setting period size\n");
		return 1;
	}

	/* Write the parameters to the driver. */
	rc = snd_pcm_hw_params(*handle, params);
	if (rc < 0) {
		printf("Unable to set hw parameters: %s\n", snd_strerror(rc));
		return 1;
	}

	/* 16-bit, mono. */
	size = frames * 2;
	buffer = (unsigned char *) malloc(size);

	printf("Frames: %d\n", (int) frames);
	printf("Allocating buffer of size: %d\n", size);

	return 0;

}

static char print_char(char c)
{
	if (isprint(c))
		return c;
	else
		return '.';
}

void get_alsa(snd_pcm_t *handle, double *fft_input)
{
	int rc;

	rc = snd_pcm_readi(handle, buffer, frames);

	if (rc == -EPIPE) {
		printf("Overrun occurred\n");
		snd_pcm_prepare(handle);
	} else if (rc < 0) {
		printf("Error from read: %s\n", snd_strerror(rc));
	} else if (rc != (int) frames) {
		printf("Short read. Only read %d frames, expected %d\n", rc, (int) frames);
	} else if (rc == (int) frames) {
		int i;

		/* Convert to double. */
		for (i = 0; i < rc; i++) {
			int16_t left = *(int16_t *) &buffer[i * 4];

			fft_input[i] = (double) left;
		}
	}
}


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
	double x1[SAMPLE_SIZE];
	double y1[100];
	int i;
	const int n = 20;
	fftw_complex *in, *out;
	fftw_plan p;
	snd_pcm_t *handle;
	double *fft_input;


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

	init_alsa(&handle);
	fft_input = malloc(sizeof(double) * SAMPLE_SIZE);

	out = fftw_malloc(sizeof(fftw_complex) * SAMPLE_SIZE);
	p = fftw_plan_dft_r2c_1d(SAMPLE_SIZE, fft_input, out, FFTW_ESTIMATE);
	while (1) {
		int i;
		double max = 0;
		int max_i[10] = { 0 };

		get_alsa(handle, fft_input);

		fftw_execute(p);

		for (i = 0; i < SAMPLE_SIZE; i++) {
			double abs = cabs(out[i]);
			if (abs > max) {
				max = abs;
				max_i[9] = max_i[8];
				max_i[8] = max_i[7];
				max_i[7] = max_i[6];
				max_i[6] = max_i[5];
				max_i[5] = max_i[4];
				max_i[4] = max_i[3];
				max_i[3] = max_i[2];
				max_i[2] = max_i[1];
				max_i[1] = max_i[0];
				max_i[0] = i;
			}
		}

		printf("max was %d %d %d %d %d\n",
		       max_i[0],
		       max_i[1],
		       max_i[2],
		       max_i[3],
		       max_i[4]);

	}
	fftw_destroy_plan(p);
	fftw_free(out);

	return 0;
}

