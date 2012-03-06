/**
 * Computes the discrete fourier transform (DFT).
 */

#define WANT_SDL

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
#include <signal.h>

#include <fftw3.h>
#include <alsa/asoundlib.h>

#define SAMPLE_SIZE     (1024 * 8)

#ifdef WANT_SDL
#include <SDL.h>

#define SCREEN_WIDTH    640
#define SCREEN_HEIGHT   480
#define SCREEN_BPP      32
#define HISTOGRAM_BINS  500
#define BIN_WIDTH       1
#endif

struct note_frequency {
	const char *name;
	int freq;
	const char *shifted_name;
};

/*
 * I don't know why, but the values I am getting from the Fourier
 * transform are not what I expect. When I play a 'D' note, I get the
 * frequency for 'A'. I always get the frequency for the note two and a
 * half steps below what I expect. Until I figure it out, have a field
 * for the shifted name so that I display the correct name.
 */
struct note_frequency frequencies[] =
{
	{ "C",  262, "F",  },
	{ "C#", 277, "F#", },
	{ "D",  294, "G",  },
	{ "D#", 311, "G#", },
	{ "E",  330, "A",  },
	{ "F",  349, "A#", },
	{ "F#", 370, "B",  },
	{ "G",  392, "C",  },
	{ "G#", 415, "C#", },
	{ "A",  440, "D",  },
	{ "A#", 466, "D#", },
	{ "B",  494, "E",  },
};

#define MIN_FREQ 250
#define MAX_FREQ 500
#define NUM_FREQS (sizeof(frequencies) / sizeof(*frequencies))

snd_pcm_uframes_t frames = SAMPLE_SIZE;

unsigned char *buffer;

int quit = 0;

static void handler(int sig, siginfo_t *si, void *unused)
{
	printf("Got SIGINT, exiting...\n");
	quit = 1;
}

#ifdef WANT_SDL
int init_sdl(SDL_Surface **screen)
{
	assert(screen);

	if (SDL_Init(SDL_INIT_EVERYTHING) == -1) {
		fprintf(stderr, "SDL_Init failed.\n");
		return 1;
	}

	*screen = SDL_SetVideoMode(SCREEN_WIDTH, SCREEN_HEIGHT,
				  SCREEN_BPP,
				  SDL_SWSURFACE);

	if (*screen == NULL) {
		fprintf(stderr, "SDL_SetVideoMode failed.\n");
		return 1;
	}

	SDL_WM_SetCaption("Guitar Tuner", NULL);

	return 0;
}

void put_pixel(SDL_Surface *screen, int x, int y, Uint8 r, Uint8 g, Uint8 b)
{
	Uint32 *pixmem32;
	Uint32 color;
	int pitch;

	color = SDL_MapRGB(screen->format, r, g, b);

	pitch = y * screen->pitch / screen->format->BytesPerPixel;
	pixmem32 = (Uint32 *) screen->pixels + pitch + x;
	*pixmem32 = color;
}

int draw(SDL_Surface *screen)
{
	if (SDL_Flip(screen) == -1) {
		fprintf(stderr, "SDL_Flip failed.\n");
		return 1;
	}

	return 0;
}

int draw_hist(SDL_Surface *screen, int bin, double value)
{
	int y, starty;

	//printf("bin=%d, value=%f\n", bin, value);

	starty = SCREEN_HEIGHT - value;
	if (starty < 0)
		starty = 0;

	for (y = starty; y != SCREEN_HEIGHT; y++)
		put_pixel(screen, bin, y, abs(SCREEN_HEIGHT - y)/2, 255, abs(SCREEN_HEIGHT - y)/2);
}
#endif


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
	rc = snd_pcm_hw_params_malloc(&params);
	if (rc < 0) {
		printf("Unable to set allocate hardware params: %s\n",
		       snd_strerror(rc));
		return 1;
	}

	rc = snd_pcm_hw_params_any(*handle, params);
	if (rc < 0) {
		printf("Unable to initialize hardware params: %s\n",
		       snd_strerror(rc));
		return 1;
	}

	/* Interleaved. */
	rc = snd_pcm_hw_params_set_access(*handle, params,
					  SND_PCM_ACCESS_RW_INTERLEAVED);
	if (rc < 0) {
		printf("Unable to set interleaved mode: %s\n",
		       snd_strerror(rc));
		return 1;
	}

	/* Signed 16-bit little-endian format. */
	rc = snd_pcm_hw_params_set_format(*handle, params,
					  SND_PCM_FORMAT_S16_LE);
	if (rc < 0) {
		printf("Unable to set format: %s\n",
		       snd_strerror(rc));
		return 1;
	}

	/* Mono. */
	rc = snd_pcm_hw_params_set_channels(*handle, params, 1);
	if (rc < 0) {
		printf("Unable to set to mono: %s\n",
		       snd_strerror(rc));
		return 1;
	}

	/* 44.1 Khz */
	rc = snd_pcm_hw_params_set_rate_near(*handle, params, &exact_rate, 0);
	if (rc < 0) {
		printf("Error setting rate: %s\n",
		       snd_strerror(rc));
		return 1;
	}

	if (rate != exact_rate) {
		printf("%d Hz not supported. Using %d Hz instead", rate,
		       exact_rate);
	}

	/* Set period size to 32 frames. */
	rc = snd_pcm_hw_params_set_period_size_near(*handle, params, &frames,
						    &dir);
	if (rc < 0) {
		printf("Error setting period size: %s\n",
		       snd_strerror(rc));
		return 1;
	}

	/* Write the parameters to the driver. */
	rc = snd_pcm_hw_params(*handle, params);
	if (rc < 0) {
		printf("Unable to set hw parameters: %s\n",
		       snd_strerror(rc));
		return 1;
	}

	/* Prepare interface. */
	rc = snd_pcm_prepare(*handle);
	if (rc < 0) {
		printf("Unable to prepare audio interface: %s\n",
		       snd_strerror(rc));
		return 1;
	}

	snd_pcm_hw_params_free(params);

	/* 16-bit, mono. */
	size = frames * 2;
	buffer = malloc(size);

	printf("Frames: %d\n", (int) frames);
	printf("Allocating buffer of size: %d %d\n", size, exact_rate);

	return 0;

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
			int16_t left = *(int16_t *) &buffer[i * 2];

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

void print_closest_frequency(int val)
{
	struct note_frequency *closest = NULL;
	int diff, closest_diff = 9999;
	int i;

	if (val < 20) {
		printf("Frequency too low!\n");
		return;
	}

	while (val < MIN_FREQ)
		val *= 2;
	while (val > MAX_FREQ)
		val /= 2;

	for (i = 0; i < NUM_FREQS; i++) {
		diff = abs(frequencies[i].freq - val);
		if (diff < closest_diff) {
			closest_diff = diff;
			closest = &frequencies[i];
		}
	}

	if (closest)
		printf("=> %s <= (diff: %d)\n", closest->shifted_name,
		       closest_diff);
	else
		printf("closest = NULL!\n", closest);
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
	struct sigaction sa;
	int rc;
#ifdef WANT_SDL
	SDL_Surface *screen;
	SDL_Event event;
	double histogram[HISTOGRAM_BINS];
	int j;
#endif

	sa.sa_flags = SA_SIGINFO;
	sigemptyset(&sa.sa_mask);
	sa.sa_sigaction = handler;

	rc = sigaction(SIGINT, &sa, NULL);
	if (rc == -1) {
		perror("sigaction");
		exit(EXIT_FAILURE);
	}

#if 0

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
#endif

#ifdef WANT_SDL
	if (init_sdl(&screen) != 0) {
		fprintf(stderr, "init_sdl failed.\n");
		return 1;
	}
#endif

	init_alsa(&handle);
	fft_input = malloc(sizeof(double) * SAMPLE_SIZE);

	out = fftw_malloc(sizeof(fftw_complex) * SAMPLE_SIZE);
	p = fftw_plan_dft_r2c_1d(SAMPLE_SIZE, fft_input, out, FFTW_ESTIMATE);
	while (!quit) {
		int i;
		double max = 0;
		int max_i[10] = { 0 };

		get_alsa(handle, fft_input);

		fftw_execute(p);

		for (i = 0; i < SAMPLE_SIZE / 2; i++) {
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

		print_closest_frequency(max_i[0]);

#ifdef WANT_SDL
		SDL_FillRect(screen, NULL, SDL_MapRGB(screen->format, 0, 0, 0));

		for (i = 0; i < HISTOGRAM_BINS; i++) {
			histogram[i] = 0;
			for (j = 0; j < BIN_WIDTH; j++)
				histogram[i] += cabs(out[i * BIN_WIDTH + j]);
			histogram[i] /= 50000;
			histogram[i] /= BIN_WIDTH;

			draw_hist(screen, i, histogram[i]);
		}

		/*
		for (i = 0; i < HISTOGRAM_BINS; i++) {
			histogram[i] = 0;
			for (j = 0; j < BIN_WIDTH; j++)
				histogram[i] += fft_input[i * BIN_WIDTH + j];
			histogram[i] /= BIN_WIDTH;
			histogram[i] /= 5;

			draw_hist(screen, i, histogram[i]);
		}
		*/

		if (draw(screen) != 0) {
			fprintf(stderr, "draw failed\n");
			return 1;
		}

		while (SDL_PollEvent(&event)) {

			switch (event.type) {
			case SDL_QUIT:
				quit = 1;
				break;
			case SDL_KEYDOWN:
			case SDL_KEYUP:
				switch (event.key.keysym.sym) {
				case SDLK_q:
					quit = 1;
					break;
				default:
					break;
				}
				break;
			default:
				break;
			}
		}
#endif

	}
	fftw_destroy_plan(p);
	fftw_free(out);
	free(fft_input);

	snd_pcm_close(handle);
	free(buffer);

	SDL_Quit();

	return 0;
}

