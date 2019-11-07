#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cstdint>
#include <cstdlib>
#include "world/matlabfunctions.h"
#include "world/common.h"
#include "world/fft.h"
#include "world/dio.h"
#include "sndfile.h"

// from https://github.com/lazybeaver/xorshift

class XorShift64Star {
private:
	uint64_t state;

public:
	XorShift64Star(uint64_t seed) : state(seed) {
	}

	uint64_t next() {
		this->state ^= (this->state >> 12);
		this->state ^= (this->state << 25);
		this->state ^= (this->state >> 27);
		return this->state * 2685821657736338717;
	}

	double nextDouble() {
		return (double)this->next() / (double)UINT64_MAX;
	}
};

void exampleInterp1() {
	double x [] = { 0,    2,          6 };
	double y [] = {10,   20,          0 };
	double xi[] = { 0, 1, 2, 3, 4, 5, 6 };
	double yi[] = { 0, 0, 0, 0, 0, 0, 0 };
	//     result  10 15 20 15 10  5  0

	int x_length = sizeof(x) / sizeof(x[0]);
	int xi_length = sizeof(xi) / sizeof(xi[0]);
	interp1(x, y, x_length, xi, xi_length, yi);
	for (int i=0; i<xi_length; i++) {
		printf("%4.1f\n", yi[i]);
	}

	// Output:
	// 10.0
	// 15.0
	// 20.0
	// 15.0
	// 10.0
	//  5.0
	//  0.0
}

void exampleInterp1Q() {
	double x = 0;
	double shift = 2;
	double y [] = {10,   20,    0 };
	double xi[] = { 0, 1, 2, 3, 4 };
	double yi[] = { 0, 0, 0, 0, 0 };
	//     result  10 15 20 10  0

	int x_length = sizeof(y) / sizeof(y[0]);
	int xi_length = sizeof(xi) / sizeof(xi[0]);
	interp1Q(x, shift, y, x_length, xi, xi_length, yi);
	for (int i=0; i<xi_length; i++) {
		printf("%4.1f\n", yi[i]);
	}

	// Output:
	// 10.0
	// 15.0
	// 20.0
	// 10.0
	//  0.0
}

void exampleDCCorrection() {
	double input [7] = {0, 2, 1, 0, 1, 2, 8};
	double output[7] = {0, 0, 0, 0, 0, 0, 0};
	int current_f0 = 64; // index=4
	int fs = 256;
	int fft_size = 16;
	DCCorrection(input, current_f0, fs, fft_size, output);
	int output_length = sizeof(output) / sizeof(output[0]);
	for (int i=0; i<output_length; i++) {
		printf("%4.1f\n", output[i]);
	}
}

void exampleFFT() {
	int fft_size = 8;
	double y[8] = {1, -1, .5, -.5, .25, -.25, .125, 0};
	fft_complex y_spectrum[8];
	fft_plan forwardFFT = fft_plan_dft_r2c_1d(fft_size, y, y_spectrum, FFT_FORWARD);
	fft_execute(forwardFFT);
	for (int i=0; i<fft_size; i++) {
		printf("%+4.3f +%4.3f\n", y_spectrum[i][0], y_spectrum[i][1]);
	}
	printf("\n");

	double y2[8] = {0};
	fft_plan backwardFFT = fft_plan_dft_c2r_1d(fft_size, y_spectrum, y2, FFT_BACKWARD);
	fft_execute(backwardFFT);
	for (int i=0; i<fft_size; i++) {
		printf("%+4.3f\n", y2[i]);
	}
	printf("\n");
}

#define bufSize = 100

void exampleDio(const char *file, int noiseDB, bool reverse) {
	SF_INFO info;
	SNDFILE *f = sf_open(file, SFM_READ, &info);
	if (f == NULL) {
		fprintf(stderr, "failed to open wav file\n");
		return;
	}
	if (info.channels != 1) {
		fprintf(stderr, "invalid channel count\n");
		sf_close(f);
		return;
	}
	double *x = new double[info.frames];
	sf_count_t n = sf_readf_double(f, x, info.frames);
	sf_close(f);
	if (n != info.frames) {
		fprintf(stderr, "failed to read wav file\n");
		delete[] x;
		return;
	}

	if (reverse) {
		for (int i=0; i<info.frames/2; i++) {
			int j = info.frames - 1 - i;
			double v = x[j];
			x[j] = x[i];
			x[i] = v;
		}
	}

	double noiseAmp = pow(2.0, noiseDB);
	XorShift64Star r(1);
	for (int i=0; i<info.frames; i++) {
		double v = r.nextDouble();
		x[i] += (v - .5) * 2.0 * noiseAmp;
	}

	int fs = info.samplerate;
	int x_length = info.frames;
	DioOption option;
	InitializeDioOption(&option);
	int f0_length = GetSamplesForDIO(fs, x_length, option.frame_period);
	double *temporalPositions = new double[f0_length];
	double *f0 = new double[f0_length];

	Dio(x, x_length, fs, &option, temporalPositions, f0);

	for (int i=0; i<f0_length; i++) {
		printf("%05.3f: %05.1f\n", temporalPositions[i], f0[i]);
	}

	delete[] x;
	delete[] f0;
	delete[] temporalPositions;
}

int main(int argc, char *argv[]) {
	if (argc < 5) {
		printf("Usage: make-testdata <function name> <wav file> <noise [DB]> <reverse>\n");
		return 1;
	}
	const char *fn = argv[1];
	const char *wav = argv[2];
	int noiseDB = atoi(argv[3]);
	int reverse = atoi(argv[4]) != 0;
	if (strcmp(fn, "interp1") == 0) {
		exampleInterp1();
		return 0;
	}
	if (strcmp(fn, "interp1Q") == 0) {
		exampleInterp1Q();
		return 0;
	}
	if (strcmp(fn, "DCCorrection") == 0) {
		exampleDCCorrection();
		return 0;
	}
	if (strcmp(fn, "FFT") == 0) {
		exampleFFT();
		return 0;
	}
	if (strcmp(fn, "Dio") == 0) {
		exampleDio(wav, noiseDB, reverse);
		return 0;
	}
	printf("function name \"%s\" is undefined\n", fn);
	return 1;
}
