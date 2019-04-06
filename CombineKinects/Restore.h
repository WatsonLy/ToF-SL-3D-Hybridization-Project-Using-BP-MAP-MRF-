#pragma once
#include <cstdio>

#include <iostream>

#include <algorithm>

#include <assert.h>

#include <cstring>

#include "Restore.h"

#include "image.h"

#include "misc.h"

#include "pnmfile.h"

#define ITER 5        // number of BP iterations at each scale
#define LEVELS 1     // number of scales

#define DISC_K 200.0F       // truncation of discontinuity cost
#define DATA_K 10000.0F     // truncation of data cost
#define LAMBDA1 0.05F         // weighting of first data cost
#define LAMBDA2 0.05F		// weighting of the second data cost


#define INF 1E10     // large cost
#define VALUES 4000   // number of possible graylevel values

class Restore
{
public:

	image<uint16_t> *restore_ms(image<uint16_t> *img1, image<uint16_t> *img2);
	Restore();
	~Restore();

private:

	void bp_cb(image<float[VALUES]> *u, image<float[VALUES]> *d,
		image<float[VALUES]> *l, image<float[VALUES]> *r,
		image<float[VALUES]> *data,
		int iter, image<uint16_t> *img1, image<uint16_t> *img2);

	image<uint16_t>* output(image<float[VALUES]> *u, image<float[VALUES]> *d,
		image<float[VALUES]> *l, image<float[VALUES]> *r,
		image<float[VALUES]> *data);

	image<float[VALUES]>* comp_data(image<uint16_t> *img1, image<uint16_t> *img2);

	void msg(float s1[VALUES], float s2[VALUES],
		float s3[VALUES], float s4[VALUES],
		float dst[VALUES], int &low, int &up);

	float* dt(float *f, int n, int &low, int &up);

};

