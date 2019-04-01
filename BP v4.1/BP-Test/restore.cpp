/*
Copyright (C) 2006 Pedro Felzenszwalb

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
*/


#include <cstdio>

#include <iostream>

#include <algorithm>

#include <assert.h>

#include <cstring>
#include <vector>

#include "image.h"

#include "misc.h"

#include "pnmfile.h"


#define ITER 1       // number of BP iterations at each scale
#define ITER0 2
#define LEVELS 1     // number of scales


#define DISC_K 200.0F       // truncation of discontinuity cost
#define DATA_K 10000.0F     // truncation of data cost
#define LAMBDA 0.05F         // weighting of data cost


#define INF 1E10     // large cost

#define VALUES 256   // number of possible graylevel values

#define TWIDTH 32	// Width of tile
#define THEIGHT 32	// Height of tile


// dt of 1d function

static float *dt(float *f, int n) {
  
	float *d = new float[n];
  
	int *v = new int[n];
  
	float *z = new float[n+1];
  
	int k = 0;
  
	v[0] = 0;
  
	z[0] = -INF;
  
	z[1] = +INF;

  
	for (int q = 1; q <= n-1; q++) 
	   {
    float s  = ((f[q]+square(q))-(f[v[k]]+square(v[k]))) / (2*(q-v[k]));
    
		while (s <= z[k]) 
		{
      k--;
      
			s  = ((f[q]+square(q))-(f[v[k]]+square(v[k]))) / (2*(q-v[k]));
    
		}
    
		k++;
    
		v[k] = q;
    
		z[k] = s;
    
		z[k+1] = +INF;
  
	   }
  
	k = 0;
  
	for (int q = 0; q <= n-1; q++) 
	{
    while (z[k+1] < q)
      
		k++;
    
	     d[q] = square(q-v[k]) + f[v[k]];
  
	}
  
	delete [] v;
  
	delete [] z;
  
	return d;

}


// compute message
void msg(float s1[VALUES], float s2[VALUES], float s3[VALUES], float s4[VALUES],
 float dst[VALUES]) {
 
// aggregate and find min
float minimum = INF;
  
for (int value = 0; value < VALUES; value++) {
	dst[value] = s1[value] + s2[value] + s3[value] + s4[value];
    if (dst[value] < minimum)
      minimum = dst[value];
  }

// dt
  float *tmp = dt(dst, VALUES);

  // truncate and store in destination vector
  minimum += DISC_K;
  for (int value = 0; value < VALUES; value++)
    dst[value] = std::min(tmp[value], minimum);

  // normalize
  float val = 0;
  for (int value = 0; value < VALUES; value++) 
    val += dst[value];

  val /= VALUES;
  for (int value = 0; value < VALUES; value++) 
    dst[value] -= val;

  delete tmp;
}

void msg1(float s1[VALUES], float s2[VALUES], float s3[VALUES], 
	float dst[VALUES]) {

	// aggregate and find min
	float minimum = INF;

	for (int value = 0; value < VALUES; value++) {
		dst[value] = s1[value] + s2[value] + s3[value];
		if (dst[value] < minimum)
			minimum = dst[value];
	}

	// dt
	float *tmp = dt(dst, VALUES);

	// truncate and store in destination vector
	minimum += DISC_K;
	for (int value = 0; value < VALUES; value++)
		dst[value] = std::min(tmp[value], minimum);

	// normalize
	float val = 0;
	for (int value = 0; value < VALUES; value++)
		val += dst[value];

	val /= VALUES;
	for (int value = 0; value < VALUES; value++)
		dst[value] -= val;

	delete tmp;
}

void msg2(float s1[VALUES], float s2[VALUES], float dst[VALUES]) {

	// aggregate and find min
	float minimum = INF;

	for (int value = 0; value < VALUES; value++) {
		dst[value] = s1[value] + s2[value];
		if (dst[value] < minimum)
			minimum = dst[value];
	}

	// dt
	float *tmp = dt(dst, VALUES);

	// truncate and store in destination vector
	minimum += DISC_K;
	for (int value = 0; value < VALUES; value++)
		dst[value] = std::min(tmp[value], minimum);

	// normalize
	float val = 0;
	for (int value = 0; value < VALUES; value++)
		val += dst[value];

	val /= VALUES;
	for (int value = 0; value < VALUES; value++)
		dst[value] -= val;

	delete tmp;
}

// computation of data costs
image<float[VALUES]> *comp_data(image<uchar> *img) {
  int width = img->width();
  int height = img->height();
  image<float[VALUES]> *data = new image<float[VALUES]>(width, height);

  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      for (int value = 0; value < VALUES; value++) {
	float val = square((float)(imRef(img, x, y)-value));
	imRef(data, x, y)[value] = LAMBDA * std::min(val, DATA_K);
      }
    }
  }

  return data;
}

// generate output from current messages
image<uchar> *output(image<float[VALUES]> *u, image<float[VALUES]> *d, 
		     image<float[VALUES]> *l, image<float[VALUES]> *r, 
		     image<float[VALUES]> *data) {
  int width = data->width();
  int height = data->height();
  image<uchar> *out = new image<uchar>(width, height);

  for (int y = 1; y < height-1; y++) {
    for (int x = 1; x < width-1; x++) {
      // keep track of best value for current pixel
      int best = 0;
      float best_val = INF;
      for (int value = 0; value < VALUES; value++) {
		float val = 
		  imRef(u, x, y+1)[value] +
		  imRef(d, x, y-1)[value] +
		  imRef(l, x+1, y)[value] +
		  imRef(r, x-1, y)[value] +
		  imRef(data, x, y)[value];
		if (val < best_val) {
		  best_val = val;
		  best = value;
		}
      }
      imRef(out, x, y) = best;
    }
  }

  return out;
}

// belief propagation using checkerboard update scheme
void bp_cb(image<float[VALUES]> *u, image<float[VALUES]> *d,
	   image<float[VALUES]> *l, image<float[VALUES]> *r,
	   image<float[VALUES]> *data,
	   int iter) {
  int width = data->width();  
  int height = data->height();

  for (int t = 0; t < ITER; t++) {
    std::cout << "iter " << t << "\n";

    for (int y = 1; y < height-1; y++) {
      for (int x = ((y+t) % 2) + 1; x < width-1; x+=2) {
	msg(imRef(u, x, y+1),imRef(l, x+1, y),imRef(r, x-1, y),
	    imRef(data, x, y), imRef(u, x, y));

	msg(imRef(d, x, y-1),imRef(l, x+1, y),imRef(r, x-1, y),
	    imRef(data, x, y), imRef(d, x, y));

	msg(imRef(u, x, y+1),imRef(d, x, y-1),imRef(r, x-1, y),
	    imRef(data, x, y), imRef(r, x, y));

	msg(imRef(u, x, y+1),imRef(d, x, y-1),imRef(l, x+1, y),
	    imRef(data, x, y), imRef(l, x, y));
      }
    }
  }
}

// multiscale belief propagation for image restoration
image<uchar> *restore_ms(image<uchar> *img) {


	//image<float[VALUES]> *data;
	//image<uchar> * tiles;

	//Data costs
	//data = comp_data(img);
	// BP
	//bp_cb(u, d, l, r, data, ITER);    

	int width = img->width();
	int height = img->height();
	// Number of tiles in each direction ignores the outer ring
	//Number of tiles in the image
	int numWidthTiles = (((width - 2) / TWIDTH) + 1);
	int numHeightTiles = (((height - 2) / THEIGHT) + 1);
	// Remaining number of rows/cols in last tile
	int remWidth = width-2 - ((numWidthTiles - 1)*TWIDTH);
	int remHeight = height-2 - ((numHeightTiles - 1)*THEIGHT);

	// The actual width and height of the tile
	image<int> *wTile, *hTile;

	// the tile image splits the image object into several image objects with width <= TWIDTH
	//image<image<uchar>> *tile = new image<image<uchar>>(numWidthTiles, numHeightTiles);

	wTile = new image<int>(numWidthTiles, 1);
	hTile = new image<int>(numHeightTiles, 1);

	int cumTileWidth = 0, cumTileHeight = 0;

	/*for (int q = 0; q < numHeightTiles; q++) {


		for (int p = 0; p < numWidthTiles; p++) {



			cumTileWidth += imRef(wTile, p, 0);

			// the rawTile image splits each tile image objects into image objects that receive raw pixel data.
			image<uchar> *rawTile = new image<uchar>(imRef(wTile, p, 0), imRef(hTile, q, 0));

			for (int r = 0; r < imRef(hTile, q, 0); r++) {
				
				// The image row/col corresponding to the tile pixel location
				int imgCol = cumTileHeight + r;
				for (int s = 0; s < imRef(wTile, p, 0); s++) {

					// The image row/col corresponding to the tile pixel location
					int imgRow = cumTileWidth + s;

					imRef(rawTile, s, r) = imRef(img, imgRow, imgCol);
				}
			}
			imRef(tile, p, q) = *rawTile;

			//delete rawTile;

		}

		cumTileHeight += imRef(hTile, q, 0);
		cumTileWidth = 0;
	}
	*/
	int x = 0;
	for (int q = 0; q < numHeightTiles; q++) {
		if (q == 0 || q == numHeightTiles - 1) imRef(hTile, q, 0) = 1;
		else if (q == numHeightTiles - 2) imRef(hTile, q, 0) = remHeight;
		else imRef(hTile, q, 0) = THEIGHT;

		std::cout << "Row: " << imRef(hTile, q, 0);
		std::cout << ("\n");
		x += imRef(hTile, q, 0);
	}
	x = 0;
	for (int p = 0; p < numWidthTiles; p++) {
		if (p == 0 || p == numWidthTiles - 1) imRef(wTile, p, 0) = 1;
		else if (p == numWidthTiles - 2) imRef(wTile, p, 0) = remWidth;
		else imRef(wTile, p, 0) = TWIDTH;

		std::cout << "Column: " << imRef(wTile, p, 0);
		std::cout << ("\n");
		x += imRef(wTile, p, 0);
	}
	std::cout << "Cumulative total columns: " << x << "\n";

	image<float[VALUES]> *data = new image<float[VALUES]>(width, height);
	image<float[VALUES]> *uBound;
	image<float[VALUES]> *dBound;
	image<float[VALUES]> *lBound;
	image<float[VALUES]> *rBound;
	image<float[VALUES]> *u;
	image<float[VALUES]> *d;
	image<float[VALUES]> *l;
	image<float[VALUES]> *r;
	image<uchar> *out = new image<uchar>(width, height);


	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			for (int value = 0; value < VALUES; value++) {
				float val = square((float)(imRef(img, x, y) - value));
				imRef(data, x, y)[value] = LAMBDA * std::min(val, DATA_K);
			}
		}
	}

	//std::vector<std::array<float,NULL> >;
	uBound = new image<float[VALUES]>(width, numHeightTiles, false);
	dBound = new image<float[VALUES]>(width, numHeightTiles, false);
	lBound = new image<float[VALUES]>(numWidthTiles, height, false);
	rBound = new image<float[VALUES]>(numWidthTiles, height, false);

	for (int t1 = 0; t1 < ITER; t1++) {
		std::cout << "iter " << t1 << "\n";

		cumTileHeight = 0;
		cumTileWidth = 0;

		//int c = 0;
		for (int col = 1; col < numHeightTiles-1; col++) {
			for (int row = 1; row < numWidthTiles-1; row++) {

				u = new image<float[VALUES]>(imRef(wTile, row, 0), imRef(hTile, col, 0), false);
				d = new image<float[VALUES]>(imRef(wTile, row, 0), imRef(hTile, col, 0), false);
				l = new image<float[VALUES]>(imRef(wTile, row, 0), imRef(hTile, col, 0), false);
				r = new image<float[VALUES]>(imRef(wTile, row, 0), imRef(hTile, col, 0), false);

				//std::cout << "U Width: " << imRef(wTile, row, 0) << "\n" << "U Height: " << imRef(hTile, col, 0) << "\n";
				//if (imRef(wTile, row, 0) < TWIDTH) {
				//	c += 1;
				//}
				//std::cout << "Number of Rows: " << c + 1 << "\n";
				//std::cout << "Row: " << col << "\n";
				//std::cout << "Column: " << row << "\n";

				cumTileWidth += imRef(wTile, row, 0);

				for (int t0 = 0; t0 < ITER0; t0++) {
					///std::cout << "iter0 " << t0 << "\n";
					for (int y = 0; y < imRef(hTile, col, 0); y++) {
						for (int x = ((y + t0) % 2); x < imRef(wTile, row, 0); x += 2) {

							if (y == 0 && x == 0) {
								//std::cout << "height: " << imRef(hTile, col, 0) << "\n";
								//std::cout << "width: " << imRef(wTile, row, 0) << "\n";
							}

							char node = '0';
							/*
								[ 7 ] [ 8 ] [ 9 ]
								[ 4 ] [ 5 ] [ 6 ]
								[ 1 ] [ 2 ] [ 3 ]
								Bottom Right Corner = '3'	Right Edge = '6'
								Bottom Left Corner = '1'	Left Edge = '4'
								Top Right Corner = '7'		Top Edge = '8'
								Top Left Corner = '9'		Bottom Edge = '2'
							*/
							if (x == 0 && y != 0 && y != (imRef(hTile, col, 0) - 1)) {
								node = '4';
							}
							else if ((x == imRef(wTile, col, 0) - 1) && y != 0 && y != (imRef(hTile, col, 0) - 1)) {
								node = '6';
							}
							else if (x != 0 && x != imRef(wTile, col, 0) - 1 && y == 0) {
								node = '8';
							}
							else if (x != 0 && x != imRef(wTile, col, 0) - 1 && y == imRef(hTile, col, 0) - 1) {
								node = '2';
							}
							else if (x == imRef(wTile, col, 0) - 1 && y == imRef(hTile, col, 0) - 1) {
								node = '3';
							}
							else if (x == 0 && y == imRef(hTile, col, 0) - 1) {
								node = '1';
							}
							else if (x == imRef(wTile, col, 0) - 1 && y == 0) {
								node = '9';
							}
							else if (x == 0 && y == 0) {
								node = '7';
							}

							if (t0 == 0) {
								//if (x == 0 && y != 0 && y != (imRef(hTile, col, 0) - 1)) {
									//std::cout << imRef(r, x - 1, y) << "\n" << imRef(rBound, row - 1, col);
									//imRef(r, x, y) = imRef(rBound, row - 1, col);
								//}
								switch (node) {
								case '1':
									msg1(imRef(d, x, y - 1), imRef(l, x + 1, y),
										imRef(data, cumTileWidth+x, cumTileHeight+y), imRef(lBound, row, cumTileHeight+y));

									msg1(imRef(d, x, y - 1), imRef(l, x + 1, y),
										imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(dBound, cumTileWidth + x, col));
									break;
								case '2':
									msg(imRef(uBound, cumTileWidth + x, col + 1), imRef(l, x + 1, y), imRef(r, x - 1, y),
										imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(u, x, y));

									msg(imRef(uBound, cumTileWidth + x, col + 1), imRef(d, x, y - 1), imRef(r, x - 1, y),
										imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(r, x, y));

									msg(imRef(uBound, cumTileWidth + x, col + 1), imRef(d, x, y - 1), imRef(l, x + 1, y),
										imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(l, x, y));
									break;
								case '3':
									msg(imRef(uBound, cumTileWidth + x, col + 1), imRef(lBound, x + 1, cumTileHeight+y), imRef(r, x - 1, y),
										imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(u, x, y));

									msg(imRef(uBound, cumTileWidth + x, col + 1), imRef(lBound, x + 1, cumTileHeight + y), imRef(d, x, y - 1),
										imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(l, x, y));
									break;
								case '4':
									msg(imRef(u, x, y + 1), imRef(l, x + 1, y), imRef(rBound, row - 1, cumTileHeight + y),
										imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(u, x, y));

									msg(imRef(d, x, y - 1), imRef(l, x + 1, y), imRef(rBound, row - 1, cumTileHeight + y),
										imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(d, x, y));

									msg(imRef(u, x, y + 1), imRef(d, x, y - 1), imRef(rBound, row - 1, cumTileHeight + y),
										imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(r, x, y));
									break;
								case '6':
									msg(imRef(u, x, y + 1), imRef(lBound, row + 1, cumTileHeight + y), imRef(r, x - 1, y),
										imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(u, x, y));

									msg(imRef(d, x, y - 1), imRef(lBound, row + 1, cumTileHeight + y), imRef(r, x - 1, y),
										imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(d, x, y));

									msg(imRef(u, x, y + 1), imRef(d, x, y - 1), imRef(lBound, row + 1, cumTileHeight + y),
										imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(l, x, y));
									break;
								case '7':
									msg(imRef(dBound, cumTileWidth + x, col - 1), imRef(l, x + 1, y), imRef(rBound, row - 1, cumTileHeight + y),
										imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(d, x, y));

									msg(imRef(u, x, y + 1), imRef(dBound, cumTileWidth + x, col - 1), imRef(rBound, row - 1, cumTileHeight + y),
										imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(r, x, y));
									break;
								case '8':
									msg(imRef(dBound, cumTileWidth + x, col - 1), imRef(l, x + 1, y), imRef(r, x - 1, y),
										imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(d, x, y));

									msg(imRef(u, x, y + 1), imRef(dBound, cumTileWidth + x, col - 1), imRef(r, x - 1, y),
										imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(r, x, y));

									msg(imRef(u, x, y + 1), imRef(dBound, cumTileWidth + x, col - 1), imRef(l, x + 1, y),
										imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(l, x, y));
									break;
								case '9':
									msg(imRef(dBound, cumTileWidth + x, col - 1), imRef(r, x - 1, y), imRef(lBound, x + 1, cumTileHeight + y),
										imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(d, x, y));

									msg(imRef(u, x, y + 1), imRef(dBound, cumTileWidth + x, col - 1), imRef(lBound, x + 1, cumTileHeight + y),
										imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(l, x, y));
									break;
								default:
									msg(imRef(u, x, y + 1), imRef(l, x + 1, y), imRef(r, x - 1, y),
										imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(u, x, y));

									msg(imRef(d, x, y - 1), imRef(l, x + 1, y), imRef(r, x - 1, y),
										imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(d, x, y));

									msg(imRef(u, x, y + 1), imRef(d, x, y - 1), imRef(r, x - 1, y),
										imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(r, x, y));

									msg(imRef(u, x, y + 1), imRef(d, x, y - 1), imRef(l, x + 1, y),
										imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(l, x, y));
									break;

								}
							}
							else if (t0 == ITER0 - 1) {
								switch (node) {
								case '1':
									msg1(imRef(d, x, y - 1), imRef(l, x + 1, y),
										imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(lBound, row, cumTileHeight + y));

									msg1(imRef(d, x, y - 1), imRef(l, x + 1, y),
										imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(dBound, cumTileWidth + x, col));
									break;
								case '2':
									msg(imRef(d, x, y - 1), imRef(l, x + 1, y), imRef(r, x - 1, y),
										imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(dBound, cumTileWidth + x, col));
									break;
								case '3':
									msg1(imRef(d, x, y - 1), imRef(r, x - 1, y),
										imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(rBound, row, cumTileHeight + y));

									msg1(imRef(d, x, y - 1), imRef(r, x - 1, y),
										imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(dBound, cumTileWidth + x, col));
									break;
								case '4':
									msg(imRef(u, x, y + 1), imRef(d, x, y - 1), imRef(l, x + 1, y),
										imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(lBound, row, cumTileHeight + y));
									break;
								case '6':
									msg(imRef(u, x, y + 1), imRef(d, x, y - 1), imRef(r, x - 1, y),
										imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(rBound, row, cumTileHeight + y));
									break;
								case '7':
									msg1(imRef(u, x, y + 1), imRef(l, x + 1, y),
										imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(lBound, row, cumTileHeight + y));

									msg1(imRef(u, x, y + 1), imRef(l, x + 1, y),
										imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(uBound, cumTileWidth + x, col));
									break;
								case '8':
									msg(imRef(u, x, y + 1), imRef(l, x + 1, y), imRef(r, x - 1, y),
										imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(uBound, cumTileWidth + x, col));
									break;
								case '9':
									msg1(imRef(u, x, y + 1), imRef(r, x - 1, y),
										imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(uBound, cumTileWidth + x, col));

									msg1(imRef(u, x, y + 1), imRef(r, x - 1, y),
										imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(rBound, row, cumTileHeight + y));
									break;
								default:
									msg(imRef(u, x, y + 1), imRef(l, x + 1, y), imRef(r, x - 1, y),
										imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(u, x, y));

									msg(imRef(d, x, y - 1), imRef(l, x + 1, y), imRef(r, x - 1, y),
										imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(d, x, y));

									msg(imRef(u, x, y + 1), imRef(d, x, y - 1), imRef(r, x - 1, y),
										imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(r, x, y));

									msg(imRef(u, x, y + 1), imRef(d, x, y - 1), imRef(l, x + 1, y),
										imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(l, x, y));
									break;

								}
							}
							else {
								switch (node) {
							case '1':
								msg2(imRef(l, x + 1, y), imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(u, row, col));

								msg2(imRef(d, x, y - 1), imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(r, row, col));
								break;
							case '2':
								msg1(imRef(l, x + 1, y), imRef(r, x - 1, y),
									imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(u, x, y));

								msg1(imRef(d, x, y - 1), imRef(r, x - 1, y),
									imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(r, x, y));

								msg1(imRef(d, x, y - 1), imRef(l, x + 1, y),
									imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(l, x, y));
								break;
							case '3':
								msg2(imRef(r, x - 1, y), imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(u, row, col));

								msg2(imRef(d, x, y - 1), imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(l, row, col));
								break;
							case '4':
								msg1(imRef(u, x, y + 1), imRef(l, x + 1, y),
									imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(u, x, y));

								msg1(imRef(d, x, y - 1), imRef(l, x + 1, y),
									imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(d, x, y));

								msg1(imRef(u, x, y + 1), imRef(d, x, y - 1),
									imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(r, x, y));
								break;
							case '6':
								msg1(imRef(u, x, y + 1), imRef(r, x - 1, y),
									imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(u, x, y));

								msg1(imRef(d, x, y - 1), imRef(r, x - 1, y),
									imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(d, x, y));

								msg1(imRef(u, x, y + 1), imRef(d, x, y - 1),
									imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(l, x, y));
								break;
							case '7':
								msg2(imRef(l, x + 1, y), imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(d, x, y));

								msg2(imRef(u, x, y + 1), imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(r, x, y));
								break;
							case '8':
								msg1(imRef(l, x + 1, y), imRef(r, x - 1, y),
									imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(d, x, y));

								msg1(imRef(u, x, y + 1), imRef(r, x - 1, y),
									imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(r, x, y));

								msg1(imRef(u, x, y + 1), imRef(l, x + 1, y),
									imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(l, x, y));
								break;
							case '9':
								msg2(imRef(r, x - 1, y), imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(d, x, y));

								msg2(imRef(u, x, y + 1), imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(l, x, y));
								break;
							default:
								msg(imRef(u, x, y + 1), imRef(l, x + 1, y), imRef(r, x - 1, y),
									imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(u, x, y));

								msg(imRef(d, x, y - 1), imRef(l, x + 1, y), imRef(r, x - 1, y),
									imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(d, x, y));

								msg(imRef(u, x, y + 1), imRef(d, x, y - 1), imRef(r, x - 1, y),
									imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(r, x, y));

								msg(imRef(u, x, y + 1), imRef(d, x, y - 1), imRef(l, x + 1, y),
									imRef(data, cumTileWidth + x, cumTileHeight + y), imRef(l, x, y));
								break;

							}
							}
							/*msg(imRef(u, x, y + 1), imRef(l, x + 1, y), imRef(r, x - 1, y),
								imRef(data, x, y), imRef(u, x, y));

							msg(imRef(d, x, y - 1), imRef(l, x + 1, y), imRef(r, x - 1, y),
								imRef(data, x, y), imRef(d, x, y));

							msg(imRef(u, x, y + 1), imRef(d, x, y - 1), imRef(r, x - 1, y),
								imRef(data, x, y), imRef(r, x, y));

							msg(imRef(u, x, y + 1), imRef(d, x, y - 1), imRef(l, x + 1, y),
								imRef(data, x, y), imRef(l, x, y));*/

								

						}
					}
				}
				for (int y = 1; y < imRef(hTile, col, 0)-1; y++) {
					for (int x = 1; x < imRef(wTile, row, 0)-1; x ++) {
						int best = 0;
						float best_val = INF;
						for (int value = 0; value < VALUES; value++) {
							float val =
								imRef(u, x, y + 1)[value] +
								imRef(d, x, y - 1)[value] +
								imRef(l, x + 1, y)[value] +
								imRef(r, x - 1, y)[value] +
								imRef(data, cumTileWidth + x, cumTileHeight + y)[value];
							if (val < best_val) {
								best_val = val;
								best = value;
							}
						}
						imRef(out, cumTileWidth+x, cumTileHeight+y) = best;
					}
				}

				delete u;
				delete d;
				delete l;
				delete r;
			}
			cumTileHeight += imRef(hTile, col, 0);
			cumTileWidth = 0;

		}
		
	}
	delete uBound;
	delete dBound;
	delete rBound;
	delete lBound;
	//u = new image<float[VALUES]>(2,2, false);
	//d = new image<float[VALUES]>(2, 2, false);
	//l = new image<float[VALUES]>(2, 2, false);
	//r = new image<float[VALUES]>(2, 2, false);

	//image<uchar> *out = output(u, d, l, r, data);
	//delete u;
	//delete d;
	//delete l;
	//delete r;
	delete data;
	//delete tile;

	return out;
}

int main(int argc, char **argv) {
  image<uchar> *img, *out, *edges;

 if (argc != 3) {
    std::cerr << "usage: " << argv[0] << " in(pgm) out(pgm)\n";
    exit(1);
  }
 

  // load input
  img = loadPGM(argv[1]);

  // restore
  out = restore_ms(img);

  // save output
  savePGM(out, argv[2]);
  
  delete img;
  delete out;
  return 0;
}
