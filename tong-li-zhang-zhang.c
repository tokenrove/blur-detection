/* Wavelet-based Blur Detection
 *
 * Based on the algorithm given in [1].
 *
 * [1] Hanghang Tong, Mingjing Li, Hongjiang Zhang, Changshui Zhang,
 * "Blur detection for digital images using wavelet transform", 2004.
 *
 * [2] Harish Ramakrishnan, "Detection and Estimation of Image Blur",
 * 2010.
 *
 * Copyright 2010 Julian Squires <julian@cipht.net>
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include "image.h"
#include "tong-li-zhang-zhang.h"

/**** BLUR DETECTION ****/

void detect_blur(float **emap, int w, int h, float threshold, float *da_ratio, float *blur_extent)
{
	int n_edge, n_da, n_rg, n_brg;

	assert(LEVELS == 3);
	n_edge = n_da = n_rg = n_brg = 0;
	for(int l = 0; l < h; l += 2)
		for(int k = 0; k < w; k += 2) {
			float p[LEVELS];
			for(int i = 1; i <= LEVELS; i++) {
				int x,y;
				x = k>>i; if(x >= (w>>i)) x = (w>>i)-1;
				y = l>>i; if(y >= (h>>i)) y = (h>>i)-1;
				p[i-1] = emap[i][x+y*(w>>i)];
			}

			if(p[0] > threshold || p[1] > threshold || p[2] > threshold) {
				n_edge++;
				if(p[0] > p[1] && p[1] > p[2])
					n_da++;
				else if((p[0] < p[1] && p[1] < p[2]) || /* roof or gradual step */
					(p[1] > p[0] && p[1] > p[2])) { /* roof only */
					n_rg++;
					if(p[0] < threshold) n_brg++;
				}
			}
		}

	*da_ratio = (float)n_da/(n_edge+EPSILON_F);
	*blur_extent = (float)(n_brg+EPSILON_F)/(n_rg+EPSILON_F);
}


/**** EDGE DETECTION ****/

/* Per nomenclature of [1] and [2], (k,l) are coordinates on the
 * image, i is the decomposition level. */
#define LH dst[(l+(h>>i))*w + (k>>i)]
#define HL dst[l*w          + (k>>i)+(w>>i)]
#define HH dst[(l+(h>>i))*w + k+(w>>i)]
#define LL dst[l*w          + k]

void construct_edge_map(float *dst, float **emap, int w, int h)
{
	for(int i = 1; i <= LEVELS; i++)
		for(int l = 0; l < h>>i; l++)
			for(int k = 0; k < w>>i; k++)
				emap[i][k+l*(w>>i)] = sqrtf(LH*LH+HL*HL+HH*HH);
}



/**** HAAR TRANSFORM ****/

/* Note that if an image's dimension is odd, we discard the odd pixel
 * rather than dealing with padding.  This would cause degenerate
 * behavior on single pixel wide/high images, but I don't think we
 * care. */

// XXX we should be able to do this in-place
static void inner_haar_transform(float *dst, int w, int h, int stride)
{
	float *tmp;

	tmp = malloc(sizeof(float)*stride*h);
	memcpy(tmp, dst, sizeof(float)*stride*h);

	// for each row,
	// 1D haar transform
	for(int i = 0; i < h; i++) {
		for(int j = 0; j < (w&~1); j+=2) {
			dst[i*stride+j/2+w/2] = 0.5*(tmp[i*stride+j+1]-tmp[i*stride+j]);
			dst[i*stride+j/2] = 0.5*(tmp[i*stride+j]+tmp[i*stride+j+1]);
		}
	}
	memcpy(tmp, dst, sizeof(float)*stride*h);
        // column-wise
	for(int i = 0; i < w; i++)
		for(int j = 0; j < (h&~1); j+=2) {
			dst[(h/2+j/2)*stride+i] = 0.5*(tmp[stride*(1+j)+i]-tmp[stride*j+i]);
			dst[(j/2)*stride+i] = 0.5*(tmp[stride*(1+j)+i]+tmp[j*stride+i]);
		}

	free(tmp);
}

void haar_transform(float *data, int w, int h)
{
	for(int i = 0; i < LEVELS; i++)
		inner_haar_transform(data, w>>i, h>>i, w);
}
