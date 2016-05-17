#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include "image.h"
#include "tong-li-zhang-zhang.h"

/**** TUNABLE VALUES ****/

/*
 * THRESHOLD could be defined based on mean intensity of
 * Laplacian-xform of data set, as in [2], but currently we use a
 * value close to that given in [1], which is based on the expected
 * limitations of human vision.
 *
 * Note that changing this value radically changes the blur extent
 * result.
 */
static float threshold = 0.011764; /* -t= */

/* [1] recommends using a MIN_ZERO of 0.05, but since we don't do
 * proper NMS yet (see suppress_nonmaxima()), we set this a bit lower
 * in the expectation that we won't have pruned as many edge points.
 */
static float min_zero = 0.01;	/* -z= */

/**** MAIN OPERATION ****/

static int operate_on_image(char *path)
{
	void *image;
	unsigned char *luma;
	float *out, **emap, da_ratio, blur_extent;
	int n, w, h;

	image = image_open(path, LUMA_COLORSPACE);
	if(image == NULL) return 1;

	luma = image_luma(image);
	w = image_width(image);
	h = image_height(image);
	n = w*h;
	out = malloc(sizeof(float)*n);
	assert(out != NULL);
	for(int i = 0; i < n; i++) out[i] = luma[i]/255.0;
	haar_transform(out, w, h);

	emap = malloc(sizeof(float*)*(LEVELS+1));
	memset(emap, 0, sizeof(float*)*(LEVELS+1));
	for(int i = 1; i <= LEVELS; i++) emap[i] = malloc(sizeof(float)*((w>>i)*(h>>i)));
	construct_edge_map(out, emap, w, h);
	/* Note: we don't perform non-maxima suppression at this
	 * point, because it seems to yield worse results compared to
	 * simply adjusting the threshold. */

	detect_blur(emap, w, h, threshold, &da_ratio, &blur_extent);
	printf("%s -- da_ratio: %f  blur_extent: %f\n", (da_ratio > min_zero) ? "sharp":"blurred",
	       da_ratio, blur_extent);

	for(int i = 1; i <= LEVELS; i++) free(emap[i]);
	free(emap);
	free(out);
	return 0;
}


int main(int argc, char **argv)
{
	int status, i;

	for(status = 0, i = 1; i < argc; i++) {
		if(argv[i][0] == '-') {
			if(argv[i][1] == 't')
				sscanf(argv[i], "-t=%f", &threshold);
			else if(argv[i][1] == 'z')
				sscanf(argv[i], "-z=%f", &min_zero);
			continue;
		}

		status |= operate_on_image(argv[i]);
	}

	return status;
}
