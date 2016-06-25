/* Detect out-of-focus images
 *
 * Per the algorithm given by Lim, Yen, and Wu in "Detection of
 * Out-Of-Focus Digital Photographs".
 *
 * Potential optimizations:
 *   - fitting operations to cache size
 *   - performing filters on GPU (OpenCL)
 *   - online/incremental mean and variance computation
 *   - change image.c to only keep BLOCK_HEIGHT scanlines at a time (reduce resident set)
 *   - see also notes in IIR filter section
 *
 * Many things can be tuned to change the behavior of this algorithm.  For example:
 *   - smoothing of the rule of thirds weighting matrix
 *   - sharpness threshold
 *   - IIR filter coefficients
 *   - highpass filter could be replaced with Canny-Deriche
 *   - sky hue values
 *   - block size (defined in image.h presently) -- may affect the extent to which
 *     large features are incorporated in local sharpness metrics
 *
 * Copyright 2010 Julian Squires <julian@cipht.net>.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "image.h"
#include "tong-li-zhang-zhang.h"

/**** TUNABLE VALUES ****
 * See also BLOCK_WIDTH and BLOCK_HEIGHT in image.h.
 */

static enum { IIR_FILTER, TONG_LI_ZHANG_ZHANG, MARICHAL_MA_ZHANG } strategy = IIR_FILTER;

static float sharp_threshold = 0.32; /* -s= threshold at which a block is considered sharp */
#define SKY_HUE_LOW 3.1	    /* hue, in radians, to detect (ignored) sky blocks */
#define SKY_HUE_HIGH 3.8
static float sky_threshold = 0.7; /* -k= threshold of sky hue to ignore */

/** IIR filter strategy tuning values */
static int sharp_stride = 1; /* -r= increase to skip rows/columns for performance */
/* -b= how much bandpass energy to consider part of feature range when
 * calculating sharpness */
static float min_bandpass_level = 0.12;
static float sharp_scale = 15.0;	/* -x= scaling factor on sharpness values */


/**** DATA STRUCTURES ****/

struct block_metrics
{
	float sharpness;
	float reflectivity;	/* variance of sharpness samples */
	float brightness;
	float mean_saturation;
	float mean_hue;
	float sky_ratio;
};

struct figures_of_merit
{
	float composition;	/* sharpness weighted by rule of thirds */
	float brightness_idx;
	float saturation_idx;
	float density;
	float median_sharpness;
};

struct image
{
	char *path;
	void *handle;
	int width_in_blocks, height_in_blocks;
	struct block_metrics *blocks;
	struct figures_of_merit merit;
	int decision;
} image;

float *hue, *saturation, *intensity;
float *column_sharpness, *row_sharpness;
float **emap;			/* for TLZZ metric */

/**** GLOBAL MERITS ****/

struct {
	int w, h;
	float *weight;
} thirds;

static float manhattan_distance(float px, float py, float qx, float qy, int w, int h)
{
	return (fabsf(px-qx)/w)+(fabsf(py-qy)/h);
}

static float min_manhattan_distance(float x, float y, float *p, int n, int w, int h)
{
	float min = INFINITY, j;
	for(int i = 0; i < n; i++) {
		j = manhattan_distance(x, y, p[2*i], p[2*i+1], w, h);
		if(j < min) min = j;
	}
	return min;
}

static void rebuild_thirds_matrix(int w, int h)
{
	float power_points[4*2] = {
		1.0/3,1.0/3,  2.0/3,1.0/3,
		1.0/3,2.0/3,  2.0/3,2.0/3 };

	if(thirds.weight != NULL) free(thirds.weight);

	thirds.w = w;
	thirds.h = h;
	thirds.weight = malloc(sizeof(float)*w*h);

	for(int i = 0; i < 4; i++) {
		power_points[2*i] *= w;
		power_points[2*i+1] *= h;
		power_points[2*i] += 0.5;
		power_points[2*i+1] += 0.5;
	}

	for(int i = 0; i < h; i++)
		for(int j = 0; j < w; j++)
			thirds.weight[j+i*w] = 1.0 - min_manhattan_distance((float)j+0.5, (float)i+0.5, power_points, 4, w, h);
}


static int float_cmp(const void *a_, const void *b_)
{
	float a = *(float *)a_, b = *(float *)b_;
	return (a < b) ? -1 : (a > b) ? 1 : 0;
}

static void compute_global_merits(struct image *image)
{
	float composition_sum, mean_bright_sharp, mean_bright_blur, mean_sat_sharp, mean_sat_blur;
	float *sorted_sharpness;
	int n, sharp_count;

	sorted_sharpness = malloc(sizeof(float)*image->width_in_blocks*image->height_in_blocks);
	composition_sum = mean_bright_sharp = mean_bright_blur = mean_sat_sharp = mean_sat_blur = 0.0;
	n = sharp_count = 0;
	for(int i = 0; i < image->width_in_blocks*image->height_in_blocks; i++) {
		if(image->blocks[i].sky_ratio > sky_threshold) continue;
		n++;
		sorted_sharpness[i] = image->blocks[i].sharpness;
		composition_sum += image->blocks[i].sharpness * thirds.weight[i];
		if(image->blocks[i].sharpness > sharp_threshold) {
			sharp_count++;
			mean_bright_sharp += image->blocks[i].brightness;
			mean_sat_sharp += image->blocks[i].mean_saturation;
		} else {
			mean_bright_blur += image->blocks[i].brightness;
			mean_sat_blur += image->blocks[i].mean_saturation;
		}
	}

	qsort(sorted_sharpness, n, sizeof(float), float_cmp);
	if(n%2 == 0)
		image->merit.median_sharpness = (sorted_sharpness[n/2]+sorted_sharpness[n/2-1])/2.0;
	else
		image->merit.median_sharpness = sorted_sharpness[n/2];
	free(sorted_sharpness);

	mean_bright_sharp /= n;
	mean_bright_blur /= n;
	mean_sat_sharp /= n;
	mean_bright_blur /= n;
	image->merit.composition = composition_sum/n;
	image->merit.brightness_idx = mean_bright_sharp-mean_bright_blur;
	image->merit.saturation_idx = mean_sat_sharp-mean_sat_blur;
	image->merit.density = (float)sharp_count/n;
}


/**** BLOCK METRICS ****/

/* IIR filters
 * Coefficients per Shaked and Tastl, 2004.
 *
 * Notes and Caveats:
 *
 * - on some architectures it might be faster to compute position in the ring
 *    buffer without branches, so be sure to compare and profile.
 *
 * - this implementation depends on having a few padding zeros at the
 *   beginning of whatever input is fed to it.  If you repurpose this
 *   code, be sure to carefully look at how FILTER_PAD is used!
 *
 * - per http://www-users.cs.york.ac.uk/~fisher/mkfilter/trad.html, premultiply
 *   by gain to remove three multiplications from each filter; right now we use
 *   the alpha/beta arrays for flexibility (easy tweaking of values) but observe
 *   their symmetry: they can be embedded in the code, replaced by
 *   adds/subtracts.
 */

enum { FILTER_PAD = 4 };

/* butterworth filter, matlab butter(3, 0.1) */
static float bp_alpha[] = { 2.3741, -1.9294, 0.5321 };
static float bp_beta[] = { 0.0029, 0.0087, 0.0087, 0.0029 };
static float bp_ring[4];

static float bandpass(float *m, int x)
{
	float alpha, beta;
	beta = m[x]*bp_beta[0] + m[x-1]*bp_beta[1] + m[x-2]*bp_beta[2] + m[x-3]*bp_beta[3];
	alpha = bp_ring[(x+3)&3]*bp_alpha[0] + bp_ring[(x+2)&3]*bp_alpha[1] + bp_ring[(x+1)&3]*bp_alpha[2];
	return bp_ring[x&3] = alpha + beta;
}

/* butterworth highpass, 0.75 corner, originally calculated in Octave with
 *     butter(3, 0.75, 'high');
 * Recalculated at higher precision by
 *   http://www-users.cs.york.ac.uk/~fisher/mkfilter/trad.html
 * (same parameters) */
static float hp_alpha[] = { -1.4590290622, -0.9103690003, -0.1978251873 };
#define HP_GAIN 3.155634919e-02
static float hp_beta[] = { HP_GAIN, -3*HP_GAIN, 3*HP_GAIN, -HP_GAIN };
static float hp_ring[4];

static float highpass(float *m, int x)
{
	float alpha, beta;
	beta = m[x]*hp_beta[0] + m[x-1]*hp_beta[1] + m[x-2]*hp_beta[2] + m[x-3]*hp_beta[3];
	alpha = hp_ring[(x+3)&3]*hp_alpha[0] + hp_ring[(x+2)&3]*hp_alpha[1] + hp_ring[(x+1)&3]*hp_alpha[2];
	return hp_ring[x&3] = alpha + beta;
}

static float sharpness(float *m, int n)
{
	float sum, bp, hp, r;
	int lastx = -1;
	float lastr = 0.0;

	memset(bp_ring, 0, sizeof(bp_ring));
	memset(hp_ring, 0, sizeof(hp_ring));
	sum = 0.0;
	lastx = -1;
	for(int x = 0; x < n; x++) {
		bp = bandpass(m, x);
		hp = highpass(m, x); /* to keep the ring up to date  */
		if(bp < min_bandpass_level) continue;
		r = hp/bp;
		r *= r;
		if(lastx != -1)
			sum += (x-lastx) * (r + lastr);
		lastx = x; lastr = r;
	}
	sum *= 0.5;
	return sharp_scale*sum;
}

static void copy_column(float *d, float *s)
{
	for(int i = 0; i < BLOCK_HEIGHT; i++, d++, s += BLOCK_WIDTH) *d = *s;
}

static void copy_row(float *d, float *s)
{
	memcpy(d, s, BLOCK_WIDTH*sizeof(float));
}


/**** LIM-YEN-WU ****/

static void calculate_block_reflectivity(struct block_metrics *block, int n)
{
	int m = 0;
	float sum_of_squares = 0.0;

	for(int i = 0; i < BLOCK_WIDTH; i += sharp_stride, m++) {
		sum_of_squares += powf(column_sharpness[i]-block->sharpness, 2.0);
	}
	for(int i = 0; i < BLOCK_HEIGHT; i += sharp_stride, m++) {
		sum_of_squares += powf(row_sharpness[i]-block->sharpness, 2.0);
	}
	assert(m == n);
	block->reflectivity = sum_of_squares / (n-1);
}

static void calculate_block_sharpness_values(struct block_metrics *block)
{
	int n = 0;
	float pixels[((BLOCK_WIDTH > BLOCK_HEIGHT) ? BLOCK_WIDTH : BLOCK_HEIGHT)+FILTER_PAD] = {0};

	for(int i = 0; i < BLOCK_WIDTH; i += sharp_stride, n++) {
		copy_column(pixels+FILTER_PAD, intensity+i);
		column_sharpness[i] = sharpness(pixels+FILTER_PAD, BLOCK_HEIGHT);
		block->sharpness += column_sharpness[i];
	}

	for(int i = 0; i < BLOCK_HEIGHT; i += sharp_stride, n++) {
		copy_row(pixels+FILTER_PAD, intensity+i*BLOCK_WIDTH);
		row_sharpness[i] = sharpness(pixels+FILTER_PAD, BLOCK_WIDTH);
		block->sharpness += row_sharpness[i];
	}

	block->sharpness /= n;
	calculate_block_reflectivity(block, n);
}

static void calculate_block_means(int n, struct block_metrics *block)
{
	float sky_histogram = 0;
	for(int i = 0; i < n; i++) {
		block->mean_hue += hue[i];
		if(hue[i] > SKY_HUE_LOW && hue[i] < SKY_HUE_HIGH) sky_histogram++;
		block->mean_saturation += saturation[i];
		block->brightness += intensity[i];
	}
	block->sky_ratio = sky_histogram/n;
	block->mean_hue /= n;
	block->mean_saturation /= n;
	block->brightness /= n;
}

/**** TONG-LI-ZHANG-ZHANG ****/

/*
 * THRESHOLD could be defined based on mean intensity of
 * Laplacian-xform of data set, as in [2], but currently we use a
 * value close to that given in [1], which is based on the expected
 * limitations of human vision.
 *
 * Note that changing this value radically changes the blur extent
 * result.
 */
static float tlzz_threshold = 0.011764; /* -t= */

/* [1] recommends using a MIN_ZERO of 0.05, but since we don't do
 * proper NMS yet (see suppress_nonmaxima()), we set this a bit lower
 * in the expectation that we won't have pruned as many edge points.
 */
static float min_zero = 0.01;	/* -z= */

static void calculate_TLZZ_sharpness_metric(struct block_metrics *block)
{
	float da_ratio, blur_extent;

	haar_transform(intensity, BLOCK_WIDTH, BLOCK_HEIGHT);
	construct_edge_map(intensity, emap, BLOCK_WIDTH, BLOCK_HEIGHT);
	// XXX add arguments
	detect_blur(emap, BLOCK_WIDTH, BLOCK_HEIGHT, tlzz_threshold, &da_ratio, &blur_extent);
	block->sharpness = 1.0-blur_extent;
}

/**** MARICHAL-MA-ZHANG ****/

#include <jpeglib.h>

static int min_dct_value = 1;	/* -d= */
static float max_histogram_value = 0.005; /* -h= */

static float weights[] = {	/* diagonal weighting */
	8,7,6,5,4,3,2,1,
	1,8,7,6,5,4,3,2,
	2,1,8,7,6,5,4,3,
	3,2,1,8,7,6,5,4,
	4,3,2,1,8,7,6,5,
	5,4,3,2,1,8,7,6,
	6,5,4,3,2,1,8,7,
	7,6,5,4,3,2,1,8
};
static float total_weight = 344;

static inline void update_histogram(JCOEF *block, int *histogram)
{
	for(int k = 0; k < DCTSIZE2; k++, block++)
		if(abs(*block) > min_dct_value) histogram[k]++;
}

static float compute_blur(int *histogram)
{
	float blur = 0.0;
	for(int k = 0; k < DCTSIZE2; k++)
		if(histogram[k] < max_histogram_value*histogram[0])
			blur += weights[k];
	blur /= total_weight;
	return blur;
}

void calculate_MMZ_sharpness_metric(struct image *image, struct block_metrics *block, int y, int x)
{
	int histogram[DCTSIZE2] = {0};
	JBLOCKARRAY cs;
	struct jpeg_decompress_struct *cinfo;
	jvirt_barray_ptr *coeffp;

	cinfo = image_cinfo(image->handle);
	coeffp = image_coeffp(image->handle);

	/* Note: only looking at the luma; assuming it's the first component. */
	for(int i = 0; i < BLOCK_HEIGHT/DCTSIZE; i++) {
		cs = cinfo->mem->access_virt_barray((j_common_ptr)cinfo, coeffp[0], i+y*BLOCK_HEIGHT/DCTSIZE, 1, FALSE);
		for(int j = 0; j < BLOCK_WIDTH/DCTSIZE; j++)
			update_histogram(cs[0][j+x*BLOCK_WIDTH/DCTSIZE], histogram);
	}

	block->sharpness = 1.0 - compute_blur(histogram);
}

static void calculate_blockwise_merits(struct image *image, int i, int j)
{
	struct block_metrics *block;
	int n;

	block = &image->blocks[j+i*image->width_in_blocks];
	memset(block, 0, sizeof(struct block_metrics));
	n = image_hsi_of_block(image->handle, j*BLOCK_WIDTH, i*BLOCK_HEIGHT, hue, saturation, intensity);
	if(n != BLOCK_WIDTH*BLOCK_HEIGHT) {
		fprintf(stderr, "calculate_blockwise_merits: panic!  Not setup to deal with a non-block element (%d).\n", n);
		exit(1);
	}

	calculate_block_means(n, block);
	switch(strategy) {
	case TONG_LI_ZHANG_ZHANG:
		calculate_TLZZ_sharpness_metric(block);
		break;
	case MARICHAL_MA_ZHANG:
		calculate_MMZ_sharpness_metric(image, block, i, j);
		break;
	case IIR_FILTER:
	default:
		calculate_block_sharpness_values(block);
	}
}


/**** MAIN OPERATION ****/

static void output_metadata(struct image *image)
{
	// XXX embed as IPTC etc
	printf("%s -- composition: %f  density: %f  median: %f  brightdx: %f  satidx: %f\n",
	       image->decision ? "well-focused":"ill-focused",
	       image->merit.composition, image->merit.density, image->merit.median_sharpness,
	       image->merit.brightness_idx, image->merit.saturation_idx);
}

static int operate_on_image(char *path)
{
	struct image image;

	image.path = path;
	image.handle = image_open(path, RGB_COLORSPACE);
	if(image.handle == NULL) return 1;

	/* Note: we throw away non-block size pieces rather than adjusting to
	 * them.  With an appropriate block size, the effect shouldn't be too
	 * drastic, and allows more optimization of the analysis functions. */
	image.width_in_blocks = (int)floor((double)image_width(image.handle)/BLOCK_WIDTH);
	image.height_in_blocks = (int)floor((double)image_height(image.handle)/BLOCK_HEIGHT);
	image.blocks = malloc(sizeof(struct block_metrics)*image.width_in_blocks*image.height_in_blocks);
	if(image.blocks == NULL) return 1;
	memset(image.blocks, 0, sizeof(struct block_metrics)*image.width_in_blocks*image.height_in_blocks);

	for(int i = 0; i < image.height_in_blocks; i++)
		for(int j = 0; j < image.width_in_blocks; j++)
			calculate_blockwise_merits(&image, i, j);

	if(thirds.w != image.width_in_blocks || thirds.h != image.height_in_blocks)
		rebuild_thirds_matrix(image.width_in_blocks, image.height_in_blocks);

	compute_global_merits(&image);
	image.decision = 0;
	if(image.merit.density > 0.6 ||
	   (image.merit.composition > -0.1 &&
	    image.merit.median_sharpness > 0.03 &&
	    image.merit.brightness_idx > -0.15 &&
	    image.merit.saturation_idx > -0.15 &&
	    image.merit.density > 0.1))
		image.decision = 1;

	output_metadata(&image);
	image_close(image.handle);
	return 0;
}

int main(int argc, char **argv)
{
	int status, i;

	thirds.w = thirds.h = 0; thirds.weight = NULL;

	hue = malloc(sizeof(float) * BLOCK_WIDTH * BLOCK_HEIGHT);
	saturation = malloc(sizeof(float) * BLOCK_WIDTH * BLOCK_HEIGHT);
	intensity = malloc(sizeof(float) * BLOCK_WIDTH * BLOCK_HEIGHT);
	column_sharpness = malloc(sizeof(float) * BLOCK_WIDTH);
	row_sharpness = malloc(sizeof(float) * BLOCK_HEIGHT);
	if(!(hue && saturation && intensity && column_sharpness && row_sharpness)) return 1;
	memset(column_sharpness, 0, sizeof(float) * BLOCK_WIDTH);
	memset(row_sharpness, 0, sizeof(float) * BLOCK_HEIGHT);

	emap = malloc(sizeof(float*)*(LEVELS+1));
	memset(emap, 0, sizeof(float*)*(LEVELS+1));
	for(int i = 1; i <= LEVELS; i++) emap[i] = malloc(sizeof(float)*((BLOCK_WIDTH>>i)*(BLOCK_HEIGHT>>i)));

	for(i = 1, status = 0; i < argc; i++) {
		if(argv[i][0] == '-') {
			if(argv[i][1] == 's')
				sscanf(argv[i], "-s=%f", &sharp_threshold);
			else if(argv[i][1] == 'k')
				sscanf(argv[i], "-k=%f", &sky_threshold);
			else if(argv[i][1] == 'r')
				sscanf(argv[i], "-r=%i", &sharp_stride);
			else if(argv[i][1] == 'b')
				sscanf(argv[i], "-b=%f", &min_bandpass_level);
			else if(argv[i][1] == 'x')
				sscanf(argv[i], "-x=%f", &sharp_scale);
			else if(!strcmp(argv[i], "-a=filter"))
				strategy = IIR_FILTER;
			else if(!strcmp(argv[i], "-a=tlzz"))
				strategy = TONG_LI_ZHANG_ZHANG;
			else if(argv[i][1] == 't')
				sscanf(argv[i], "-t=%f", &tlzz_threshold);
			else if(argv[i][1] == 'z')
				sscanf(argv[i], "-z=%f", &min_zero);
			else if(!strcmp(argv[i], "-a=mmz"))
				strategy = MARICHAL_MA_ZHANG;
			else if(argv[i][1] == 'd')
				sscanf(argv[i], "-d=%d", &min_dct_value);
			else if(argv[i][1] == 'h')
				sscanf(argv[i], "-h=%f", &max_histogram_value);
			continue;
		}

		status |= operate_on_image(argv[i]);
	}

	for(int i = 1; i <= LEVELS; i++) free(emap[i]);
	free(emap);

	free(hue); free(saturation); free(intensity); free(column_sharpness); free(row_sharpness);
	if(thirds.weight) free(thirds.weight);
	return status;
}
