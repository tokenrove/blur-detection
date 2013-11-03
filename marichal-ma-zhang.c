/* Fast blur detection using JPEG DCT coefficients
 *
 * Based on "Blur Determination in the Compressed Domain Using DCT
 * Information" by Xavier Marichal, Wei-Ying Ma, and Hong-Jiang Zhang.
 *
 * Tweak MIN_DCT_VALUE and MAX_HISTOGRAM_VALUE to adjust
 * effectiveness.  I reduced these values from those given in the
 * paper because I find the original to be less effective on large
 * JPEGs.
 *
 * Copyright 2010 Julian Squires <julian@cipht.net>
 */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
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


static int operate_on_image(char *path)
{
        struct jpeg_error_mgr jerr;
	struct jpeg_decompress_struct cinfo;
	jvirt_barray_ptr *coeffp;
	JBLOCKARRAY cs;
	FILE *in;
	int histogram[DCTSIZE2] = {0};

        cinfo.err = jpeg_std_error(&jerr);
        jpeg_create_decompress(&cinfo);
	if((in = fopen(path, "rb")) == NULL) {
		fprintf(stderr, "%s: Couldn't open.\n", path);
		jpeg_destroy_decompress(&cinfo);
		return 0;
	}
	jpeg_stdio_src(&cinfo, in);
	jpeg_read_header(&cinfo, TRUE);
	// XXX might be a little faster if we ask for grayscale
	coeffp = jpeg_read_coefficients(&cinfo);

	/* Note: only looking at the luma; assuming it's the first component. */
	for(int i = 0; i < cinfo.comp_info[0].height_in_blocks; i++) {
		cs = cinfo.mem->access_virt_barray((j_common_ptr)&cinfo, coeffp[0], i, 1, FALSE);
		for(int j = 0; j < cinfo.comp_info[0].width_in_blocks; j++)
			update_histogram(cs[0][j], histogram);
	}

	printf("%f\n", compute_blur(histogram));
	// output metadata XXX should be in IPTC etc

	// XXX also need to destroy coeffp?
	jpeg_destroy_decompress(&cinfo);
	return 0;
}

int main(int argc, char **argv)
{
	int status, i;

	for(status = 0, i = 1; i < argc; i++) {
		if(argv[i][0] == '-') {
			if(argv[i][1] == 'd')
				sscanf(argv[i], "-d=%d", &min_dct_value);
			else if(argv[i][1] == 'h')
				sscanf(argv[i], "-h=%f", &max_histogram_value);
			continue;
		}
		status |= operate_on_image(argv[i]);
	}

	return status;
}
