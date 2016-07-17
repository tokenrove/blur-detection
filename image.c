/* Image loading abstraction.
 *
 * Presently only for JPEGs; added in case we need to also support TIFF or various RAW.
 *
 * Note that we could probably save some miniscule amount of memory
 * reusing the JPEG decompression object or similar.
 *
 * Copyright 2010 Julian Squires <julian@cipht.net>.
 */

#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "image.h"
#include <jpeglib.h>

struct image
{
	char *path;
	FILE *in;
	unsigned char *rgb, *luma;
	struct jpeg_decompress_struct cinfo;
        struct jpeg_error_mgr jerr;
	jvirt_barray_ptr *coeffp;
};

void *image_cinfo(void *image) { return &(((struct image *)image)->cinfo); }
void *image_coeffp(void *image) { return ((struct image *)image)->coeffp; }

static size_t bytes_per_scanline(struct image *image)
{
	return image->cinfo.output_width*image->cinfo.output_components;
}

unsigned int image_width(void *image)
{
	return ((struct image *)image)->cinfo.output_width;
}

unsigned int image_height(void *image)
{
	return ((struct image *)image)->cinfo.output_height;
}

unsigned int image_components(void *image)
{
	return ((struct image *)image)->cinfo.output_components;
}

unsigned char *image_rgb(void *image) { return ((struct image *)image)->rgb; }
unsigned char *image_luma(void *image) { return ((struct image *)image)->luma; }


static void populate_pixels(struct image *image, unsigned char **dst)
{
	unsigned char *p;

	p = *dst = malloc(image_height(image)*bytes_per_scanline(image));
	assert(p != NULL);

	for(int n; image->cinfo.output_scanline < image->cinfo.output_height; p += n*bytes_per_scanline(image))
		n = jpeg_read_scanlines(&image->cinfo, &p, 1);
}

/* Note: not for use in long-lived programs; error cases can leak memory. */
void *image_open(char *path, colorspace_t colorspace)
{
	struct image *image;

	image = malloc(sizeof(struct image));
	if(image == NULL) return NULL;
	memset(image, 0, sizeof(struct image));

        image->cinfo.err = jpeg_std_error(&image->jerr);
        jpeg_create_decompress(&image->cinfo);
	if((image->in = fopen(path, "rb")) == NULL) {
		fprintf(stderr, "%s: Couldn't open.\n", path);
		jpeg_destroy_decompress(&image->cinfo);
		free(image);
		return NULL;
	}
	image->path = path;
	jpeg_stdio_src(&image->cinfo, image->in);
	jpeg_read_header(&image->cinfo, TRUE);
	if(colorspace == LUMA_COLORSPACE)
		image->cinfo.out_color_space = JCS_GRAYSCALE;
	jpeg_start_decompress(&image->cinfo);

	if((colorspace == RGB_COLORSPACE && image->cinfo.output_components != 3) ||
	   (colorspace == LUMA_COLORSPACE && image->cinfo.output_components != 1)) {
		fprintf(stderr, "%s: Non-RGB/YCC files are unsupported (got %d components instead of %d).\n",
			path, image->cinfo.output_components, (colorspace == RGB_COLORSPACE) ? 3 : 1);
		jpeg_destroy_decompress(&image->cinfo);
		fclose(image->in);
		free(image);
		return NULL;
	}

	image->rgb = NULL;
	image->luma = NULL;
	populate_pixels(image, (colorspace == LUMA_COLORSPACE) ? &image->luma : &image->rgb);
	jpeg_finish_decompress(&image->cinfo);

	/* XXX a bit of a hack; a better solution would be using
	 * libjpeg's buffered image mode */
	fclose(image->in);
	if((image->in = fopen(path, "rb")) == NULL)
		return NULL;
	jpeg_stdio_src(&image->cinfo, image->in);
	jpeg_read_header(&image->cinfo, TRUE);
	image->coeffp = jpeg_read_coefficients(&image->cinfo);

	return image;
}

void image_close(void *image_)
{
	struct image *image = image_;

	if(image->rgb != NULL) free(image->rgb);
	if(image->luma != NULL) free(image->luma);
	jpeg_destroy_decompress(&image->cinfo);
	fclose(image->in);
	memset(image, 0, sizeof(struct image));
	free(image);
}

#define M_PI_F 3.1415926535897932384626433832795f

inline static float minf(float a, float b, float c)
{
	return (a < b) ? ((a < c) ? a : c) : ((b < c) ? b : c);
}

int image_hsi_of_block(void *image_, int x, int y, float *hue, float *saturation, float *intensity)
{
	struct image *image = image_;
	int i, j, n, w, h;
	unsigned char *rgb;

	assert(3 == image_components(image_));

	w = BLOCK_WIDTH; h = BLOCK_HEIGHT;
	assert(x+w <= image_width(image_));
	assert(y+h <= image_height(image_));

	rgb = image->rgb;
	rgb += image_components(image)*(x + (y * image_width(image)));
	for(i = 0, n = 0; i < h; i++, rgb += 3*(image_width(image_)-w)) {
		for(j = 0; j < w; j++, n++, rgb += 3, hue++, saturation++, intensity++) {
			float r, g, b, sum;

			sum = rgb[0] + rgb[1] + rgb[2] + EPSILON_F;
			r = (float)rgb[0] / sum;
			g = (float)rgb[1] / sum;
			b = (float)rgb[2] / sum;
			*hue = acosf((0.5 * ((r-g) + (r-b))) /
				     (EPSILON_F + sqrtf(((r-g)*(r-g)) + ((r-b) * (g-b)))));
			if(b > g) *hue = (2.0*M_PI_F) - *hue;
			*saturation = 1.0 - 3.0 * minf(r, g, b);
			*intensity = (float)(rgb[0]+rgb[1]+rgb[2]) / 765.0;
		}
	}

	return n;
}
