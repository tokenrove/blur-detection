
#ifndef IMAGE_H
#define IMAGE_H

/* Tunable throughout, but I don't see why we wouldn't go with
 * something that's an even number of DCT blocks. */
enum { BLOCK_WIDTH=64, BLOCK_HEIGHT=64 }; /* in pixels */
#define EPSILON_F 5.960465e-8
typedef enum { RGB_COLORSPACE, LUMA_COLORSPACE } colorspace_t;

extern void *image_open(char *path, colorspace_t colorspace);
extern void image_close(void *image);
extern unsigned int image_width(void *image);
extern unsigned int image_height(void *image);
extern unsigned int image_components(void *image);
/* BLOCK_WIDTH * BLOCK_HEIGHT.  Return value is actual number of samples. */
extern int image_hsi_of_block(void *image, int x, int y, float *hue, float *saturation, float *intensity);
extern unsigned char *image_luma(void *image);
extern unsigned char *image_rgb(void *image);
extern void *image_cinfo(void *image);
extern void *image_coeffp(void *image);

#endif
