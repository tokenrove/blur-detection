#ifndef TONG_LI_ZHANG_ZHANG_H
#define TONG_LI_ZHANG_ZHANG_H

/* Note that detect_blur() must be changed if the number of levels has
 * changed. */
enum { LEVELS = 3 };

extern void detect_blur(float **emap, int w, int h, float threshold, float *da_ratio, float *blur_extent);
extern void construct_edge_map(float *dst, float **emap, int w, int h);
extern void haar_transform(float *data, int w, int h);

#endif
