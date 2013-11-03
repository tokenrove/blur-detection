/**** DEBUG OUTPUT ****/

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "image.h"

void debug_dump_octave_matrix(char *name, float *m, int w, int h)
{
	FILE *out;

	out = fopen(name, "w");
	for(int i = 0; i < h; i++, fputc('\n', out))
		for(int j = 0; j < w; j++) fprintf(out, "%f%s", m[j+i*w], (j == w-1) ? "":"\t");
	fclose(out);
}

void debug_dump_ppm(char *path_, unsigned char *rgb, int w, int h, unsigned char *borders)
{
	FILE *out;
	unsigned char c;
	char *path, *extension = ".dump.ppm";

	path = malloc(strlen(path_) + strlen(extension));
	strcpy(path, path_);
	strcat(path, extension);
	out = fopen(path, "wb");
	free(path);
	fprintf(out, "P6 %d %d 255\n", w, h);
	for(int i = 0; i < h; i++)
		for(int j = 0; j < w; j++) {
			if(((i%BLOCK_HEIGHT) == 0 || (i%BLOCK_HEIGHT) == BLOCK_HEIGHT-1 ||
			    (j%BLOCK_WIDTH) == 0 || (j%BLOCK_WIDTH) == BLOCK_WIDTH-1) &&
			   (i/BLOCK_HEIGHT) < (h/BLOCK_HEIGHT)) {
				c = borders[(j/BLOCK_WIDTH)+(i/BLOCK_HEIGHT)*(w/BLOCK_WIDTH)];
				fputc(c, out); fputc(c, out); fputc(c, out);
				continue;
			}
			fputc(rgb[3*(j+i*w)], out);
			fputc(rgb[1+3*(j+i*w)], out);
			fputc(rgb[2+3*(j+i*w)], out);
		}
	fclose(out);
}

void debug_dump_pgm(char *path_, unsigned char *luma, int w, int h)
{
	FILE *out;
	char *path, *extension = ".dump.pgm";

	path = malloc(strlen(path_) + strlen(extension));
	strcpy(path, path_);
	strcat(path, extension);
	out = fopen(path, "wb");
	free(path);
	fprintf(out, "P5 %d %d 255\n", w, h);
	for(int i = 0; i < h; i++) for(int j = 0; j < w; j++) fputc(luma[j+i*w], out);
	fclose(out);
}
