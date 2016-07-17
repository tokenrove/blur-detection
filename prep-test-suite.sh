#!/bin/sh
#
# This script generates synthetic blurred versions of existing test
# corpus files.

set -eu

mkdir -p t/corpus/blur
convert t/corpus/sharp/lena_rgb.jpg -blur 0x3 t/corpus/blur/lena_rgb_blur_0x3.jpg
touch t/.prepared
