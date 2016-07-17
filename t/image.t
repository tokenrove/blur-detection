#!/usr/bin/env bats

CMD=${CMD:-"valgrind --error-exitcode=42 --leak-check=full ./blur-detection"}

@test "#6: Ensure whole image is read" {
    $CMD t/corpus/hsi-bug-lena-wide.jpg | grep -q well-focused
}
