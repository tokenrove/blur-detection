#!/usr/bin/env bats
#
# Tests known focused/blurred results from a corpus

CMD=${CMD:-"valgrind --error-exitcode=42 --leak-check=full ./blur-detection"}

expect_sharp() {
    $CMD "$1" | grep -q well-focused
}

expect_blurred() {
    $CMD "$1" | grep -q ill-focused
}

@test "known sharp" {
    for i in t/corpus/sharp/*; do
        expect_sharp $i || exit 1
    done
}

@test "known blurred" {
    for i in t/corpus/blur/*; do
        expect_blurred $i
    done
}
