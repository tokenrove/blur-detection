#!/usr/bin/env bats

CMD=${CMD:-"valgrind --error-exitcode=42 --leak-check=full ./blur-detection"}

@test "Bug: uninitialized values when sky threshold triggered" {
    $CMD -k=0 t/corpus/sky-threshold-bug-blurred-circle.jpg
}
