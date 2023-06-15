#!/bin/bash

host=$(hostname);
my_path="$(dirname -- ${BASH_SOURCE[0]})"

"$my_path/../bin/$host/run_abc" "$@"
