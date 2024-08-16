#!/usr/bin/env bash

# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u


echo "########################"
echo "#### hictk metadata ####"

# readlink -f is not available on macos...
function readlink_py {
  set -eu
  python3 -c 'import os, sys; print(os.path.realpath(sys.argv[1]))' "$1"
}

export function readlink_py

status=0

if [ $# -ne 1 ]; then
  2>&1 echo "Usage: $0 path_to_hictk"
  status=1
fi

hictk_bin="$1"

data_dir="$(readlink_py "$(dirname "$0")/../data/")"
script_dir="$(readlink_py "$(dirname "$0")")"


test_files=(
  "$data_dir/cooler/cooler_test_file.cool"
  "$data_dir/cooler/multires_cooler_test_file.mcool"
  "$data_dir/cooler/single_cell_cooler_test_file.scool"
  "$data_dir/hic/4DNFIZ1ZVXC8.hic8"
  "$data_dir/hic/4DNFIZ1ZVXC8.hic9"
)

export PATH="$PATH:$script_dir"

if [ $status -ne 0 ]; then
  exit $status
fi

if ! check_test_files_exist.sh "${test_files[@]}"; then
  exit 1
fi


for file in "${test_files[@]}"; do
  for fmt in json toml yaml; do
    test_script="is_valid_${fmt}.py"
    echo "Testing $(basename "$file") [$(echo "$fmt" | tr '[:lower:]' '[:upper:]')]"
    if ! "$hictk_bin" metadata "$file" -f "$fmt" | "$test_script" -q; then
      status=1
    fi

    if ! "$hictk_bin" metadata "$file" -f "$fmt" --recursive | "$test_script" -q; then
      status=1
    fi
  done
done

if [ "$status" -eq 0 ]; then
  printf '\n### PASS ###\n'
else
  printf '\n### FAIL ###\n'
fi

exit "$status"
