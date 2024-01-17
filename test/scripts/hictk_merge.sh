#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u


echo "#####################"
echo "#### hictk merge ####"

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

input_cooler="$data_dir/integration_tests/4DNFIZ1ZVXC8.mcool"
input_hic="$data_dir/hic/4DNFIZ1ZVXC8.hic9"
resolution=100000

export PATH="$PATH:$script_dir"

if ! command -v cooler &> /dev/null; then
  2>&1 echo "Unable to find cooler in your PATH"
  status=1
fi

# Try to detect the error outlined below as early as possible:
# https://github.com/open2c/cooler/pull/298
cooler --help > /dev/null

if [ $status -ne 0 ]; then
  exit $status
fi

if ! check_test_files_exist.sh "$input_cooler"; then
  exit 1
fi

outdir="$(mktemp -d -t hictk-tmp-XXXXXXXXXX)"
trap 'rm -rf -- "$outdir"' EXIT

cooler merge "$outdir/expected.cool" "$input_cooler::/resolutions/$resolution" "$input_cooler::/resolutions/$resolution"

# Test merrging coolers
"$hictk_bin" merge "$input_cooler::/resolutions/$resolution" \
                   "$input_cooler::/resolutions/$resolution" \
                   -o "$outdir/out.cool" \
                   --chunk-size=9999
if ! compare_matrix_files.sh "$hictk_bin" "$outdir/expected.cool" "$outdir/out.cool" "$resolution"; then
  status=1
fi

# Test merging .hic
"$hictk_bin" merge "$input_hic" \
                   "$input_hic" \
                   -o "$outdir/out.hic" \
                   --resolution "$resolution" \
                   --chunk-size=9999
if ! compare_matrix_files.sh "$hictk_bin" "$outdir/expected.cool" "$outdir/out.hic" "$resolution"; then
  status=1
fi

if [ "$status" -eq 0 ]; then
  printf '\n### PASS ###\n'
else
  printf '\n### FAIL ###\n'
fi

exit "$status"
