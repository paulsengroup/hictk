#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

echo "##########################"
echo "#### hictk load (4DN) ####"

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

data_dir="$(readlink_py "$(dirname "$0")/../data/integration_tests")"
script_dir="$(readlink_py "$(dirname "$0")")"

pairs="$data_dir/4DNFIKNWM36K.subset.pairs.xz"
ref_cooler_fixed_bins="$data_dir/4DNFIKNWM36K.subset.fixed-bins.cool"
ref_cooler_variable_bins="$data_dir/4DNFIKNWM36K.subset.variable-bins.cool"

export PATH="$PATH:$script_dir"

if ! command -v xz &> /dev/null; then
  2>&1 echo "Unable to find xz in your PATH"
  status=1
fi

if [ $status -ne 0 ]; then
  exit $status
fi

if ! check_test_files_exist.sh "$pairs" "$ref_cooler_fixed_bins" "$ref_cooler_variable_bins"; then
  exit 1
fi

outdir="$(mktemp -d -t hictk-tmp-XXXXXXXXXX)"
trap 'rm -rf -- "$outdir"' EXIT

resolution=10000
batch_size=999999

cooler dump -t chroms "$ref_cooler_fixed_bins" > "$outdir/chrom.sizes"

# Test cooler with fixed bin size
xzcat "$pairs" |
  "$hictk_bin" load \
    -f 4dn \
    --chunk-size "$batch_size" \
    --bin-size "$resolution" \
    --tmpdir "$outdir" \
    "$outdir/chrom.sizes" \
    "$outdir/out.cool"

if ! compare_matrix_files.sh "$hictk_bin" "$outdir/out.cool" "$ref_cooler_fixed_bins" "$resolution"; then
  status=1
fi

# Test cooler with variable bin size
cooler dump -t bins "$ref_cooler_variable_bins" > "$outdir/bins.bed"

xzcat "$pairs" |
  "$hictk_bin" load \
    -f 4dn \
    --chunk-size "$batch_size" \
    --bin-table "$outdir/bins.bed" \
    --force \
    --tmpdir "$outdir" \
    "$outdir/chrom.sizes" \
    "$outdir/out.cool"

if ! compare_matrix_files.sh "$hictk_bin" "$outdir/out.cool" "$ref_cooler_variable_bins"; then
  status=1
fi


# Test hic with fixed bin size
xzcat "$pairs" |
  "$hictk_bin" load \
    -f 4dn \
    --chunk-size "$batch_size" \
    --bin-size "$resolution" \
    --tmpdir "$outdir" \
    "$outdir/chrom.sizes" \
    "$outdir/out.hic"

if ! compare_matrix_files.sh "$hictk_bin" "$outdir/out.hic" "$ref_cooler_fixed_bins" "$resolution"; then
  status=1
fi

if [ "$status" -eq 0 ]; then
  printf '\n### PASS ###\n'
else
  printf '\n### FAIL ###\n'
fi

exit "$status"
