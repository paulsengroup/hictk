#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

echo "############################"
echo "#### hictk dump (trans) ####"

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

ref_cooler="$data_dir/integration_tests/4DNFIZ1ZVXC8.mcool"
ref_hic8="$data_dir/hic/4DNFIZ1ZVXC8.hic8"
ref_hic9="$data_dir/hic/4DNFIZ1ZVXC8.hic9"

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

if ! check_test_files_exist.sh "$ref_cooler" "$ref_hic8" "$ref_hic9"; then
  exit 1
fi

outdir="$(mktemp -d -t hictk-tmp-XXXXXXXXXX)"
trap 'rm -rf -- "$outdir"' EXIT

# Test chrom-chrom matrix
cooler dump --join "$ref_cooler::/resolutions/100000" --range chr2L --range2 chrX > "$outdir/expected.pixels"
"$hictk_bin" dump --join "$ref_cooler::/resolutions/100000" --range chr2L --range2 chrX > "$outdir/out.cooler.pixels"
"$hictk_bin" dump --join --resolution 100000 "$ref_hic8" --range chr2L --range2 chrX > "$outdir/out.hic8.pixels"
"$hictk_bin" dump --join --resolution 100000 "$ref_hic9" --range chr2L --range2 chrX > "$outdir/out.hic9.pixels"

if ! compare_plain_files.sh "$outdir/expected.pixels" "$outdir/out.cooler.pixels"; then
  status=1
fi

if ! compare_plain_files.sh "$outdir/expected.pixels" "$outdir/out.hic8.pixels"; then
  status=1
fi

if ! compare_plain_files.sh "$outdir/expected.pixels" "$outdir/out.hic9.pixels"; then
  status=1
fi

# Test --trans-only matrix (sorted)
cooler dump --join "$ref_cooler::/resolutions/100000" | awk -F '\t' '$1!=$4' | tee "$outdir/expected.pixels" > /dev/null
"$hictk_bin" dump --join "$ref_cooler::/resolutions/100000" --trans-only | tee "$outdir/out.cooler.pixels" > /dev/null
"$hictk_bin" dump --join --resolution 100000 "$ref_hic8" --trans-only | tee "$outdir/out.hic8.pixels" > /dev/null
"$hictk_bin" dump --join --resolution 100000 "$ref_hic9" --trans-only | tee "$outdir/out.hic9.pixels" > /dev/null

if ! compare_plain_files.sh "$outdir/expected.pixels" "$outdir/out.cooler.pixels"; then
  status=1
fi

if ! compare_plain_files.sh "$outdir/expected.pixels" "$outdir/out.hic8.pixels"; then
  status=1
fi

if ! compare_plain_files.sh "$outdir/expected.pixels" "$outdir/out.hic9.pixels"; then
  status=1
fi

# Test --trans-only matrix (unsorted)
"$hictk_bin" dump --join "$ref_cooler::/resolutions/100000" --trans-only --unsorted | sort -V | tee "$outdir/out.cooler.pixels" > /dev/null
"$hictk_bin" dump --join --resolution 100000 "$ref_hic8" --trans-only --unsorted | sort -V | tee "$outdir/out.hic8.pixels" > /dev/null
"$hictk_bin" dump --join --resolution 100000 "$ref_hic9" --trans-only --unsorted | sort -V | tee "$outdir/out.hic9.pixels" > /dev/null

if ! compare_plain_files.sh "$outdir/expected.pixels" "$outdir/out.cooler.pixels"; then
  status=1
fi

if ! compare_plain_files.sh "$outdir/expected.pixels" "$outdir/out.hic8.pixels"; then
  status=1
fi

if ! compare_plain_files.sh "$outdir/expected.pixels" "$outdir/out.hic9.pixels"; then
  status=1
fi

if [ "$status" -eq 0 ]; then
  printf '\n### PASS ###\n'
else
  printf '\n### FAIL ###\n'
fi

exit "$status"
