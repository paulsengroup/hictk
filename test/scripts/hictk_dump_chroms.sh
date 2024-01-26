#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

echo "#############################"
echo "#### hictk dump (chroms) ####"

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
ref_hic="$data_dir/hic/4DNFIZ1ZVXC8.hic8"

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

if ! check_test_files_exist.sh "$ref_cooler"; then
  exit 1
fi

outdir="$(mktemp -d -t hictk-tmp-XXXXXXXXXX)"
trap 'rm -rf -- "$outdir"' EXIT

cooler dump -t chroms "$ref_cooler::/resolutions/100000" > "$outdir/expected.chrom.sizes"
"$hictk_bin" dump -t chroms "$ref_cooler::/resolutions/100000" > "$outdir/out.cooler.chrom.sizes"
"$hictk_bin" dump -t chroms "$ref_hic" > "$outdir/out.hic.chrom.sizes"

if ! compare_plain_files.sh "$outdir/expected.chrom.sizes" "$outdir/out.cooler.chrom.sizes"; then
  status=1
fi

if ! compare_plain_files.sh "$outdir/expected.chrom.sizes" "$outdir/out.hic.chrom.sizes"; then
  status=1
fi

if [ "$status" -eq 0 ]; then
  printf '\n### PASS ###\n'
else
  printf '\n### FAIL ###\n'
fi

exit "$status"
