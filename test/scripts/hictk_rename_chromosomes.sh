#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

echo "##################################"
echo "#### hictk rename-chromosomes ####"

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

input_cooler="$data_dir/cooler/cooler_test_file.cool"
input_mcool="$data_dir/cooler/multires_cooler_test_file.mcool"
input_scool="$data_dir/cooler/single_cell_cooler_test_file.scool"

export PATH="$PATH:$script_dir"

if [ $status -ne 0 ]; then
  exit $status
fi

if ! check_test_files_exist.sh "$input_cooler" "$input_mcool" "$input_scool"; then
  exit 1
fi

outdir="$(mktemp -d -t hictk-tmp-XXXXXXXXXX)"
trap 'rm -rf -- "$outdir"' EXIT

# Test adding chr prefix
cp "$input_cooler" "$outdir/out1.cool"
"$hictk_bin" rename-chroms "$outdir/out1.cool" --add-chr-prefix
if ! "$hictk_bin" dump -t chroms "$outdir/out1.cool" | grep -q chr ; then
  status=1
fi

# Test removing chr prefix
cp "$outdir/out1.cool" "$outdir/out2.cool"
"$hictk_bin" rename-chroms "$outdir/out2.cool" --remove-chr-prefix
if ! "$hictk_bin" dump -t chroms "$outdir/out2.cool" | grep -vq chr ; then
  status=1
fi


# Test renaming using custom name mappings
printf '1\tABC\n' > "$outdir/mappings.txt"
cp "$input_cooler" "$outdir/out3.cool"
"$hictk_bin" rename-chroms "$outdir/out3.cool" --name-mappings "$outdir/mappings.txt"
if ! "$hictk_bin" dump -t chroms "$outdir/out3.cool" | grep -q ABC ; then
  status=1
fi


# Test mcool
cp "$input_mcool" "$outdir/out4.mcool"
"$hictk_bin" rename-chroms "$outdir/out4.mcool" --name-mappings "$outdir/mappings.txt"
if ! "$hictk_bin" dump -t chroms "$outdir/out4.mcool" | grep -q ABC ; then
  status=1
fi


# Test scool
cp "$input_scool" "$outdir/out5.scool"
"$hictk_bin" rename-chroms "$outdir/out5.scool" --name-mappings "$outdir/mappings.txt"
if ! "$hictk_bin" dump -t chroms "$outdir/out5.scool" | grep -q ABC ; then
  status=1
fi


if [ "$status" -eq 0 ]; then
  printf '\n### PASS ###\n'
else
  printf '\n### FAIL ###\n'
fi

exit "$status"
