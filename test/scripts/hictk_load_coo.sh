#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

echo "##########################"
echo "#### hictk load (COO) ####"

# readlink -f is not available on macos...
function readlink_py {
  set -eu
  python3 -c 'import os, sys; print(os.path.realpath(sys.argv[1]))' "$1"
}

function check_files_exist {
  set -eu
  status=0
  for f in "$@"; do
    if [ ! -f "$f" ]; then
      2>&1 echo "Unable to find test file \"$f\""
      status=1
    fi
  done

  return "$status"
}

function compare_files {
  set -o pipefail
  set -eu

  hictk="$1"
  resolution="${4}"
  f1="$2"
  f2="$3"

  2>&1 echo "Comparing $f1 with $f2..."
  if diff <("$hictk" dump --join "$f1"   \
                          --resolution   \
                          "$resolution") \
          <("$hictk" dump --join "$f2"   \
                          --resolution   \
                          "$resolution"); then
    2>&1 echo "Files are identical"
    return 0
  else
    2>&1 echo "Files differ"
    return 1
  fi
}

function shuffle {
  if command -v shuf &> /dev/null; then
    shuf
  else
    sort -R
  fi
}

export function readlink_py shuffle

status=0

if [ $# -ne 2 ]; then
  2>&1 echo "Usage: $0 path_to_hictk [un]sorted"
  status=1
fi

hictk_bin="$1"
if [[ "$2" == 'sorted' ]]; then
  sorted=true
else
  sorted=false
fi

data_dir="$(readlink_py "$(dirname "$0")/../data/integration_tests")"
script_dir="$(readlink_py "$(dirname "$0")")"

ref_cooler="$data_dir/4DNFIKNWM36K.subset.fixed-bins.cool"
resolution=10000
batch_size=999999

export PATH="$PATH:$script_dir"

if [ $status -ne 0 ]; then
  exit $status
fi

if ! check_files_exist "$ref_cooler"; then
  exit 1
fi

outdir="$(mktemp -d -t hictk-tmp-XXXXXXXXXX)"
trap 'rm -rf -- "$outdir"' EXIT

"$hictk_bin" dump -t chroms "$ref_cooler" > "$outdir/chrom.sizes"

if [[ "$sorted" == true ]]; then
  "$hictk_bin" dump -t pixels "$ref_cooler" |
    "$hictk_bin" load \
      -f coo \
      --assume-sorted \
      --batch-size "$batch_size" \
      --bin-size "$resolution" \
      --tmpdir "$outdir" \
      "$outdir/chrom.sizes" \
      "$outdir/out.cool"
else
  "$hictk_bin" dump -t pixels "$ref_cooler" |
     shuffle |
    "$hictk_bin" load \
      -f coo \
      --assume-unsorted \
      --batch-size "$batch_size" \
      --bin-size "$resolution" \
      --tmpdir "$outdir" \
      "$outdir/chrom.sizes" \
      "$outdir/out.cool"
fi

if ! compare_files "$hictk_bin" "$outdir/out.cool" "$ref_cooler" "$resolution"; then
  status=1
fi


if [[ "$sorted" == false ]]; then
  "$hictk_bin" dump -t pixels "$ref_cooler" |
    shuffle |
    "$hictk_bin" load \
      -f coo \
      --assume-unsorted \
      --batch-size "$batch_size" \
      --bin-size "$resolution" \
      --tmpdir "$outdir" \
      "$outdir/chrom.sizes" \
      "$outdir/out.hic"

  if ! compare_files "$hictk_bin" "$outdir/out.hic" "$ref_cooler" "$resolution"; then
    status=1
  fi
fi

if [ "$status" -eq 0 ]; then
  printf '\n### PASS ###\n'
else
  printf '\n### FAIL ###\n'
fi

exit "$status"
