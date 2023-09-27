#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

echo "###########################"
echo "#### hictk dump (resolutions) ####"

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
  set -e

  2>&1 echo "Comparing $1 with $2..."
  if diff "$1" "$2"; then
    2>&1 echo "Files are identical"
    return 0
  else
    2>&1 echo "Files differ"
    return 1
  fi
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

mclr="$data_dir/integration_tests/4DNFIZ1ZVXC8.mcool"
sclr="$data_dir/cooler/single_cell_cooler_test_file.scool"
hic="$data_dir/hic/4DNFIZ1ZVXC8.hic8"

expected_res=(
  1000
  5000
  10000
  25000
  50000
  100000
  250000
  500000
  1000000
  2500000
)

export PATH="$PATH:$script_dir"

if ! check_files_exist "$mclr" "$sclr" "$hic"; then
  exit 1
fi

outdir="$(mktemp -d -t hictk-tmp-XXXXXXXXXX)"
trap 'rm -rf -- "$outdir"' EXIT

printf "%s\n" "${expected_res[@]}" > "$outdir/expected.txt"
printf "%d\n" "100000" > "$outdir/expected.100000.txt"
"$hictk_bin" dump -t resolutions "$mclr" > "$outdir/mcool.res.txt"
"$hictk_bin" dump -t resolutions "$sclr" > "$outdir/scool.res.100000.txt"
"$hictk_bin" dump -t resolutions "$mclr::/resolutions/100000" > "$outdir/cool.res.100000.txt"
"$hictk_bin" dump -t resolutions "$hic" > "$outdir/hic.res.txt"
"$hictk_bin" dump -t resolutions "$hic" --resolution 100000 > "$outdir/hic.res.100000.txt"


for f in "$outdir/"*.res.txt; do
  if ! compare_files "$outdir/expected.txt" "$f"; then
    status=1
  fi
done

for f in "$outdir/"*.res.100000.txt; do
  if ! compare_files "$outdir/expected.100000.txt" "$f"; then
    status=1
  fi
done

if [ "$status" -eq 0 ]; then
  printf '\n### PASS ###\n'
else
  printf '\n### FAIL ###\n'
fi

exit "$status"
