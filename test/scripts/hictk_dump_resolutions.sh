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

if ! check_test_files_exist.sh "$mclr" "$sclr" "$hic"; then
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
  if ! compare_plain_files.sh "$outdir/expected.txt" "$f"; then
    status=1
  fi
done

for f in "$outdir/"*.res.100000.txt; do
  if ! compare_plain_files.sh "$outdir/expected.100000.txt" "$f"; then
    status=1
  fi
done

if [ "$status" -eq 0 ]; then
  printf '\n### PASS ###\n'
else
  printf '\n### FAIL ###\n'
fi

exit "$status"
