#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

echo "###########################"
echo "#### hictk dump (normalizations) ####"

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

expected_norms_hic=(
  KR
  SCALE
  VC
  VC_SQRT
)

expected_norms_cooler=(
  KR
  SCALE
  VC
  VC_SQRT
  weight
)

export PATH="$PATH:$script_dir"

if ! check_test_files_exist.sh "$mclr" "$sclr" "$hic"; then
  exit 1
fi

outdir="$(mktemp -d -t hictk-tmp-XXXXXXXXXX)"
trap 'rm -rf -- "$outdir"' EXIT

printf "%s\n" "${expected_norms_hic[@]}" > "$outdir/expected.hic.txt"
printf "%s\n" "${expected_norms_cooler[@]}" > "$outdir/expected.cool.txt"

"$hictk_bin" dump -t normalizations "$mclr" > "$outdir/mcool.norms.txt"
"$hictk_bin" dump -t normalizations "$sclr" > "$outdir/scool.norms.empty.txt"
"$hictk_bin" dump -t normalizations "$mclr::/resolutions/100000" > "$outdir/cool.norms.txt"

"$hictk_bin" dump -t normalizations "$hic" > "$outdir/hic.norms.txt"
"$hictk_bin" dump -t normalizations "$hic" --resolution 100000 > "$outdir/hic.norms.txt"

for f in "$outdir/"*cool*.norms.txt; do
  if ! compare_plain_files.sh "$outdir/expected.cool.txt" "$f"; then
    status=1
  fi
done

for f in "$outdir/"*hic*.norms.txt; do
  if ! compare_plain_files.sh "$outdir/expected.hic.txt" "$f"; then
    status=1
  fi
done

for f in "$outdir/"*.norms.empty.txt; do
  2>&1 echo "Checking that $f is empty..."
  if [ -s "$f" ]; then
    2>&1 echo "File is NOT empty"
    status=1
  else
    2>&1 echo "File is empty"
  fi
done

if [ "$status" -eq 0 ]; then
  printf '\n### PASS ###\n'
else
  printf '\n### FAIL ###\n'
fi

exit "$status"
