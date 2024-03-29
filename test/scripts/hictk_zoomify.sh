#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

echo "#######################"
echo "#### hictk zoomify ####"

# readlink -f is not available on macos...
function readlink_py {
  set -eu
  python3 -c 'import os, sys; print(os.path.realpath(sys.argv[1]))' "$1"
}

function check_resolution_available {
  set -eu

  hictk_bin_="$1"
  file="$2"
  resolution="$3"

  if ! "$hictk_bin_" dump -t resolutions "$file" |
    awk "END{exit !(NR == 1 && /${resolution}/)}"; then
      echo 1>&2 "ERROR: Unable to find resolution ${resolution}!"
      return 1
  fi
}

export function readlink_py check_resolution_available

status=0

if [ $# -ne 1 ]; then
  2>&1 echo "Usage: $0 path_to_hictk"
  status=1
fi

hictk_bin="$1"
hictk_bin_opt="$(which hictk 2> /dev/null || true)"
if [ -z "$hictk_bin_opt" ]; then
  hictk_bin_opt="$hictk_bin"
fi

data_dir="$(readlink_py "$(dirname "$0")/../data/integration_tests")"
script_dir="$(readlink_py "$(dirname "$0")")"

ref_cooler="$data_dir/4DNFIZ1ZVXC8.mcool"
ref_hic="$data_dir/4DNFIZ1ZVXC8.500000.hic"
resolutions=(50000 100000 250000 2500000)

export PATH="$PATH:$script_dir"

if ! check_test_files_exist.sh "$ref_cooler" "$ref_hic"; then
  exit 1
fi

outdir="$(mktemp -d -t hictk-tmp-XXXXXXXXXX)"
trap 'rm -rf -- "$outdir"' EXIT

# Test cooler (multiple resolutions)
"$hictk_bin" zoomify \
  -t $(nproc.sh) \
  "$ref_cooler::/resolutions/${resolutions[0]}" \
  --compression-lvl 1 \
  "$outdir/out.mcool"

for res in "${resolutions[@]}"; do
  if ! compare_matrix_files.sh "$hictk_bin_opt" "$outdir/out.mcool" "$ref_cooler" "$res"; then
    status=1
  fi
done

# Test cooler (single resolution)
"$hictk_bin" zoomify \
  "$ref_cooler::/resolutions/${resolutions[0]}" \
  "$outdir/out.cool" \
  -t $(nproc.sh) \
  --no-copy-base-resolution \
  --compression-lvl 1 \
  --resolutions "${resolutions[1]}"

if ! check_resolution_available "$hictk_bin_opt" "$outdir/out.cool" "${resolutions[1]}"; then
  status=1
fi

if ! compare_matrix_files.sh "$hictk_bin_opt" "$outdir/out.cool" "$ref_cooler" "${resolutions[1]}"; then
  status=1
fi

# Test hic (multiple resolutions)
"$hictk_bin" zoomify \
  "$ref_hic" \
  "$outdir/out.hic" \
  -t $(nproc.sh) \
  --compression-lvl 1

for res in "${resolutions[@]}"; do
  if ! compare_matrix_files.sh "$hictk_bin_opt" "$outdir/out.hic" "$ref_cooler" "$res"; then
    status=1
  fi
done

# Test hic (single resolution)
"$hictk_bin" zoomify \
  "$ref_hic" \
  "$outdir/out.hic" \
  -t $(nproc.sh) \
  --compression-lvl 1 \
  --no-copy-base-resolution \
  --resolutions "${resolutions[1]}" \
  --force


if ! check_resolution_available "$hictk_bin_opt" "$outdir/out.hic" "${resolutions[1]}" ; then
  status=1
fi

if ! compare_matrix_files.sh "$hictk_bin_opt" "$outdir/out.hic" "$ref_cooler" "${resolutions[1]}"; then
  status=1
fi

if [ "$status" -eq 0 ]; then
  printf '\n### PASS ###\n'
else
  printf '\n### FAIL ###\n'
fi

exit "$status"
