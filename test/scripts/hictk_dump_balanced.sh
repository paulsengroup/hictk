#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

echo "###############################"
echo "#### hictk dump (balanced) ####"

# readlink -f is not available on macos...
function readlink_py {
  set -eu
  python3 -c 'import os, sys; print(os.path.realpath(sys.argv[1]))' "$1"
}

function truncate_counts {
  set -o pipefail
  set -e

  awk -F '\t' '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%.3f\n", $1, $2, $3, $4, $5, $6, $7}'

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

cooler dump --balanced --na-rep nan --join "$ref_cooler::/resolutions/100000" -r chr2L | cut -f 1-6,8 | truncate_counts > "$outdir/expected.pixels"
"$hictk_bin" dump --join --balance "weight" "$ref_cooler::/resolutions/100000" -r chr2L | truncate_counts > "$outdir/out.cooler.pixels"

if ! compare_plain_files.sh "$outdir/expected.pixels" "$outdir/out.cooler.pixels"; then
  status=1
fi

if [ "$status" -eq 0 ]; then
  printf '\n### PASS ###\n'
else
  printf '\n### FAIL ###\n'
fi

exit "$status"
