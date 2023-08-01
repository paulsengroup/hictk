#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

echo "########################"
echo "#### hictk validate ####"

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

status=0

if [ $# -ne 1 ]; then
  2>&1 echo "Usage: $0 path_to_hictk"
  status=1
fi

hictk_bin="$1"

data_dir="$(readlink_py "$(dirname "$0")/../data/")"
script_dir="$(readlink_py "$(dirname "$0")")"

valid_hic="$data_dir/hic/4DNFIZ1ZVXC8.hic8"
valid_cooler="$data_dir/cooler/ENCFF993FGR.2500000.cool"
valid_mcool="$data_dir/cooler/multires_cooler_test_file.mcool"
valid_scool="$data_dir/cooler/single_cell_cooler_test_file.scool"
invalid_cooler="$data_dir/cooler/invalid_coolers/4DNFI9GMP2J8.1000000.cool"

export PATH="$PATH:$script_dir"

if [ $status -ne 0 ]; then
  exit $status
fi

if ! check_files_exist "$valid_hic" "$valid_cooler" "$valid_mcool" "$valid_scool" "$invalid_cooler"; then
  exit 1
fi

outdir="$(mktemp -d -t hictk-tmp-XXXXXXXXXX)"
trap 'rm -rf -- "$outdir"' EXIT

echo "# Validating $valid_hic..."
if ! "$hictk_bin" validate "$valid_hic" >> "$outdir/out.txt"; then
  status=1
fi

echo "# Validating $valid_cooler..."
if ! "$hictk_bin" validate --validate-index "$valid_cooler" >> "$outdir/out.txt"; then
  status=1
fi

echo "# Validating $valid_mcool..."
if ! "$hictk_bin" validate --validate-index "$valid_mcool" >> "$outdir/out.txt"; then
  status=1
fi

echo "# Validating $valid_scool..."
if ! "$hictk_bin" validate --validate-index "$valid_scool" >> "$outdir/out.txt"; then
  status=1
fi

echo "# Validating $invalid_cooler..."
if "$hictk_bin" validate --validate-index "$invalid_cooler" >> "$outdir/out.txt"; then
  status=1
fi

if [ "$status" -eq 0 ]; then
  printf '\n### PASS ###\n'
else
  cat "$outdir/out.txt"
  printf '\n### FAIL ###\n'
fi

exit "$status"
