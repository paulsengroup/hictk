#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

echo "###########################"
echo "#### hictk dump (cells) ####"

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

expected_cells=(
  GSM2687248_41669_ACAGTG-R1-DpnII.100000.cool
  GSM2687249_41670_GGCTAC-R1-DpnII.100000.cool
  GSM2687250_41671_TTAGGC-R1-DpnII.100000.cool
  GSM2687251_41672_AGTTCC-R1-DpnII.100000.cool
  GSM2687252_41673_CCGTCC-R1-DpnII.100000.cool
)

export PATH="$PATH:$script_dir"

if ! check_files_exist "$mclr" "$sclr" "$hic"; then
  exit 1
fi

outdir="$(mktemp -d -t hictk-tmp-XXXXXXXXXX)"
trap 'rm -rf -- "$outdir"' EXIT

printf "%s\n" "${expected_cells[@]}" > "$outdir/expected.txt"


"$hictk_bin" dump -t cells "$sclr" > "$outdir/scool.cells.txt"

if ! compare_files "$outdir/expected.txt" "$outdir/scool.cells.txt"; then
  status=1
fi

if ! "$hictk_bin" dump -t cells "$mclr" &> /dev/null; then
  2>&1 echo "hictk dump -t resolution $mclr: OK"
else
  2>&1 echo "hictk dump -t resolution $mclr: FAIL"
  status=1
fi

if ! "$hictk_bin" dump -t cells "$hic" &> /dev/null; then
  2>&1 echo "hictk dump -t resolution $hic: OK"
else
  2>&1 echo "hictk dump -t resolution $hic: FAIL"
  status=1
fi

if ! "$hictk_bin" dump -t cells "$hic" --resolution 100000 &> /dev/null; then
  2>&1 echo "hictk dump -t resolution $hic --resolution 1000000: OK"
else
  2>&1 echo "hictk dump -t resolution $hic --resolution 1000000: FAIL"
  status=1
fi


if [ "$status" -eq 0 ]; then
  printf '\n### PASS ###\n'
else
  printf '\n### FAIL ###\n'
fi

exit "$status"
