#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

echo "##################################"
echo "#### hictk convert (cool2hic) ####"

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
  resolution="$4"
  hic="$2"
  clr="$3::/resolutions/$resolution"

  2>&1 echo "Comparing $hic with $clr..."
  if diff <("$hictk" dump --join "$hic"   \
                          --resolution    \
                          "$resolution")  \
          <("$hictk" dump --join "$clr"); then
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

ref_cool="$data_dir/integration_tests/4DNFIZ1ZVXC8.mcool"

export PATH="$PATH:$script_dir"

if ! check_files_exist "$ref_cool"; then
  exit 1
fi

outdir="$(mktemp -d -t hictk-tmp-XXXXXXXXXX)"
trap 'rm -rf -- "$outdir"' EXIT

resolutions=(100000 2500000)

"$hictk_bin" convert \
             "$ref_cool" \
             "$outdir/out.hic" \
             --resolutions ${resolutions[*]} \
             --threads "$(nproc)" \
             --batch-size 100000

for resolution in "${resolutions[@]}"; do
  if ! compare_files "$hictk_bin" "$outdir/out.hic" "$ref_cool" "$resolution"; then
    status=1
  fi
done

if [ "$status" -eq 0 ]; then
  printf '\n### PASS ###\n'
else
  printf '\n### FAIL ###\n'
fi

exit "$status"
