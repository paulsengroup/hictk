#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

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

function compare_coolers {
  set -o pipefail
  set -eu

  hictk="$1"
  resolution="$4"
  clr="$2::/resolutions/$resolution"
  hic="$3"

  2>&1 echo "Comparing $clr with $hic..."
  if diff <("$hictk" dump --join "$clr")  \
          <("$hictk" dump --join "$hic"   \
                          --resolution    \
                          "$resolution");  then
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

hic="$data_dir/hic/4DNFIZ1ZVXC8.hic9"

export PATH="$PATH:$script_dir"

if ! check_files_exist "$hic"; then
  exit 1
fi

outdir="$(mktemp -d -t hictk-tmp-XXXXXXXXXX)"
trap 'rm -rf -- "$outdir"' EXIT

resolutions=(50000 2500000)

"$hictk_bin" convert \
             "$hic" \
             "$outdir/out.mcool" \
             --resolutions ${resolutions[*]}

for resolution in "${resolutions[@]}"; do
  if ! compare_coolers "$hictk_bin" "$outdir/out.mcool" "$hic" "$resolution"; then
    status=1
  fi
done

if [ "$status" -eq 0 ]; then
  printf '\n### PASS ###\n'
else
  printf '\n### FAIL ###\n'
fi

exit "$status"
