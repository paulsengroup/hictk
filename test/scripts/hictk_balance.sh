#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

echo "##################################"
echo "#### hictk balance            ####"

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

function dump_interactions {
  set -o pipefail
  set -eu

  hictk="$1"
  resolution="$3"
  f="$2"

  if [[ "$f"  == *.hic ]]; then
    weight=WEIGHT
  else
    weight=weight
  fi

  "$hictk" dump "$f"                \
                --balance="$weight" \
                --resolution        \
                "$resolution"       |
                cut -f 3
}

function absolute_error {
  set -o pipefail
  set -eu

  f1="$1"
  f2="$2"

  # shellcheck disable=SC2016
  cmd='function abs(v) {
          return v < 0 ? -v : v
       }
       ($2!=$1 && abs($1 - $2) > 1.0e-5) { print $0 }
       '

  # Fail if the absolute error is > 1.0e-5
  if paste "$f1" "$f2" | awk -F '\t' "$cmd" | grep . ; then
    return 1
  else
    return 0
  fi
}

function compare_matrices {
  set -o pipefail
  set -eu

  hictk="$1"
  resolution="$4"
  f1="$2"
  f2="$3"

  2>&1 echo "Comparing $f1 with $f2..."
  if absolute_error \
      <(dump_interactions "$hictk" "$f1" "$resolution") \
      <(dump_interactions "$hictk" "$f2" "$resolution"); then
    2>&1 echo "Files are identical"
    return 0
  else
    2>&1 echo "Files differ"
    return 1
  fi
}

export function readlink_py

status=0

if [ $# -ne 2 ]; then
  2>&1 echo "Usage: $0 path_to_hictk juicer_tools.jar"
  status=1
fi

hictk_bin="$1"
juicer_tools_jar="$2"

data_dir="$(readlink_py "$(dirname "$0")/../data/")"
script_dir="$(readlink_py "$(dirname "$0")")"

ref_cool="$data_dir/cooler/ENCFF993FGR.2500000.cool"
ref_hic="$data_dir/hic/ENCFF993FGR.hic"

export PATH="$PATH:$script_dir"

if ! check_files_exist "$ref_cool" "$ref_hic" "$juicer_tools_jar"; then
  exit 1
fi

outdir="$(mktemp -d -t hictk-tmp-XXXXXXXXXX)"
trap 'rm -rf -- "$outdir"' EXIT

cp "$ref_cool" "$ref_hic" "$outdir"

"$hictk_bin" balance "$outdir/"*.cool -t $(nproc) --chunk-size=100 --mode=cis --force
if ! compare_matrices "$hictk_bin" "$outdir/"*.cool "$ref_cool" 2500000; then
  status=1
fi

"$hictk_bin" balance "$outdir/"*.hic -t $(nproc) --chunk-size=100 --mode=cis --force --juicer-tools-jar "$juicer_tools_jar"
if ! compare_matrices "$hictk_bin" "$outdir/"*.hic "$ref_cool" 2500000; then
  status=1
fi

"$hictk_bin" balance "$outdir/"*.cool -t $(nproc) --in-memory --mode=cis --force
if ! compare_matrices "$hictk_bin" "$outdir/"*.cool "$ref_cool" 2500000; then
  status=1
fi

if [ "$status" -eq 0 ]; then
  printf '\n### PASS ###\n'
else
  printf '\n### FAIL ###\n'
fi

exit "$status"
