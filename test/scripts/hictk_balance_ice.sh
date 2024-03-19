#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

echo "##################################"
echo "#### hictk balance (ICE)      ####"

# readlink -f is not available on macos...
function readlink_py {
  set -eu
  python3 -c 'import os, sys; print(os.path.realpath(sys.argv[1]))' "$1"
}

function dump_interactions {
  set -o pipefail
  set -eu

  hictk="$1"
  resolution="$3"
  f="$2"

  "$hictk" dump "$f"             \
                --balance=weight \
                --resolution     \
                "$resolution"    |
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

data_dir="$(readlink_py "$(dirname "$0")/../data/")"
script_dir="$(readlink_py "$(dirname "$0")")"

ref_cool="$data_dir/cooler/ENCFF993FGR.2500000.cool"
ref_hic="$data_dir/hic/ENCFF993FGR.2500000.hic"

export PATH="$PATH:$script_dir"

if ! check_test_files_exist.sh "$ref_cool" "$ref_hic"; then
  exit 1
fi

outdir="$(mktemp -d -t hictk-tmp-XXXXXXXXXX)"
trap 'rm -rf -- "$outdir"' EXIT

cp "$ref_cool" "$ref_hic" "$outdir"

"$hictk_bin" balance ice "$outdir/"*.cool   \
                         -t $(nproc.sh)     \
                         --chunk-size=100   \
                         --mode=cis         \
                         --tmpdir="$outdir" \
                         --force
if ! compare_matrices "$hictk_bin_opt" "$outdir/"*.cool "$ref_cool" 2500000; then
  status=1
fi

"$hictk_bin" balance ice "$outdir/"*.hic    \
                         -t $(nproc.sh)     \
                         --chunk-size=100   \
                         --mode=cis         \
                         --tmpdir="$outdir" \
                         --name=weight      \
                         --force
if ! compare_matrices "$hictk_bin_opt" "$outdir/"*.hic "$ref_cool" 2500000; then
  status=1
fi

"$hictk_bin" balance ice "$outdir/"*.cool   \
                         -t $(nproc.sh)     \
                         --in-memory        \
                         --mode=cis         \
                         --tmpdir="$outdir" \
                         --force
if ! compare_matrices "$hictk_bin_opt" "$outdir/"*.cool "$ref_cool" 2500000; then
  status=1
fi

if [ "$status" -eq 0 ]; then
  printf '\n### PASS ###\n'
else
  printf '\n### FAIL ###\n'
fi

exit "$status"
