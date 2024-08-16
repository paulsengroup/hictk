#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

echo "##############################"
echo "#### hictk dump (weights) ####"

# readlink -f is not available on macos...
function readlink_py {
  set -eu
  python3 -c 'import os, sys; print(os.path.realpath(sys.argv[1]))' "$1"
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
       ($2!=$1 && abs($1 - $2) > 1.0e-4) { print $0 }
       '

  # Fail if the absolute error is > 1.0e-5
  if paste "$f1" "$f2" | awk -F '\t' "$cmd" | grep . ; then
    return 1
  else
    return 0
  fi
}


export function readlink_py absolute_error

status=0

if [ $# -ne 1 ]; then
  2>&1 echo "Usage: $0 path_to_hictk"
  status=1
fi

hictk_bin="$1"

data_dir="$(readlink_py "$(dirname "$0")/../data/")"
script_dir="$(readlink_py "$(dirname "$0")")"

ref_cooler="$data_dir/cooler/ENCFF993FGR.2500000.cool"
ref_hic="$data_dir/hic/ENCFF993FGR.2500000.hic"

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

if ! check_test_files_exist.sh "$ref_cooler" "$ref_hic"; then
  exit 1
fi

outdir="$(mktemp -d -t hictk-tmp-XXXXXXXXXX)"
trap 'rm -rf -- "$outdir"' EXIT

cooler dump -t bins "$ref_cooler" | cut -f 4 > "$outdir/expected1.weights"

"$hictk_bin" dump -t weights "$ref_cooler" | tail -n+2 | cut -f 1 > "$outdir/out.cooler.weights"
"$hictk_bin" dump -t weights --resolution 2500000 "$ref_hic" | tail -n+2 | cut -f 1 > "$outdir/out.hic.weights"

if ! absolute_error "$outdir/expected1.weights" "$outdir/out.cooler.weights"; then
  status=1
fi

if ! absolute_error "$outdir/expected1.weights" "$outdir/out.hic.weights"; then
  status=1
fi

if [ "$status" -eq 0 ]; then
  printf '\n### PASS ###\n'
else
  printf '\n### FAIL ###\n'
fi

exit "$status"
