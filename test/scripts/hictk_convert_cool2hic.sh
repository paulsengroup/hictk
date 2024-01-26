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

ref_cool="$data_dir/integration_tests/4DNFIZ1ZVXC8.mcool"

export PATH="$PATH:$script_dir"

if ! check_test_files_exist.sh "$ref_cool"; then
  exit 1
fi

outdir="$(mktemp -d -t hictk-tmp-XXXXXXXXXX)"
trap 'rm -rf -- "$outdir"' EXIT

resolutions=(100000 2500000)

"$hictk_bin" convert \
             "$ref_cool" \
             "$outdir/out.hic" \
             --resolutions ${resolutions[*]} \
             --threads "$(nproc.sh)" \
             --chunk-size 100000 \
             --compression-lvl 1

for resolution in "${resolutions[@]}"; do
  if ! compare_matrix_files.sh "$hictk_bin_opt" "$outdir/out.hic" "$ref_cool" "$resolution"; then
    status=1
  fi
done

"$hictk_bin" dump -t normalizations "$ref_cool" | sed 's/weight/ICE/' | sort > "$outdir/normalizations.mcool"
"$hictk_bin" dump -t normalizations "$outdir/out.hic" | sort > "$outdir/normalizations.hic"

if ! compare_plain_files.sh "$outdir/normalizations.mcool" "$outdir/normalizations.hic"; then
  status=1
fi

if [ "$status" -eq 0 ]; then
  printf '\n### PASS ###\n'
else
  printf '\n### FAIL ###\n'
fi

exit "$status"
