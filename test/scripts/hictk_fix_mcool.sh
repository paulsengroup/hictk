#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

echo "#########################"
echo "#### hictk fix-mcool ####"

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
hictk_bin_opt="$(which hictk)"
if [ -z "$hictk_bin_opt" ]; then
  hictk_bin_opt="$hictk_bin"
fi

data_dir="$(readlink_py "$(dirname "$0")/../data/")"
script_dir="$(readlink_py "$(dirname "$0")")"

invalid_mcool="$data_dir/cooler/invalid_coolers/corrupted_index.mcool"

export PATH="$PATH:$script_dir"

if [ $status -ne 0 ]; then
  exit $status
fi

if ! check_test_files_exist.sh "$invalid_mcool"; then
  exit 1
fi

outdir="$(mktemp -d -t hictk-tmp-XXXXXXXXXX)"
trap 'rm -rf -- "$outdir"' EXIT

"$hictk_bin" fix-mcool "$invalid_mcool" "$outdir/out.mcool" --check-base-resolution --compression-lvl 1

if ! "$hictk_bin_opt" validate --validate-index "$outdir/out.mcool"; then
  status=1
fi

if [ "$status" -eq 0 ]; then
  printf '\n### PASS ###\n'
else
  cat "$outdir/out.txt"
  printf '\n### FAIL ###\n'
fi

exit "$status"
