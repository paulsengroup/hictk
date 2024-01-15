#!/usr/bin/env bash

# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

if [ $# -lt 3 ]; then
  2>&1 echo "Usage: $0 path_to_hictk file1 file2 [resolution]"
  exit 1
fi

resolution="${4-0}"
f1="$2"
f2="$3"


function dump_table {
  set -o pipefail
  set -eu

  hictk="$1"
  f="$2"
  table="$3"
  resolution="$4"

  args=()

  if [ "$table" = pixels ]; then
    args+=(--join)
  fi

  if [ ! "$resolution" = 0 ]; then
    args+=(--resolution "$resolution")
  fi

  "$hictk" dump "${args[@]}" --table "$table" "$f"
}


function compare_chromosomes {
  set -o pipefail
  set -eu

  hictk="$1"
  f1="$2"
  f2="$3"

  diff <(dump_table "$hictk" "$f1" chroms 0) \
       <(dump_table "$hictk" "$f2" chroms 0);
}

function compare_bins {
  set -o pipefail
  set -eu

  hictk="$1"
  f1="$2"
  f2="$3"
  resolution="$4"

  diff <(dump_table "$hictk" "$f1" bins "$resolution") \
       <(dump_table "$hictk" "$f2" bins "$resolution");
}

function compare_pixels {
  set -o pipefail
  set -eu

  hictk="$1"
  f1="$2"
  f2="$3"
  resolution="$4"

  diff <(dump_table "$hictk" "$f1" pixels "$resolution") \
       <(dump_table "$hictk" "$f2" pixels "$resolution");
}


status=0
2>&1 echo "Comparing $f1 with $f2..."
if ! compare_chromosomes "$hictk" "$f1" "$f2" "$resolution"; then
  status=1
fi

if ! compare_bins "$hictk" "$f1" "$f2" "$resolution"; then
  status=1
fi

if ! compare_pixels "$hictk" "$f1" "$f2" "$resolution"; then
  status=1
fi

if [ "$status" -eq 0 ]; then
  2>&1 echo "Files are identical"
  exit 0
else
  2>&1 echo "Files differ"
  exit 1
fi
