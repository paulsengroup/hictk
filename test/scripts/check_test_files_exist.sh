#!/usr/bin/env bash

# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

if [ $# -lt 1 ]; then
  2>&1 echo "Usage: $0 files..."
  exit 1
fi


status=0
for f in "$@"; do
  if [ ! -f "$f" ]; then
    2>&1 echo "Unable to find test file \"$f\""
    status=1
  fi
done

exit "$status"
