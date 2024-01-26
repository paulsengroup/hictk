#!/usr/bin/env bash

# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

if [ $# -ne 2 ]; then
  2>&1 echo "Usage: $0 file1 file2"
  exit 1
fi

set -o pipefail
set -e

2>&1 echo "Comparing $1 with $2..."
if diff "$1" "$2"; then
  2>&1 echo "Files are identical"
  exit 0
else
  2>&1 echo "Files differ"
  exit 1
fi
