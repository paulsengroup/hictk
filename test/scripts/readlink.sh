#!/usr/bin/env bash

# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

if [ $# -ne 1 ]; then
  2>&1 echo "Usage: $0 file"
  exit 1
fi

python3 -c 'import os, sys; print(os.path.realpath(sys.argv[1]))' "$1"
