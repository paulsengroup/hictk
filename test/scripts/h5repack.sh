#!/usr/bin/env bash

# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -u

tmpdir="$(mktemp -d)"

# shellcheck disable=SC2064
trap "rm -rf '$tmpdir'" EXIT

for f in "$@"; do
  dest="$tmpdir/$(basename "$f")"
  h5repack -f SHUF -f GZIP=9 -L "$f" "$dest"

  old_size="$(stat -c%s "$f")"
  new_size="$(stat -c%s "$dest")"

  if [ "$old_size" -gt "$new_size" ]; then
    2>&1 echo "$f: space saved $((old_size - new_size))"
    mv "$dest" "$f"
  else
    2>&1 echo "$f: size increase after repack: leaving original file untouched..."
    rm "$dest"
  fi
done
