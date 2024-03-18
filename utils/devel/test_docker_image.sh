#!/usr/bin/env bash

# Copyright (c) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -eu
set -o pipefail

if [ $# -ne 1 ]; then
  2>&1 echo "Usage: $0 hictk:latest"
  exit 1
fi

IMG="$1"

tmpdir="$(mktemp -d)"
trap "rm -rf '$tmpdir'" EXIT

cat > "$tmpdir/runme.sh" <<- 'EOM'

set -eu
cd /tmp/hictk

apt-get update
apt-get install -y --no-install-recommends \
  curl \
  python3 \
  xz-utils

tar -xf test/data/hictk_test_data.tar.xz

tmpdir="$(mktemp -d)"
trap "rm -rf '$tmpdir'" EXIT

whereis -b hictk

test/scripts/hictk_convert_hic2cool.sh "$(which hictk)"
test/scripts/hictk_convert_cool2hic.sh "$(which hictk)"

EOM

chmod 755 "$tmpdir/runme.sh"

sudo docker run --rm --entrypoint=/bin/bash \
  -v "$tmpdir/runme.sh:/tmp/runme.sh:ro" \
  -v "$PWD/test/scripts:/tmp/hictk/test/scripts:ro" \
  -v "$PWD/test/data/hictk_test_data.tar.xz:/tmp/hictk/test/data/hictk_test_data.tar.xz:ro" \
  "$IMG" \
  /tmp/runme.sh
