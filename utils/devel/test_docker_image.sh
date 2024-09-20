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

whereis -b hictk

apt-get update
apt-get install -q -y --no-install-recommends \
  python3 \
  python3-pip \
  python3-venv \
  tar \
  which \
  zstd

cd /tmp/hictk

tar -xf test/data/hictk_test_data.tar.zst

python3 -m venv venv --upgrade
venv/bin/pip install test/integration

venv/bin/hictk_integration_suite \
  "$(which hictk)" \
  test/integration/config.toml \
  --data-dir test/data \
  --threads "$(nproc)"

EOM

chmod 755 "$tmpdir/runme.sh"

sudo docker run --rm --entrypoint=/bin/bash \
  -v "$tmpdir/runme.sh:/tmp/runme.sh:ro" \
  -v "$PWD/test/integration:/tmp/hictk/test/integration:ro" \
  -v "$PWD/test/data/hictk_test_data.tar.zst:/tmp/hictk/test/data/hictk_test_data.tar.zst:ro" \
  "$IMG" \
  /tmp/runme.sh
