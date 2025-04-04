#!/usr/bin/env bash

# Copyright (c) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -eu
set -o pipefail

if [ $# -ne 1 ]; then
  1>&2 echo "Usage: $0 hictk:latest"
  exit 1
fi

IMG="$1"

tmpdir="$(mktemp -d)"
# shellcheck disable=SC2064
trap "rm -rf '$tmpdir'" EXIT

cat > "$tmpdir/runme.sh" <<- 'EOM'

set -eu

whereis -b hictk

apt-get update

python="$(apt-cache search '^python3\.[0-9]+$' | cut -f 1 -d ' ' | sort -V | tail -n 1)"

apt-get install -q -y --no-install-recommends \
  "$python" \
  python3-pip \
  "$python"-venv \
  tar \
  zstd

cd /tmp/hictk

tar -xf test/data/hictk_test_data.tar.zst

"$python" -m venv venv --upgrade
venv/bin/pip install test/integration

venv/bin/hictk_integration_suite \
  "$(which hictk)" \
  test/integration/config.toml \
  --data-dir test/data \
  --threads "$(nproc)"

EOM

chmod 755 "$tmpdir/runme.sh"

sudo docker run --rm --entrypoint=/bin/bash \
  --volume "$tmpdir/runme.sh:/tmp/runme.sh:ro" \
  --volume "$PWD/test/integration:/tmp/hictk/test/integration:ro" \
  --volume "$PWD/test/data/hictk_test_data.tar.zst:/tmp/hictk/test/data/hictk_test_data.tar.zst:ro" \
  --env 'HICTK_NO_TELEMETRY=1' \
  "$IMG" \
  /tmp/runme.sh
