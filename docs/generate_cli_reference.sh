#!/usr/bin/env bash

# Copyright (c) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -eu
set -o pipefail

if [ $# -ne 1 ]; then
  2>&1 echo "Usage: $0 path/to/hictk"
  exit 1
fi

hictk="$1"

subcommands=(
  balance
  'balance ice'
  'balance scale'
  'balance vc'
  convert
  dump
  fix-mcool
  load
  merge
  metadata
  rename-chromosomes
  validate
  zoomify
)


cat << EOT
..
   Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
   SPDX-License-Identifier: MIT

CLI Reference
#############

For an up-to-date list of subcommands and CLI options refer to ``hictk --help``.

Subcommands
-----------

.. code-block:: text

EOT

"$hictk" --help |& sed "s|$hictk|hictk|g" | sed 's/^/  /' | sed '/^[[:space:]]*$/d'

for subcmd in "${subcommands[@]}"; do
  header="hictk $subcmd"
cat << EOT

$header
$(printf '\055%.0s' $(seq ${#header}))

.. code-block:: text

EOT

  # shellcheck disable=SC2086
  "$hictk" $subcmd --help |& sed "s|$hictk|hictk|g" | sed 's/^/  /' | sed '/^[[:space:]]*$/d'
done
