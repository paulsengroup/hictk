#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u


if [ $# -ne 3 ]; then
  2>&1 echo "Usage:   $0 path_to_juicer_tools_jar path_to_hic_tools_jar output_prefix"
  2>&1 echo "Example: $0 juicer_tools_1.22.01.jar hic_tools.3.30.00.jar test_file"
  exit 1
fi

juicer_tools_jar="$1"
hic_tools_jar="$2"
out_prefix="$3"

if ! command -v java &> /dev/null; then
  2>&1 echo "Unable to find java in your PATH"
fi

if ! java -jar "$juicer_tools_jar" &> /dev/null; then
  2>&1 echo "Unable to run $juicer_tools_jar in your PATH"
  exit 1
fi

if ! java -jar "$hic_tools_jar" &> /dev/null; then
  2>&1 echo "Unable to run $hic_tools_jar in your PATH"
  exit 1
fi

if ! command -v cooler &> /dev/null; then
  2>&1 echo "Unable to find cooler in your PATH"
  exit 1
fi

if ! command -v curl &> /dev/null; then
  2>&1 echo "Unable to find curl in your PATH"
  exit 1
fi

if ! command -v pigz &> /dev/null; then
  2>&1 echo "Unable to find pigz in your PATH"
  exit 1
fi

touch "$out_prefix.hic8"
touch "$out_prefix.hic9"

wd="$(mktemp -d)"
trap 'rm -rf -- "$wd"' EXIT

url='https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/1cf3518f-839a-42b9-b2c7-7f81ad5935c3/4DNFIZ1ZVXC8.mcool'
src="$wd/$(basename "$url")"
assembly='dm6'

curl -L "$url" -o "$src"

cooler dump -t pixels --join "$src::/resolutions/5000" |
  awk -F '\t' 'BEGIN{ OFS=FS } {print "1",$1,int(($2+$3)/2),"0","0",$4,int(($5+$6)/2),"1",$7}' |
  sort -k2,2 -k6,6 --parallel=8 -S 8G |
  pigz -9 > "$wd/pixels.txt.gz"


java -Xms512m -Xmx8g -jar "$juicer_tools_jar" \
  pre \
  -k KR,VC,VC_SQRT,SCALE \
  "$wd/pixels.txt.gz" "$out_prefix.hic8" "$assembly"

java -Xms512m -Xmx8g -jar "$hic_tools_jar" \
  pre \
  -k KR,VC,VC_SQRT,SCALE \
  "$wd/pixels.txt.gz" "$out_prefix.hic9" "$assembly"
