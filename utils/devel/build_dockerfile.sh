#!/usr/bin/env bash

# Copyright (c) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -u
set -o pipefail

IMAGE_NAME='hictk'

OS_NAME=ubuntu
OS_VERSION=24.04

COMPILER_NAME='clang'
COMPILER_VERSION='20'
C_COMPILER="$COMPILER_NAME-$COMPILER_VERSION"

FINAL_BASE_IMAGE_TAG="$OS_VERSION"
FINAL_BASE_IMAGE="docker.io/library/$OS_NAME:$FINAL_BASE_IMAGE_TAG"

if [ "$(uname)" == "Darwin" ]; then
  BUILD_USER="$USER"
else
  BUILD_USER='root'
fi

sudo -u "$BUILD_USER" docker pull "$FINAL_BASE_IMAGE"
FINAL_BASE_IMAGE_DIGEST="$(sudo -u "$BUILD_USER" docker inspect --format='{{index .RepoDigests 0}}' "$FINAL_BASE_IMAGE" | grep -o '[[:alnum:]:]\+$')"

BUILD_BASE_IMAGE="ghcr.io/paulsengroup/ci-docker-images/$OS_NAME-$OS_VERSION-cxx-$C_COMPILER:latest"
FINAL_BASE_IMAGE="$(echo "${FINAL_BASE_IMAGE#docker.io/library/}" | cut -d : -f 1)"

HICTK_GIT_HASH="$(git rev-parse HEAD)"
HICTK_GIT_SHORT_HASH="$(git rev-parse --short HEAD)"
HICTK_GIT_TAG="$(git for-each-ref 'refs/tags/v*.*.*' --count 1 --sort=-v:refname --format "%(refname:short)" --points-at HEAD)"
CREATION_DATE="$(date -I)"

if [[ $(git status --porcelain -uno) ]]; then
  HICTK_GIT_IS_DIRTY=1
else
  HICTK_GIT_IS_DIRTY=0
fi

IMAGE_TAG="sha-$HICTK_GIT_SHORT_HASH"
if [ $HICTK_GIT_IS_DIRTY -ne 0 ]; then
  IMAGE_TAG+='-dirty'
fi

if [ -z "$HICTK_GIT_TAG" ]; then
  HICTK_GIT_TAG="sha-$HICTK_GIT_SHORT_HASH"
fi

1>&2 echo "Building \"$IMAGE_NAME:$IMAGE_TAG\"..."

sudo -u "$BUILD_USER" docker pull "$BUILD_BASE_IMAGE"

set -x

# sudo -u "$BUILD_USER" docker buildx build --platform linux/amd64,linux/arm64 \
sudo -u "$BUILD_USER" docker buildx build \
  --build-arg "BUILD_BASE_IMAGE=$BUILD_BASE_IMAGE" \
  --build-arg "FINAL_BASE_IMAGE=ubuntu@$FINAL_BASE_IMAGE_DIGEST" \
  --build-arg "FINAL_BASE_IMAGE_TAG=$FINAL_BASE_IMAGE_TAG" \
  --build-arg "HICTK_GIT_HASH=$HICTK_GIT_HASH" \
  --build-arg "HICTK_GIT_SHORT_HASH=$HICTK_GIT_SHORT_HASH" \
  --build-arg "HICTK_GIT_TAG=$HICTK_GIT_TAG" \
  --build-arg "HICTK_GIT_IS_DIRTY=$HICTK_GIT_IS_DIRTY" \
  --build-arg "CREATION_DATE=$CREATION_DATE" \
  -t "$IMAGE_NAME:latest" \
  -t "$IMAGE_NAME:$(echo "$CREATION_DATE" | tr -d '\-' )" \
  -t "$IMAGE_NAME:$IMAGE_TAG" \
  "$(git rev-parse --show-toplevel)" \
  --progress=plain --load

set +x

# sudo apptainer build -F "${img_name}_v${ver}.sif" \
#                         "docker-daemon://${img_name}:${ver}"
