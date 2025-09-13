#!/usr/bin/env bash

# Copyright (c) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -u
set -o pipefail

IMAGE_NAME='hictk'

if [ "$(uname)" == "Darwin" ]; then
  BUILD_USER="$USER"
else
  BUILD_USER='root'
fi

COMPILER_NAME='clang'
COMPILER_VERSION='21'
C_COMPILER="$COMPILER_NAME-$COMPILER_VERSION"
if [[ "$COMPILER_NAME" == gcc ]]; then
  CXX_COMPILER="g++-$COMPILER_VERSION"
else
  CXX_COMPILER="clang++-$COMPILER_VERSION"
fi

GIT_HASH="$(git rev-parse HEAD)"
GIT_SHORT_HASH="$(git rev-parse --short HEAD)"
GIT_TAG="$(git for-each-ref 'refs/tags/v*.*.*' --count 1 --sort=-v:refname --format "%(refname:short)"  --points-at HEAD)"
CREATION_DATE="$(date -I)"

if [[ $(git status --porcelain -uno) ]]; then
  GIT_IS_DIRTY=1
else
  GIT_IS_DIRTY=0
fi

IMAGE_TAG="sha-$GIT_SHORT_HASH"
if [ $GIT_IS_DIRTY -ne 0 ]; then
  IMAGE_TAG+='-dirty'
fi

if [ -z "$GIT_TAG" ]; then
  GIT_TAG="sha-$GIT_SHORT_HASH"
fi

1>&2 echo "Building \"$IMAGE_NAME:$IMAGE_TAG\"..."

sudo -u "$BUILD_USER" docker pull docker.io/library/ubuntu:24.04
FINAL_BASE_IMAGE_DIGEST="$(sudo -u "$BUILD_USER" docker inspect --format='{{index .RepoDigests 0}}' docker.io/library/ubuntu:24.04 | grep -o '[[:alnum:]:]\+$')"

BUILD_BASE_IMAGE="ghcr.io/paulsengroup/ci-docker-images/ubuntu-24.04-cxx-$C_COMPILER:latest"

sudo -u "$BUILD_USER" docker pull "$BUILD_BASE_IMAGE"

# sudo -u "$BUILD_USER" docker buildx build --platform linux/amd64,linux/arm64 \
sudo -u "$BUILD_USER" docker buildx build --platform linux/amd64 \
  --build-arg "BUILD_BASE_IMAGE=$BUILD_BASE_IMAGE" \
  --build-arg "FINAL_BASE_IMAGE=docker.io/library/ubuntu" \
  --build-arg "FINAL_BASE_IMAGE_TAG=24.04" \
  --build-arg "FINAL_BASE_IMAGE_DIGEST=$FINAL_BASE_IMAGE_DIGEST" \
  --build-arg "C_COMPILER=$C_COMPILER" \
  --build-arg "CXX_COMPILER=$CXX_COMPILER" \
  --build-arg "HICTK_GIT_HASH=$GIT_HASH" \
  --build-arg "HICTK_GIT_SHORT_HASH=$GIT_SHORT_HASH" \
  --build-arg "HICTK_GIT_TAG=$GIT_TAG" \
  --build-arg "HICTK_GIT_IS_DIRTY=$GIT_IS_DIRTY" \
  --build-arg "CREATION_DATE=$CREATION_DATE" \
  -t "$IMAGE_NAME:latest" \
  -t "$IMAGE_NAME:$(echo "$CREATION_DATE" | tr -d '\-' )" \
  -t "$IMAGE_NAME:$IMAGE_TAG" \
  "$(git rev-parse --show-toplevel)"

 # sudo singularity build -F "${img_name}_v${ver}.sif" \
 #                           "docker-daemon://${img_name}:${ver}"
