# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

##### IMPORTANT #####
# This Dockerfile requires several build arguments to be defined through --build-arg
# See utils/devel/build_dockerfile.sh for an example of how to build this Dockerfile
#####################

ARG BUILD_BASE_IMAGE
ARG FINAL_BASE_IMAGE
ARG FINAL_BASE_IMAGE_DIGEST

FROM "$BUILD_BASE_IMAGE" AS builder

ARG src_dir='/root/hictk'
ARG build_dir='/root/hictk/build'
ARG staging_dir='/root/hictk/staging'
ARG install_dir='/usr/local'


ARG C_COMPILER
ARG CXX_COMPILER

RUN if [ -z "$C_COMPILER" ]; then echo "Missing C_COMPILER --build-arg" && exit 1; fi \
&&  if [ -z "$CXX_COMPILER" ]; then echo "Missing CXX_COMPILER --build-arg" && exit 1; fi

ENV CC="$C_COMPILER"
ENV CXX="$CXX_COMPILER"
ENV CMAKE_POLICY_VERSION_MINIMUM=3.5

# Build hictk deps using Conan
RUN mkdir -p "$src_dir"

COPY conanfile.Dockerfile.py "$src_dir/conanfile.py"
RUN conan install "$src_dir/conanfile.py"        \
             --build=missing                     \
             -pr:b="$CONAN_DEFAULT_PROFILE_PATH" \
             -pr:h="$CONAN_DEFAULT_PROFILE_PATH" \
             -s build_type=Release               \
             -s compiler.libcxx=libstdc++11      \
             -s compiler.cppstd=17               \
             --output-folder="$build_dir"        \
&& conan cache clean "*" --build                 \
&& conan cache clean "*" --download              \
&& conan cache clean "*" --source

# Copy source files
COPY LICENSE "$src_dir/"
COPY external "$src_dir/external/"
COPY cmake "$src_dir/cmake/"
COPY CMakeLists.txt "$src_dir/"
COPY src "$src_dir/src/"

ARG HICTK_GIT_HASH
ARG HICTK_GIT_SHORT_HASH
ARG HICTK_GIT_TAG
ARG HICTK_GIT_IS_DIRTY

RUN if [ -z "$HICTK_GIT_HASH" ]; then echo "Missing HICTK_GIT_HASH --build-arg" && exit 1; fi \
&&  if [ -z "$HICTK_GIT_SHORT_HASH" ]; then echo "Missing HICTK_GIT_SHORT_HASH --build-arg" && exit 1; fi \
&&  if [ -z "$HICTK_GIT_IS_DIRTY" ]; then echo "Missing HICTK_GIT_IS_DIRTY --build-arg" && exit 1; fi \
&&  if [ -z "$HICTK_GIT_TAG" ]; then echo "Missing HICTK_GIT_TAG --build-arg" && exit 1; fi

ARG CCACHE_DISABLE=1

# Configure project
RUN cmake -DCMAKE_BUILD_TYPE=Release                   \
          -DCMAKE_PREFIX_PATH="$build_dir"             \
          -DENABLE_DEVELOPER_MODE=OFF                  \
          -DHICTK_ENABLE_TESTING=OFF                   \
          -DHICTK_WITH_ARROW=OFF                       \
          -DHICTK_WITH_EIGEN=OFF                       \
          -DCMAKE_INSTALL_PREFIX="$staging_dir"        \
          -DHICTK_GIT_RETRIEVED_STATE=true             \
          -DHICTK_GIT_TAG="$HICTK_GIT_TAG"             \
          -DHICTK_GIT_IS_DIRTY="$HICTK_GIT_IS_DIRTY"   \
          -DHICTK_GIT_HEAD_SHA1="$HICTK_GIT_HASH"      \
          -DHICTK_GIT_DESCRIBE="$HICTK_GIT_SHORT_HASH" \
          -G Ninja                                     \
          -S "$src_dir"                                \
          -B "$build_dir"

# Build and install project
RUN cmake --build "$build_dir" -t hictk -j "$(nproc)"  \
&& cmake --install "$build_dir" --component Runtime    \
&& rm -r "$build_dir"

ARG FINAL_BASE_IMAGE
ARG FINAL_BASE_IMAGE_DIGEST
FROM "${FINAL_BASE_IMAGE}@${FINAL_BASE_IMAGE_DIGEST}" AS base

ARG staging_dir='/root/hictk/staging'
ARG install_dir='/usr/local'

ARG BUILD_BASE_IMAGE
ARG FINAL_BASE_IMAGE
ARG FINAL_BASE_IMAGE_DIGEST

ARG HICTK_GIT_HASH
ARG HICTK_GIT_SHORT_HASH
ARG VERSION
ARG CREATION_DATE

RUN if [ -z "$BUILD_BASE_IMAGE" ]; then echo "Missing BUILD_BASE_IMAGE --build-arg" && exit 1; fi \
&&  if [ -z "$FINAL_BASE_IMAGE" ]; then echo "Missing FINAL_BASE_IMAGE --build-arg" && exit 1; fi \
&&  if [ -z "$FINAL_BASE_IMAGE_DIGEST" ]; then echo "Missing FINAL_BASE_IMAGE_DIGEST --build-arg" && exit 1; fi \
&&  if [ -z "$HICTK_GIT_HASH" ]; then echo "Missing HICTK_GIT_HASH --build-arg" && exit 1; fi \
&&  if [ -z "$HICTK_GIT_SHORT_HASH" ]; then echo "Missing HICTK_GIT_SHORT_HASH --build-arg" && exit 1; fi \
&&  if [ -z "$CREATION_DATE" ]; then echo "Missing CREATION_DATE --build-arg" && exit 1; fi

RUN if [ "$BUILDARCH" != 'amd64' ]; then \
    apt-get update \
&&  apt-get install -q -y --no-install-recommends libatomic1 \
&&  rm -rf /var/lib/apt/lists/*; \
fi

# Export project binaries to the final build stage
COPY --from=builder "$staging_dir" "$install_dir"

WORKDIR /data
ENTRYPOINT ["/usr/local/bin/hictk"]

RUN hictk --help
RUN hictk --version

# https://github.com/opencontainers/image-spec/blob/main/annotations.md#pre-defined-annotation-keys
LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/hictk'
LABEL org.opencontainers.image.documentation='https://hictk.readthedocs.io/en/stable/'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/hictk'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title='hictk'
LABEL org.opencontainers.image.description='Blazing fast toolkit to work with .hic and .cool files'
LABEL org.opencontainers.image.base.digest="$FINAL_BASE_IMAGE_DIGEST"
LABEL org.opencontainers.image.base.name="$FINAL_BASE_IMAGE"
LABEL paulsengroup.hictk.image.build-base="$BUILD_BASE_IMAGE"

LABEL org.opencontainers.image.revision="$HICTK_GIT_HASH"
LABEL org.opencontainers.image.created="$CREATION_DATE"
LABEL org.opencontainers.image.version="${VERSION:-sha-$HICTK_GIT_SHORT_HASH}"
