# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

##### IMPORTANT #####
# This Dockerfile requires several build arguments to be defined through --build-arg
# See utils/devel/build_dockerfile.sh for an example of how to build this Dockerfile
#####################

ARG BUILD_BASE_IMAGE=ghcr.io/paulsengroup/ci-docker-images/ubuntu-24.04-cxx-clang-20:latest
ARG FINAL_BASE_IMAGE=ubuntu:24.04

FROM "${BUILD_BASE_IMAGE}" AS builder

ARG src_dir='/root/hictk'
ARG build_dir='/root/hictk/build'
ARG staging_dir='/root/hictk/staging'
ARG install_dir='/usr/local'

# Environment variables used for building
ARG CCACHE_DISABLE=1
ARG CMAKE_BUILD_TYPE=Release
ARG CMAKE_GENERATOR=Ninja
ARG CMAKE_INSTALL_PREFIX="$staging_dir"
ARG CMAKE_POLICY_VERSION_MINIMUM=3.5
ARG CMAKE_PREFIX_PATH="$build_dir"
ARG CXX_STANDARD=23

ARG HICTK_GIT_HASH=0000000000000000000000000000000000000000
ARG HICTK_GIT_IS_DIRTY=false
ARG HICTK_GIT_SHORT_HASH=00000000
ARG HICTK_GIT_TAG=unknown

RUN mkdir -p "$src_dir" "$build_dir"

# Build hictk deps using Conan
COPY conanfile.py "$src_dir/conanfile.py"
RUN conan install "$src_dir/conanfile.py"                      \
             --build=missing                                   \
             -c:a=tools.cmake.cmaketoolchain:generator=Ninja   \
             --options='hictk/*:with_arrow=False'              \
             --options='hictk/*:with_benchmark_deps=False'     \
             --options='hictk/*:with_cli_tool_deps=True'       \
             --options='hictk/*:with_eigen=False'              \
             --options='hictk/*:with_fuzzy_testing_deps=False' \
             --options='hictk/*:with_telemetry_deps=True'      \
             --options='hictk/*:with_unit_testing_deps=False'  \
             --output-folder="$build_dir"                      \
             -pr:b="$CONAN_DEFAULT_PROFILE_PATH"               \
             -pr:h="$CONAN_DEFAULT_PROFILE_PATH"               \
             --settings=build_type=Release                     \
             --settings=compiler.cppstd="$CXX_STANDARD"        \
&&  conan cache clean "*" --build                              \
&&  conan cache clean "*" --download                           \
&&  conan cache clean "*" --source

# Copy source files
COPY LICENSE "$src_dir/"
COPY external "$src_dir/external/"
COPY cmake "$src_dir/cmake/"
COPY CMakeLists.txt "$src_dir/"
COPY src "$src_dir/src/"

# Configure project
RUN cmake -DCMAKE_CXX_STANDARD="$CXX_STANDARD"                \
          -DCMAKE_LINKER_TYPE=LLD                             \
          -DENABLE_DEVELOPER_MODE=OFF                         \
          -DHICTK_ENABLE_GIT_VERSION_TRACKING=OFF             \
          -DHICTK_ENABLE_TESTING=OFF                          \
          -DHICTK_GIT_DESCRIBE="$HICTK_GIT_SHORT_HASH"        \
          -DHICTK_GIT_HEAD_SHA1="$HICTK_GIT_HASH"             \
          -DHICTK_GIT_IS_DIRTY="$HICTK_GIT_IS_DIRTY"          \
          -DHICTK_GIT_RETRIEVED_STATE=true                    \
          -DHICTK_GIT_TAG="$HICTK_GIT_TAG"                    \
          -DHICTK_WITH_ARROW=OFF                              \
          -DHICTK_WITH_EIGEN=OFF                              \
          -S "$src_dir"                                       \
          -B "$build_dir"

# Build and install project
RUN cmake --build "$build_dir" -t hictk -j "$(nproc)"  \
&& cmake --install "$build_dir" --component Runtime    \
&& rm -r "$build_dir"

ARG FINAL_BASE_IMAGE=ubuntu:24.04

FROM "${FINAL_BASE_IMAGE}" AS base

ARG staging_dir='/root/hictk/staging'
ARG install_dir='/usr/local'

ARG BUILD_BASE_IMAGE=ghcr.io/paulsengroup/ci-docker-images/ubuntu-24.04-cxx-clang-20:latest
ARG FINAL_BASE_IMAGE=ubuntu

ARG HICTK_GIT_HASH=0000000000000000000000000000000000000000
ARG HICTK_GIT_SHORT_HASH=00000000
ARG HICTK_GIT_TAG=unknown
ARG HICTK_GIT_IS_DIRTY=false
ARG VERSION=latest
ARG CREATION_DATE=2000-01-01

RUN if [ "$BUILDARCH" != 'amd64' ]; then \
    apt-get update \
&&  apt-get install -q -y --no-install-recommends libatomic1 \
&&  rm -rf /var/lib/apt/lists/*; \
fi

RUN apt-get update \
&&  apt-get install -q -y --no-install-recommends ca-certificates \
&&  rm -rf /var/lib/apt/lists/*

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
LABEL org.opencontainers.image.base.name="$FINAL_BASE_IMAGE"
LABEL paulsengroup.hictk.image.build-base="$BUILD_BASE_IMAGE"
LABEL org.opencontainers.image.revision="$HICTK_GIT_HASH"
LABEL org.opencontainers.image.created="$CREATION_DATE"
LABEL org.opencontainers.image.version="${VERSION:-sha-$HICTK_GIT_SHORT_HASH}"
