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

# Install b2 using Conan
RUN printf '[requires]\nb2/5.2.1\n[options]\nb2*:toolset=%s' \
           "$(basename "$(which "$CC")")" | cut -f 1 -d - > /tmp/conanfile.txt

RUN conan install /tmp/conanfile.txt                 \
                 --build=missing                     \
                 -pr:b="$CONAN_DEFAULT_PROFILE_PATH" \
                 -pr:h="$CONAN_DEFAULT_PROFILE_PATH"

# Build hictk deps using Conan
RUN mkdir -p "$src_dir"

COPY conanfile.py "$src_dir"
RUN conan install "$src_dir/conanfile.py"       \
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
COPY test "$src_dir/test/"

ARG GIT_HASH
ARG GIT_SHORT_HASH
ARG GIT_TAG
ARG GIT_IS_DIRTY

RUN if [ -z "$GIT_HASH" ]; then echo "Missing GIT_HASH --build-arg" && exit 1; fi \
&&  if [ -z "$GIT_SHORT_HASH" ]; then echo "Missing GIT_SHORT_HASH --build-arg" && exit 1; fi \
&&  if [ -z "$GIT_IS_DIRTY" ]; then echo "Missing GIT_IS_DIRTY --build-arg" && exit 1; fi \
&&  if [ -z "$GIT_TAG" ]; then echo "Missing GIT_TAG --build-arg" && exit 1; fi

# Configure project
RUN cmake -DCMAKE_BUILD_TYPE=Release            \
          -DCMAKE_PREFIX_PATH="$build_dir"      \
          -DENABLE_DEVELOPER_MODE=OFF           \
          -DCMAKE_INSTALL_PREFIX="$staging_dir" \
          -DGIT_RETRIEVED_STATE=true            \
          -DGIT_TAG="$GIT_TAG"                  \
          -DGIT_IS_DIRTY="$GIT_IS_DIRTY"        \
          -DGIT_HEAD_SHA1="$GIT_HASH"           \
          -DGIT_DESCRIBE="$GIT_SHORT_HASH"      \
          -G Ninja                              \
          -S "$src_dir"                         \
          -B "$build_dir"

# Build and install project
RUN cmake --build "$build_dir" -t hictk -j "$(nproc)"  \
&& cmake --install "$build_dir" \
&& rm -rf "$build_dir/include" "$build_dir/lib"

ARG FINAL_BASE_IMAGE
ARG FINAL_BASE_IMAGE_DIGEST
FROM "${FINAL_BASE_IMAGE}@${FINAL_BASE_IMAGE_DIGEST}" AS base

ARG staging_dir='/root/hictk/staging'
ARG install_dir='/usr/local'

ARG BUILD_BASE_IMAGE
ARG FINAL_BASE_IMAGE
ARG FINAL_BASE_IMAGE_DIGEST

ARG GIT_HASH
ARG GIT_SHORT_HASH
ARG VERSION
ARG CREATION_DATE

RUN if [ -z "$BUILD_BASE_IMAGE" ]; then echo "Missing BUILD_BASE_IMAGE --build-arg" && exit 1; fi \
&&  if [ -z "$FINAL_BASE_IMAGE" ]; then echo "Missing FINAL_BASE_IMAGE --build-arg" && exit 1; fi \
&&  if [ -z "$FINAL_BASE_IMAGE_DIGEST" ]; then echo "Missing FINAL_BASE_IMAGE_DIGEST --build-arg" && exit 1; fi \
&&  if [ -z "$GIT_HASH" ]; then echo "Missing GIT_HASH --build-arg" && exit 1; fi \
&&  if [ -z "$GIT_SHORT_HASH" ]; then echo "Missing GIT_SHORT_HASH --build-arg" && exit 1; fi \
&&  if [ -z "$CREATION_DATE" ]; then echo "Missing CREATION_DATE --build-arg" && exit 1; fi

# Export project binaries to the final build stage
COPY --from=builder "$staging_dir" "$install_dir"

WORKDIR /data
ENTRYPOINT ["/usr/local/bin/hictk"]

RUN hictk --help
RUN hictk --version

# https://github.com/opencontainers/image-spec/blob/main/annotations.md#pre-defined-annotation-keys
LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/hictk'
LABEL org.opencontainers.image.documentation='https://github.com/paulsengroup/hictk'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/hictk'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title='hictk'
LABEL org.opencontainers.image.description='CLI utility to convert .hic files to .cool and .mcool'
LABEL org.opencontainers.image.base.digest="$FINAL_BASE_IMAGE_DIGEST"
LABEL org.opencontainers.image.base.name="$FINAL_BASE_IMAGE"
LABEL paulsengroup.hictk.image.build-base="$BUILD_BASE_IMAGE"

LABEL org.opencontainers.image.revision="$GIT_HASH"
LABEL org.opencontainers.image.created="$CREATION_DATE"
LABEL org.opencontainers.image.version="${VERSION:-sha-$GIT_SHORT_HASH}"
