# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

from conan import ConanFile
from conan.tools.build import check_min_cppstd

required_conan_version = ">=1.53.0"

# Debugging  TODO remove


class HictkConan(ConanFile):
    name = "hictk"
    description = "Blazing fast toolkit to work with .hic and .cool files."
    license = "MIT"
    topics = ("hictk", "bioinformatics")
    homepage = "https://github.com/paulsengroup/hictk"
    url = "https://github.com/paulsengroup/hictk"
    package_type = "header-library"
    settings = "os", "arch", "compiler", "build_type"

    options = {
        "shared": [True, False],
        "fPIC": [True, False],
    }

    default_options = {
        "shared": False,
        "fPIC": True,
    }

    generators = "CMakeDeps"

    @property
    def _min_cppstd(self):
        return 17

    def requirements(self):
        self.requires("bshoshany-thread-pool/5.0.0#d94da300363f0c35b8f41b2c5490c94d")
        self.requires("bzip2/1.0.8#d00dac990f08d991998d624be81a9526")
        self.requires("cli11/2.4.2#6e507e3d37c863d547be3069be114cce")
        self.requires("concurrentqueue/1.0.4#1e48e1c712bcfd892087c9c622a51502")
        self.requires("fast_float/7.0.0#e4a4a338590ab5eaaf517c64607629d0")
        self.requires("fmt/11.0.2#b4a24d70b93466b9b508ddb7b014643b")
        self.requires("hdf5/1.14.5#51799cda2ba7acaa74c9651dea284ac4", force=True)
        self.requires("highfive/2.10.0#3d1bd25944a57fa1bc30a0a22923d528")
        self.requires("libarchive/3.7.6#11e70b88b334f684eb9f6b65f287c81e")
        self.requires("libdeflate/1.23#4994bea7cf7e93789da161fac8e26a53")
        self.requires("lz4/1.10.0#68a01ece147a441b463d8cefea68d555", force=True)
        self.requires("lzo/2.10#5725914235423c771cb1c6b607109b45")
        self.requires("nlohmann_json/3.11.3#45828be26eb619a2e04ca517bb7b828d")
        self.requires("parallel-hashmap/2.0.0#82acae64ffe2693fff5fb3f9df8e1746")
        self.requires("readerwriterqueue/1.0.6#aaa5ff6fac60c2aee591e9e51b063b83")
        self.requires("span-lite/0.11.0#519fd49fff711674cfed8cd17d4ed422")
        self.requires("spdlog/1.15.0#da21f74dd84627fa68601c4e3b9c3f00")
        self.requires("tomlplusplus/3.4.0#85dbfed71376fb8dc23cdcc0570e4727")
        self.requires("xz_utils/5.4.5#b885d1d79c9d30cff3803f7f551dbe66")
        self.requires("zstd/1.5.6#afefe79a309bc2a7b9f56c2093504c8b", force=True)
        self.requires("zlib/1.3.1#b8bc2603263cf7eccbd6e17e66b0ed76")

    def validate(self):
        if self.settings.get_safe("compiler.cppstd"):
            check_min_cppstd(self, self._min_cppstd)

    def configure(self):
        if self.settings.compiler in ["clang", "gcc"]:
            self.settings.compiler.libcxx = "libstdc++11"

        self.options["fmt"].header_only = True
        self.options["hdf5"].enable_cxx = False
        self.options["hdf5"].hl = False
        self.options["hdf5"].threadsafe = False
        self.options["hdf5"].parallel = False
        self.options["highfive"].with_boost = False
        self.options["highfive"].with_eigen = False
        self.options["highfive"].with_opencv = False
        self.options["highfive"].with_xtensor = False
        self.options["libarchive"].with_acl = False
        self.options["libarchive"].with_zlib = True
        self.options["libarchive"].with_bzip2 = True
        self.options["libarchive"].with_libxml2 = False
        self.options["libarchive"].with_expat = False
        self.options["libarchive"].with_iconv = False
        self.options["libarchive"].with_acl = False
        self.options["libarchive"].with_pcreposix = False
        self.options["libarchive"].with_cng = False
        self.options["libarchive"].with_nettle = False
        self.options["libarchive"].with_openssl = False
        self.options["libarchive"].with_libb2 = False
        self.options["libarchive"].with_lz4 = True
        self.options["libarchive"].with_lzo = True
        self.options["libarchive"].with_lzma = True
        self.options["libarchive"].with_zstd = True
        self.options["libarchive"].with_mbedtls = False
        self.options["libarchive"].with_xattr = False
        self.options["libarchive"].with_pcre2 = False
        self.options["spdlog"].header_only = True
        self.options["zstd"].build_programs = False
