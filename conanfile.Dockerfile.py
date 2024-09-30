# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

from conan import ConanFile
from conan.tools.build import check_min_cppstd

required_conan_version = ">=1.53.0"


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
        self.requires("bshoshany-thread-pool/4.1.0#be1802a8768416a6c9b1393cf0ce5e9c")
        self.requires("bzip2/1.0.8#d00dac990f08d991998d624be81a9526")
        self.requires("cli11/2.4.2#1b431bda2fb2cd3efed633899abcd8cc")
        self.requires("concurrentqueue/1.0.4#1e48e1c712bcfd892087c9c622a51502")
        self.requires("fast_float/6.1.1#e29acaa3d0543dee343abe3f6815346e")
        self.requires("fmt/10.2.1#9199a7a0611866dea5c8849a77467b25")
        self.requires("hdf5/1.14.3#31ccd8d4de83844f5db48471df1944a1")
        self.requires("highfive/2.9.0#c57477beed8b0110fadeb6da8f48bcc5")
        self.requires("libarchive/3.7.4#db39e5a5cddfd3a1884daaab4b7942fb")
        self.requires("libdeflate/1.20#3bd86e0160becf992346aa68b4938c40")
        self.requires("lz4/1.10.0#68a01ece147a441b463d8cefea68d555", force=True)
        self.requires("lzo/2.10#5725914235423c771cb1c6b607109b45")
        self.requires("nlohmann_json/3.11.3#45828be26eb619a2e04ca517bb7b828d")
        self.requires("parallel-hashmap/1.3.12#dc7755096d8a1fac7792fdd85760b6ca")
        self.requires("readerwriterqueue/1.0.6#aaa5ff6fac60c2aee591e9e51b063b83")
        self.requires("span-lite/0.11.0#519fd49fff711674cfed8cd17d4ed422")
        self.requires("spdlog/1.14.1#972bbf70be1da4bc57ea589af0efde03")
        self.requires("tomlplusplus/3.4.0#85dbfed71376fb8dc23cdcc0570e4727")
        self.requires("xz_utils/5.4.5#b885d1d79c9d30cff3803f7f551dbe66")
        self.requires("zstd/1.5.6#afefe79a309bc2a7b9f56c2093504c8b", force=True)
        self.requires("zlib/1.3.1#f52e03ae3d251dec704634230cd806a2")

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