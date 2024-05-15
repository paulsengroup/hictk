# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

from conan import ConanFile
from conan.tools.build import check_min_cppstd


required_conan_version = ">=1.53.0"


class NCHGConan(ConanFile):
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
        self.requires("catch2/3.5.4#d346ca291f8f62040fd9c1a891654711")
        self.requires("cli11/2.4.1#afacffd31f631bbb8b7c7d6425fe7a66")
        self.requires("concurrentqueue/1.0.4#1e48e1c712bcfd892087c9c622a51502")
        self.requires("eigen/3.4.0#2e192482a8acff96fe34766adca2b24c")
        self.requires("fast_float/6.1.1#e29acaa3d0543dee343abe3f6815346e")
        self.requires("fmt/10.2.1#9199a7a0611866dea5c8849a77467b25")
        self.requires("hdf5/1.14.3#31ccd8d4de83844f5db48471df1944a1")
        self.requires("highfive/2.9.0#c57477beed8b0110fadeb6da8f48bcc5")
        self.requires("libdeflate/1.19#3ea74a4549efc14d4b1202dc4bfbf602")
        self.requires("parallel-hashmap/1.3.11#1e67f4855a3f7cdeb977cc472113baf7")
        self.requires("readerwriterqueue/1.0.6#aaa5ff6fac60c2aee591e9e51b063b83")
        self.requires("span-lite/0.11.0#519fd49fff711674cfed8cd17d4ed422")
        self.requires("spdlog/1.13.0#8e88198fd5b9ee31d329431a6d0ccaa2")
        self.requires("zstd/1.5.6#67383dae85d33f43823e7751a6745ea1")

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
        self.options["hictk"].with_eigen = False
        self.options["highfive"].with_boost = False
        self.options["highfive"].with_eigen = False
        self.options["highfive"].with_opencv = False
        self.options["highfive"].with_xtensor = False
        self.options["spdlog"].header_only = True
        self.options["zstd"].build_programs = False
