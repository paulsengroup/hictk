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
        self.requires("arrow/19.0.1#e5b40539e3f5411d03603a38b5fb698a")
        self.requires("boost/1.87.0#0c087f18c4e6487235dd10480613cbb5", force=True)
        self.requires("bshoshany-thread-pool/5.0.0#d94da300363f0c35b8f41b2c5490c94d")
        self.requires("bzip2/1.0.8#00b4a4658791c1f06914e087f0e792f5")
        self.requires("catch2/3.8.1#141f4cd552b86c7278436c434473ae2f")
        self.requires("cli11/2.5.0#1b7c81ea2bff6279eb2150bbe06a200a")
        self.requires("concurrentqueue/1.0.4#1e48e1c712bcfd892087c9c622a51502")
        self.requires("eigen/3.4.0#2e192482a8acff96fe34766adca2b24c")
        self.requires("fast_float/8.0.0#edda0315516b2f1e7835972fdf5fc5ca")
        self.requires("fmt/11.2.0#579bb2cdf4a7607621beea4eb4651e0f", force=True)
        self.requires("hdf5/1.14.5#51799cda2ba7acaa74c9651dea284ac4", force=True)
        self.requires("highfive/2.10.0#c975a16d7fe3655c173f8a9aab16b416")
        self.requires("libarchive/3.7.9#7a6b3ce684024d5ac74dc9049d5c61cf")
        self.requires("libcurl/8.12.1#5bb1e5168ab52aaee0df5d556e092f47", force=True)  # otel
        self.requires("libdeflate/1.23#4994bea7cf7e93789da161fac8e26a53")
        self.requires("lz4/1.10.0#59fc63cac7f10fbe8e05c7e62c2f3504", force=True)
        self.requires("lzo/2.10#5725914235423c771cb1c6b607109b45")
        self.requires("nlohmann_json/3.12.0#2d634ab0ec8d9f56353e5ccef6d6612c", force=True)
        self.requires("opentelemetry-cpp/1.18.0#edd81c9cc34028bf0c2da87fc9756e04")
        self.requires("parallel-hashmap/2.0.0#82acae64ffe2693fff5fb3f9df8e1746")
        self.requires("pybind11/2.13.6#7d1417680344884436657a0d12212274")
        self.requires("readerwriterqueue/1.0.6#aaa5ff6fac60c2aee591e9e51b063b83")
        self.requires("span-lite/0.11.0#519fd49fff711674cfed8cd17d4ed422")
        self.requires("spdlog/1.15.1#92e99f07f134481bce4b70c1a41060e7")
        self.requires("tomlplusplus/3.4.0#85dbfed71376fb8dc23cdcc0570e4727")
        self.requires("xz_utils/5.4.5#b885d1d79c9d30cff3803f7f551dbe66")
        self.requires("zstd/1.5.7#fde461c0d847a22f16d3066774f61b11", force=True)
        self.requires("zlib/1.3.1#b8bc2603263cf7eccbd6e17e66b0ed76")

    def validate(self):
        if self.settings.get_safe("compiler.cppstd"):
            check_min_cppstd(self, self._min_cppstd)

    def configure(self):
        if self.settings.compiler in ["clang", "gcc"] and self.settings.os == "Linux":
            self.settings.compiler.libcxx = "libstdc++11"

        self.options["arrow"].compute = True
        self.options["arrow"].filesystem_layer = False
        self.options["arrow"].parquet = False
        self.options["arrow"].with_boost = True
        self.options["arrow"].with_re2 = True
        self.options["arrow"].with_thrift = False
        self.options["boost"].system_no_deprecated = True
        self.options["boost"].asio_no_deprecated = True
        self.options["boost"].filesystem_no_deprecated = True
        self.options["boost"].filesystem_version = 4
        self.options["boost"].zlib = False
        self.options["boost"].bzip2 = False
        self.options["boost"].lzma = False
        self.options["boost"].zstd = False
        self.options["boost"].without_atomic = False
        self.options["boost"].without_charconv = True
        self.options["boost"].without_chrono = True
        self.options["boost"].without_cobalt = True
        self.options["boost"].without_container = True
        self.options["boost"].without_context = False
        self.options["boost"].without_contract = True
        self.options["boost"].without_coroutine = True
        # without_date_time is set to False to workaround https://github.com/conan-io/conan-center-index/issues/26890
        self.options["boost"].without_date_time = False
        self.options["boost"].without_exception = True
        self.options["boost"].without_fiber = True
        self.options["boost"].without_filesystem = False
        self.options["boost"].without_graph = True
        self.options["boost"].without_graph_parallel = True
        self.options["boost"].without_iostreams = True
        self.options["boost"].without_json = True
        self.options["boost"].without_locale = True
        self.options["boost"].without_log = True
        self.options["boost"].without_math = True
        self.options["boost"].without_mpi = True
        self.options["boost"].without_nowide = True
        self.options["boost"].without_process = False
        self.options["boost"].without_program_options = True
        self.options["boost"].without_python = True
        self.options["boost"].without_random = True
        self.options["boost"].without_regex = True
        self.options["boost"].without_serialization = True
        self.options["boost"].without_stacktrace = True
        self.options["boost"].without_system = False
        self.options["boost"].without_test = True
        self.options["boost"].without_thread = True
        self.options["boost"].without_timer = True
        self.options["boost"].without_type_erasure = True
        self.options["boost"].without_url = True
        self.options["boost"].without_wave = True
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
        if self.settings.os != "Windows":
            # Tweaking libcurl on Windows seems to be very brittle
            self.options["libcurl"].with_dict = False
            self.options["libcurl"].with_file = False
            self.options["libcurl"].with_ftp = False
            self.options["libcurl"].with_gopher = False
            self.options["libcurl"].with_imap = False
            self.options["libcurl"].with_ldap = False
            self.options["libcurl"].with_mqtt = False
            self.options["libcurl"].with_pop3 = False
            self.options["libcurl"].with_rtsp = False
            self.options["libcurl"].with_smb = False
            self.options["libcurl"].with_smtp = False
            self.options["libcurl"].with_telnet = False
            self.options["libcurl"].with_tftp = False
            self.options["libcurl"].with_zlib = True
            self.options["libcurl"].with_zstd = True
            self.options["libcurl"].with_ntlm = False
            self.options["libcurl"].with_ntlm_wb = False
            self.options["libcurl"].with_cookies = False
            self.options["libcurl"].with_verbose_debug = False
            self.options["libcurl"].with_unix_sockets = False
            self.options["libcurl"].with_verbose_strings = False
            self.options["libcurl"].with_form_api = False
            self.options["libcurl"].with_websocket = False
        self.options["opentelemetry-cpp"].with_otlp_http_compression = True
        self.options["opentelemetry-cpp"].with_no_deprecated_code = True
        self.options["opentelemetry-cpp"].with_jaeger = False
        self.options["opentelemetry-cpp"].with_zipkin = False
        self.options["spdlog"].header_only = True
        self.options["zstd"].build_programs = False
