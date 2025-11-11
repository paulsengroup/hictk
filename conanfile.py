# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

from conan import ConanFile
from conan.tools.microsoft import is_msvc

required_conan_version = ">=2"


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
        "with_cli_tool_deps": [True, False],
        "with_benchmark_deps": [True, False],
        "with_arrow": [True, False],
        "with_eigen": [True, False],
        "with_telemetry_deps": [True, False],
        "with_unit_testing_deps": [True, False],
        "with_fuzzy_testing_deps": [True, False],
    }

    default_options = {
        "with_cli_tool_deps": True,
        "with_benchmark_deps": False,
        "with_arrow": False,
        "with_eigen": False,
        "with_telemetry_deps": True,
        "with_unit_testing_deps": True,
        "with_fuzzy_testing_deps": False,
    }

    generators = "CMakeDeps"

    @property
    def _with_abseil(self) -> bool:
        return self._with_opentelemetry

    @property
    def _with_arrow(self) -> bool:
        return self.options.with_arrow

    @property
    def _with_boost(self) -> bool:
        return (self._with_arrow and is_msvc(self)) or self.options.with_fuzzy_testing_deps

    @property
    def _with_boost_header_only(self) -> bool:
        return not self.options.with_fuzzy_testing_deps

    @property
    def _with_bzip2(self) -> bool:
        return self._with_libarchive

    @property
    def _with_catch2(self) -> bool:
        return self.options.with_unit_testing_deps

    @property
    def _with_cli11(self) -> bool:
        opts = ("with_cli_tool_deps", "with_fuzzy_testing_deps", "with_benchmark_deps")
        return any(getattr(self.options, opt) for opt in opts)

    @property
    def _with_eigen(self) -> bool:
        return self.options.with_eigen

    @property
    def _with_libarchive(self) -> bool:
        return self.options.with_cli_tool_deps

    @property
    def _with_libcurl(self) -> bool:
        return self._with_opentelemetry

    @property
    def _with_lz4(self) -> bool:
        return self._with_libarchive

    @property
    def _with_lzo(self) -> bool:
        return self._with_libarchive

    @property
    def _with_nlohmann_json(self) -> bool:
        return self._with_opentelemetry or self.options.with_cli_tool_deps

    @property
    def _with_opentelemetry(self) -> bool:
        return self.options.with_telemetry_deps

    @property
    def _with_protobuf(self) -> bool:
        return self._with_opentelemetry

    @property
    def _with_pybind11(self) -> bool:
        return self.options.with_fuzzy_testing_deps

    @property
    def _with_re2(self) -> bool:
        return self._with_arrow

    @property
    def _with_tomlplusplus(self) -> bool:
        return self.options.with_cli_tool_deps

    @property
    def _with_xz_utils(self) -> bool:
        return self._with_libarchive

    def _configure_arrow(self):
        if not self._with_arrow:
            return
        self.options["arrow"].compute = True
        self.options["arrow"].filesystem_layer = False
        self.options["arrow"].parquet = False
        self.options["arrow"].with_boost = is_msvc(self)
        self.options["arrow"].with_re2 = True
        self.options["arrow"].with_thrift = False

    def _configure_boost(self):
        if not self._with_boost and not self._with_boost_header_only:
            return

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
        self.options["boost"].without_process = not self.options.with_fuzzy_testing_deps
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

        if self._with_boost_header_only:
            self.options["boost"].header_only = not self._with_boost

    def _configure_eigen(self):
        if not self._with_eigen:
            return
        self.options["eigen"].MPL2_only = True

    def _configure_libarchive(self):
        if not self._with_libarchive:
            return
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

    def _configure_libcurl(self):
        if not self._with_libcurl:
            return
        if self.settings.os == "Windows":
            # Tweaking libcurl on Windows seems to be very brittle
            return
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

    def _configure_opentelemetry(self):
        if not self._with_opentelemetry:
            return
        self.options["opentelemetry-cpp"].with_otlp_http_compression = True
        self.options["opentelemetry-cpp"].with_no_deprecated_code = True
        self.options["opentelemetry-cpp"].with_jaeger = False
        self.options["opentelemetry-cpp"].with_zipkin = False

    def build_requirements(self):
        # these are all transitive dependencies.
        # The only reason why they are defined here is to pin their recipe versions
        if self._with_abseil:
            # opentelemetry-cpp, re2
            self.requires("abseil/20250814.0#4e0fdd34a888b97aca482e648fc27a3b", force=True)

        if self._with_bzip2:
            # libarchive
            self.requires("bzip2/1.0.8#00b4a4658791c1f06914e087f0e792f5", force=True)

        if self._with_libcurl:
            # opentelemetry-cpp
            self.requires("libcurl/8.16.0#5b85991e709eafd00ab757b79814b3ff", force=True)

        if self._with_lz4:
            # libarchive
            self.requires("lz4/1.10.0#59fc63cac7f10fbe8e05c7e62c2f3504", force=True)

        if self._with_lzo:
            # libarchive
            self.requires("lzo/2.10#f00b10acc508cea70645727d970a23e1", force=True)

        if self._with_protobuf:
            # opentelemetry-cpp
            self.requires("protobuf/6.32.1#f481fd276fc23a33b85a3ed1e898b693", force=True)

        if self._with_re2:
            # arrow
            self.requires("re2/20251105#2eec620fc641f812a5837c642ecb7d66", force=True)

        if self._with_xz_utils:
            # libarchive
            self.requires("xz_utils/5.8.1#536e4c68ef30b7308b28a1c50e663c31", force=True)

        # hdf5, libarchive, and opentelemetry-cpp
        self.requires("zlib/1.3.1#b8bc2603263cf7eccbd6e17e66b0ed76", force=True)

    def requirements(self):
        self.requires("bshoshany-thread-pool/5.0.0#d94da300363f0c35b8f41b2c5490c94d")
        self.requires("concurrentqueue/1.0.4#1e48e1c712bcfd892087c9c622a51502")
        self.requires("fast_float/8.1.0#bbf67486bec084d167da0f3e13eee534")
        self.requires("fmt/12.1.0#50abab23274d56bb8f42c94b3b9a40c7", force=True)
        self.requires("hdf5/1.14.6#0b780319690d537e6cb0683244919955", force=True)
        self.requires("highfive/3.1.1#d0c724526ebc8ce396ffa1bf7f3c7b64")
        self.requires("libdeflate/1.25#49fcd3fe6c130c2ec5a01cabb0481ded")
        self.requires("parallel-hashmap/2.0.0#82acae64ffe2693fff5fb3f9df8e1746")
        self.requires("readerwriterqueue/1.0.6#aaa5ff6fac60c2aee591e9e51b063b83")
        self.requires("span-lite/0.11.0#519fd49fff711674cfed8cd17d4ed422")
        self.requires("spdlog/1.16.0#942c2c39562ae25ba575d9c8e2bdf3b6")
        self.requires("zstd/1.5.7#fde461c0d847a22f16d3066774f61b11", force=True)

        if self._with_arrow:
            self.requires("arrow/22.0.0#e46b173ba20adc478f7926495aeed142")

        if self._with_boost or self._with_boost_header_only:
            self.requires("boost/1.89.0#010f59feedfd171b15f467b42f723d13", force=True)

        if self._with_catch2:
            self.requires("catch2/3.11.0#a560b55882ff2deb1f4bafad832a1b98")

        if self._with_cli11:
            self.requires("cli11/2.6.0#41e625280f9c960e9970c1dc7d2c97b8")

        if self._with_eigen:
            self.requires("eigen/5.0.0#f7561f543f4aafd6d2dc1f6d677e3075", force=True)

        if self._with_libarchive:
            self.requires("libarchive/3.8.1#b42b1df243cee62014e0a884a9468b7a")

        if self._with_nlohmann_json:
            self.requires("nlohmann_json/3.12.0#2d634ab0ec8d9f56353e5ccef6d6612c", force=True)

        if self._with_opentelemetry:
            self.requires("opentelemetry-cpp/1.22.0#a7deeadbd209a75ccb42aeec309dfb51")

        if self._with_pybind11:
            self.requires("pybind11/3.0.1#81ed7e3cc5c945080013ca4c5c3abbbb")

        if self._with_tomlplusplus:
            self.requires("tomlplusplus/3.4.0#85dbfed71376fb8dc23cdcc0570e4727")

    def configure(self):
        self._configure_arrow()
        self._configure_boost()
        self.options["fmt"].header_only = True
        self.options["hdf5"].enable_cxx = False
        self.options["hdf5"].hl = False
        self.options["hdf5"].threadsafe = False
        self.options["hdf5"].parallel = False
        self.options["highfive"].with_boost = False
        self.options["highfive"].with_eigen = False
        self.options["highfive"].with_opencv = False
        self.options["highfive"].with_xtensor = False
        self._configure_eigen()
        self._configure_libarchive()
        self._configure_libcurl()
        self._configure_opentelemetry()
        self.options["spdlog"].header_only = True
        self.options["zstd"].build_programs = False
