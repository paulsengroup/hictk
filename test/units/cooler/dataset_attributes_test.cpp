// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <filesystem>
#include <highfive/H5File.hpp>
#include <string>
#include <variant>

#include "hictk/cooler/dataset.hpp"
#include "hictk/cooler/group.hpp"
#include "tmpdir.hpp"

namespace hictk::cooler::test::dataset {

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)
TEST_CASE("Cooler: dataset attributes", "[dataset][short]") {
  SECTION("read") {
    const auto path = datadir / "test_read_attrs.h5";

    const RootGroup grp{HighFive::File(path.string()).getGroup("/")};
    const Dataset dset{grp, "dst"};

    CHECK(dset.has_attribute("std::string"));
    CHECK(dset.read_attribute<std::string>("std::string") == "abc");
    auto attr = dset.read_attribute("std::string");
    CHECK(std::holds_alternative<std::string>(attr));
    CHECK_THROWS(dset.read_attribute("invalid"));
    attr = dset.read_attribute("invalid", true);
    CHECK(std::holds_alternative<std::monostate>(attr));
  }

  SECTION("write") {
    const auto path = testdir() / "test_dataset_write_attr.h5";

    const RootGroup grp{HighFive::File(path.string(), HighFive::File::Truncate).getGroup("/")};
    Dataset dset{grp, "int", std::uint8_t{}};

    dset.write_attribute("attr", 123);
    CHECK(dset.read_attribute<std::int32_t>("attr") == 123);

    CHECK_THROWS(dset.write_attribute("attr", -1, false));
    dset.write_attribute("attr", -1, true);
    CHECK(dset.read_attribute<std::int32_t>("attr") == -1);
  }
}  // NOLINT(clang-analyzer-cplusplus.NewDeleteLeaks)

// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace hictk::cooler::test::dataset
