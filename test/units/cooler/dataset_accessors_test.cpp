// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>
#include <random>
#include <set>

#include "hictk/cooler/dataset.hpp"
#include "hictk/cooler/group.hpp"
#include "tmpdir.hpp"

namespace hictk::cooler::test::dataset {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler: dataset accessors", "[dataset][short]") {
  const auto path = testdir() / "test_dataset_accessors.cool";
  std::filesystem::remove(path);
  std::filesystem::copy_file(datadir / "cooler_test_file.cool", path);
  SECTION("read-only") {
    const RootGroup grp{HighFive::File(path.string()).getGroup("/")};

    Dataset dset{grp, "chroms/name"};

    CHECK(dset.size() == 20);

    CHECK(dset.get().getFile().getName() == path.string());

    CHECK(dset.file_name() == path.string());
    CHECK(dset.uri() == fmt::format(FMT_STRING("{}::/chroms/name"), path.string()));
    CHECK(dset.name() == "name");
    CHECK(dset.hdf5_path() == "/chroms/name");
    CHECK(dset.get_parent().hdf5_path() == "/");
  }
  SECTION("read-write") {
    const RootGroup grp{HighFive::File(path.string(), HighFive::File::ReadWrite).getGroup("/")};

    Dataset dset{grp, "chroms/name"};

    CHECK_NOTHROW(dset().resize(std::vector<std::size_t>{20}));
  }
}

}  // namespace hictk::cooler::test::dataset
