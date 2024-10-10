// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <filesystem>
#include <highfive/H5File.hpp>
#include <vector>

#include "hictk/cooler/dataset.hpp"
#include "hictk/cooler/group.hpp"
#include "tmpdir.hpp"

namespace hictk::cooler::test::dataset {

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)
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

// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace hictk::cooler::test::dataset
