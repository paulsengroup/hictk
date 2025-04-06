// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <iostream>

#include "hictk/balancing/ice.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/binary_buffer.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/common.hpp"
#include "hictk/cooler.hpp"
#include "hictk/default_delete_libdeflate.hpp"
#include "hictk/default_delete_zstd.hpp"
#include "hictk/expected_values_aggregator.hpp"
#include "hictk/file.hpp"
#include "hictk/filestream.hpp"
#include "hictk/fmt.hpp"
#include "hictk/genomic_interval.hpp"
#include "hictk/genomic_units.hpp"
#include "hictk/git.hpp"
#include "hictk/hash.hpp"
#include "hictk/hic.hpp"
#include "hictk/license.hpp"
#include "hictk/numeric_utils.hpp"
#include "hictk/numeric_variant.hpp"
#include "hictk/pixel.hpp"
#include "hictk/reference.hpp"
#include "hictk/static_binary_buffer.hpp"
#include "hictk/string.hpp"
#include "hictk/tmpdir.hpp"
#include "hictk/transformers.hpp"
#include "hictk/type_pretty_printer.hpp"
#include "hictk/type_traits.hpp"
#include "hictk/version.hpp"

int main() { std::cout << hictk::config::version::str_long() << "\n"; }
