// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include <cassert>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace hictk {
inline CoolerURI::CoolerURI(std::string_view p1, std::string_view p2)
    : CoolerURI(std::string{p1}, std::string{p2}) {}

inline CoolerURI::CoolerURI(std::string p1, std::string p2)
    : file_path(std::move(p1)), group_path(std::move(p2)) {}

inline CoolerURI::CoolerURI(std::pair<std::string_view, std::string_view> paths)
    : CoolerURI(paths.first, paths.second) {}

inline CoolerURI::CoolerURI(std::pair<std::string, std::string> paths)
    : CoolerURI(std::move(paths.first), std::move(paths.second)) {}

namespace internal {
constexpr std::pair<std::string_view, std::string_view> str_partition(std::string_view str,
                                                                      std::string_view delim) {
  assert(!delim.empty());

  if (str.empty()) {
    return {};
  }

  const auto p1 = str.find(delim);
  if (p1 == std::string_view::npos) {
    return {str, ""};
  }

  const auto p2 = p1 + delim.size();
  return {str.substr(0, p1), str.substr(p2)};
}
}  // namespace internal

inline CoolerURI parse_cooler_uri(std::string_view uri) {
  constexpr std::string_view separator = "::";
  auto [tok1, tok2] = internal::str_partition(uri, separator);

  if (tok1.empty()) {
    throw std::runtime_error(fmt::format(FMT_STRING("Invalid Cooler URI: \"{}\""), uri));
  }

  if (tok2.empty()) {
    return CoolerURI{std::string{uri}, "/"};
  }

  if (tok2.front() == '/') {
    return CoolerURI{std::string{tok1}, std::string{tok2}};
  }

  return CoolerURI{std::string{tok1}, "/" + std::string{tok2}};
}

}  // namespace hictk
