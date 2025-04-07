// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <nlohmann/json.hpp>
#include <string>
#include <utility>
#include <vector>

#include "hictk/tools/toml.hpp"

namespace hictk::tools::io {

namespace json {

[[nodiscard]] nlohmann::json reformat_nulls(nlohmann::json attributes);

}  // namespace json

namespace toml {

template <typename T>
[[nodiscard]] inline ::toml::array to_array(const std::vector<T>& v, bool sort = false) {
  if (sort) {
    auto sv = v;
    std::sort(sv.begin(), sv.end());
    return to_array(sv, false);
  }

  ::toml::array buff;
  for (const auto& x : v) {
    buff.push_back(x);
  }
  return buff;
}

[[nodiscard]] nlohmann::json toml_to_json(const ::toml::table& t);

[[nodiscard]] std::string format_to_json(
    const ::toml::table& attributes,
    const std::vector<std::pair<std::string, ::toml::table>>& nested_attributes);

[[nodiscard]] std::string sanitize_toml_section_title(std::string s);

[[nodiscard]] std::string format_to_toml(
    const ::toml::table& attributes,
    const std::vector<std::pair<std::string, ::toml::table>>& nested_attributes);

[[nodiscard]] std::string format_to_yaml(
    const ::toml::table& attributes,
    const std::vector<std::pair<std::string, ::toml::table>>& nested_attributes);
}  // namespace toml
}  // namespace hictk::tools::io
