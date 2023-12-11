// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// IWYU pragma: private, include "hictk/cooler.hpp"

#include <cstdint>
#include <string>
#include <string_view>
#include <variant>

#include "group.hpp"
#include "hictk/common.hpp"
#include "hictk/type_traits.hpp"

namespace hictk::cooler {

struct Attribute {
  // Variants are listed in order from the most common to the least common for perf. reasons
  // clang-format off
  using AttributeVar = std::variant<
      std::monostate,
      std::string,
      std::uint64_t,
      std::int64_t,
      double,
      std::uint32_t, std::uint16_t, std::uint8_t,
      std::int32_t, std::int16_t, std::int8_t,
      float,
      bool>;
  // clang-format on
  Attribute() = delete;

  // ParentObj e.g. DataSet, Group
  template <typename ParentObj>
  [[nodiscard]] static bool exists(ParentObj& h5obj, std::string_view key);

  template <typename T, typename ParentObj>
  static void write(ParentObj& h5obj, std::string_view key, const T& value,
                    bool overwrite_if_exists = false);
  template <typename ParentObj>
  [[nodiscard]] static auto read(const ParentObj& h5obj, std::string_view key,
                                 bool missing_ok = false) -> AttributeVar;
  template <typename T, typename ParentObj>
  [[nodiscard]] static T read(const ParentObj& h5obj, std::string_view key);

  template <typename T, typename ParentObj>
  [[nodiscard]] static std::vector<T> read_vector(const ParentObj& h5obj, std::string_view key);
  template <typename T, typename ParentObj>
  static void read_vector(const ParentObj& h5obj, std::string_view key, std::vector<T>& buff);

 private:
  template <std::size_t i = 1>  // i = 1 skips T=monostate
  [[nodiscard]] static auto read_variant(const HighFive::Attribute& attr) -> AttributeVar;
  template <typename T1, typename Tout, typename Tin = remove_cvref_t<T1>>
  [[nodiscard]] static Tout numeric_converter(T1& buff);
};
}  // namespace hictk::cooler

#include "./impl/attribute_impl.hpp"
