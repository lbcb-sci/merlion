// Copyright (c) 2021 Robert Vaser

#ifndef MERLION_STACK_HPP_
#define MERLION_STACK_HPP_

#include <cstdint>
#include <utility>
#include <vector>

#include "biosoup/nucleic_acid.hpp"
#include "biosoup/overlap.hpp"
#include "cereal/cereal.hpp"
#include "cereal/access.hpp"
#include "cereal/types/vector.hpp"
#include "cereal/types/utility.hpp"

namespace merlion {

class Stack {
 public:
  explicit Stack(const biosoup::NucleicAcid& na);

  Stack(const Stack&) = default;
  Stack& operator=(const Stack&) = default;

  Stack(Stack&&) = default;
  Stack& operator=(Stack&&) = default;

  ~Stack() = default;

  std::uint32_t id() const {
    return id_;
  }

  bool is_invalid() const {
    return is_invalid_;
  }

  void set_is_invalid() {
    is_invalid_ = true;
  }

  bool is_contained() const {
    return is_contained_;
  }

  void set_is_contained() {
    is_contained_ = true;
  }

  bool is_chimeric() const {
    return is_chimeric_;
  }

  void set_is_chimeric() {
    is_chimeric_ = true;
  }

  const std::vector<std::pair<std::uint32_t, std::uint32_t>>& layers() const {
    return layers_;
  }

  void AddLayers(
      std::vector<biosoup::Overlap>::const_iterator begin,
      std::vector<biosoup::Overlap>::const_iterator end);

  void AddLayer(const biosoup::Overlap& o);

 private:
  Stack() = default;

  template<class Archive>
  void serialize(Archive& archive) {  // NOLINT
    archive(
        CEREAL_NVP(id_),
        CEREAL_NVP(len_),
        CEREAL_NVP(is_invalid_),
        CEREAL_NVP(is_contained_),
        CEREAL_NVP(is_chimeric_),
        CEREAL_NVP(layers_));
  }

  friend cereal::access;

  std::uint32_t id_;
  std::uint32_t len_;
  bool is_invalid_;
  bool is_contained_;
  bool is_chimeric_;
  std::vector<std::pair<std::uint32_t, std::uint32_t>> layers_;
};

}  // namespace merlion

#endif  // MERLION_STACK_HPP_
