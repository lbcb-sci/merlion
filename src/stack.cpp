// Copyright (c) 2021 Robert Vaser

#include "stack.hpp"

namespace merlion {

Stack::Stack(const biosoup::NucleicAcid& na)
    : id_(na.id),
      len_(na.inflated_len),
      is_invalid_(false),
      is_contained_(false),
      is_chimeric_(false),
      layers_() {
}

void Stack::AddLayers(
    std::vector<biosoup::Overlap>::const_iterator begin,
    std::vector<biosoup::Overlap>::const_iterator end) {
  while (begin != end) {
    AddLayer(*begin);
    ++begin;
  }
}

void Stack::AddLayer(const biosoup::Overlap& o) {
  layers_.emplace_back(
      o.lhs_id == id_ ? o.lhs_begin : o.rhs_begin,
      o.lhs_id == id_ ? o.lhs_end   : o.rhs_end);
}

}  // namespace merlion
