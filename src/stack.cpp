// Copyright (c) 2021 Robert Vaser

#include <algorithm>

#include "stack.hpp"

namespace merlion {

Stack::Stack(const biosoup::NucleicAcid& na)
    : id_(na.id),
      len_(na.inflated_len),
      layers_(),
      is_chimeric_(false) {
}

void Stack::AddLayer(const biosoup::Overlap& o) {
  if (id_ == o.lhs_id) {
    layers_.emplace_back(o.lhs_begin, o.lhs_end);
  } else if (id_ == o.rhs_id) {
    layers_.emplace_back(o.rhs_begin, o.rhs_end);
  }
}

void Stack::AddLayers(
    std::vector<biosoup::Overlap>::const_iterator begin,
    std::vector<biosoup::Overlap>::const_iterator end) {
  while (begin != end) {
    AddLayer(*begin);
    ++begin;
  }
}

void Stack::SortLayers() {
  std::sort(layers_.begin(), layers_.end());
}

}  // namespace merlion
