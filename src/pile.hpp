// Copyright (c) 2021 Robert Vaser

#ifndef MERLION_PILE_HPP_
#define MERLION_PILE_HPP_

#include <cstdint>
#include <utility>
#include <vector>

#include "stack.hpp"

namespace merlion {

class Pile {
 public:
  explicit Pile(const Stack& s);

  Pile(const Pile&) = default;
  Pile& operator=(const Pile&) = default;

  Pile(Pile&&) = default;
  Pile& operator=(Pile&&) = default;

  ~Pile() = default;

  std::uint32_t id() const {
    return id_;
  }

  std::uint16_t median() const {
    return median_;
  }

  bool is_chimeric() const {
    return is_chimeric_;
  }

  void FindMedian();

  // store chimeric regions given median coverage
  void FindChimericRegions(std::uint16_t median);

 private:
  using Region = std::pair<std::uint32_t, std::uint32_t>;

  std::vector<Region> FindSlopes(double q);

  std::vector<Region> MergeRegions(const std::vector<Region>& r);

  std::uint32_t id_;
  std::vector<std::uint16_t> data_;
  std::uint16_t median_;
  bool is_chimeric_;
  std::vector<Region> chimeric_regions_;
};

}  // namespace merlion

#endif  // MERLION_PILE_HPP_
