// Copyright (c) 2021 Robert Vaser

#include <algorithm>
#include <deque>
#include <limits>

#include "pile.hpp"

namespace merlion {

constexpr std::uint32_t kPSS = 4;  // shrink 2 ^ kPSS times

constexpr double kCQ = 1.82;

template<typename T>
T Clamp(T v) {
  return v < std::numeric_limits<std::uint16_t>::max() ?
         v : std::numeric_limits<std::uint16_t>::max();
}

Pile::Pile(const Stack& s)
    : id_(s.id()),
      data_(s.len() >> kPSS),
      median_(0),
      is_chimeric_(false),
      chimeric_regions_() {
  std::vector<std::uint32_t> boundaries;
  for (const auto& it : s.layers()) {
    boundaries.emplace_back(((it.first  >> kPSS) + 1) << 1);
    boundaries.emplace_back(((it.second >> kPSS) - 1) << 1 | 1);
  }
  std::sort(boundaries.begin(), boundaries.end());

  std::uint32_t coverage = 0;
  std::uint32_t last_boundary = 0;
  for (const auto& it : boundaries) {
    if (coverage > 0) {
      for (std::uint32_t i = last_boundary; i < (it >> 1); ++i) {
        data_[i] = Clamp(data_[i] + coverage);
      }
    }
    last_boundary = it >> 1;
    coverage += it & 1 ? -1 : 1;
  }
}

void Pile::FindMedian() {
  decltype(data_) tmp(data_.begin(), data_.end());
  std::nth_element(tmp.begin(), tmp.begin() + tmp.size() / 2, tmp.end());
  median_ = tmp[tmp.size() / 2];
}

void Pile::FindChimericRegions(std::uint16_t median) {
  if (median_ < 4) {
    return;
  }

  auto slopes = FindSlopes(kCQ);
  if (slopes.empty()) {
    return;
  }

  for (std::uint32_t i = 0; i < slopes.size() - 1; ++i) {
    if (!(slopes[i].first & 1) && (slopes[i + 1].first & 1)) {
      chimeric_regions_.emplace_back(
          slopes[i].first >> 1,
          slopes[i + 1].second);
    }
  }
  chimeric_regions_ = MergeRegions(chimeric_regions_);

  auto is_chimeric_region = [&] (const Region& r) -> bool {
    for (std::uint32_t i = r.first; i <= r.second; ++i) {
      if (Clamp(data_[i] * kCQ) <= median) {
        return true;
      }
    }
    return false;
  };

  decltype(chimeric_regions_) dst;
  for (const auto& it : chimeric_regions_) {
    if (is_chimeric_region(it)) {
      dst.emplace_back(it);
    }
  }
  chimeric_regions_.swap(dst);

  if (!chimeric_regions_.empty()) {
    is_chimeric_ = true;
  }
}

std::vector<Pile::Region> Pile::FindSlopes(double q) {
  using Subpile = std::deque<std::pair<std::int32_t, std::uint16_t>>;
  auto subpile_add = [] (Subpile& s, std::uint16_t value, std::int32_t position) -> void {  // NOLINT
    while (!s.empty() && s.back().second <= value) {
      s.pop_back();
    }
    s.emplace_back(position, value);
  };
  auto subpile_update = [] (Subpile& s, std::int32_t position) {
    while (!s.empty() && s.front().first <= position) {
      s.pop_front();
    }
  };

  // find slopes
  std::vector<Region> dst;

  std::int32_t w = 847 >> kPSS;
  std::int32_t data_size = data_.size();

  Subpile left_subpile;
  std::uint32_t first_down = 0, last_down = 0;
  bool found_down = false;

  Subpile right_subpile;
  std::uint32_t first_up = 0, last_up = 0;
  bool found_up = false;

  // find slope regions
  for (std::int32_t i = 0; i < w; ++i) {
    subpile_add(right_subpile, data_[i], i);
  }
  for (std::int32_t i = 0; i < data_size; ++i) {
    if (i > 0) {
      subpile_add(left_subpile, data_[i - 1], i - 1);
    }
    subpile_update(left_subpile, i - 1 - w);

    if (i < data_size - w) {
      subpile_add(right_subpile, data_[i + w], i + w);
    }
    subpile_update(right_subpile, i);

    std::uint16_t d = Clamp(data_[i] * q);
    if (i != 0 && left_subpile.front().second > d) {
      if (found_down) {
        if (i - last_down > 1) {
          dst.emplace_back(first_down << 1 | 0, last_down);
          first_down = i;
        }
      } else {
        found_down = true;
        first_down = i;
      }
      last_down = i;
    }
    if (i != (data_size - 1) && right_subpile.front().second > d) {
      if (found_up) {
        if (i - last_up > 1) {
          dst.emplace_back(first_up << 1 | 1, last_up);
          first_up = i;
        }
      } else {
        found_up = true;
        first_up = i;
      }
      last_up = i;
    }
  }
  if (found_down) {
    dst.emplace_back(first_down << 1 | 0, last_down);
  }
  if (found_up) {
    dst.emplace_back(first_up << 1 | 1, last_up);
  }
  if (dst.empty()) {
    return dst;
  }

  // separate overlaping slopes
  while (true) {
    std::sort(dst.begin(), dst.end());

    bool is_changed = false;
    for (std::uint32_t i = 0; i < dst.size() - 1; ++i) {
      if (dst[i].second < (dst[i + 1].first >> 1)) {
        continue;
      }

      if (dst[i].first & 1) {
        right_subpile.clear();
        found_up = false;
        std::uint32_t subpile_begin = dst[i].first >> 1;
        std::uint32_t subpile_end = std::min(dst[i].second, dst[i + 1].second);

        for (std::uint32_t j = subpile_begin; j < subpile_end + 1; ++j) {
          subpile_add(right_subpile, data_[j], j);
        }
        for (std::uint32_t j = subpile_begin; j < subpile_end; ++j) {
          subpile_update(right_subpile, j);
          if (Clamp(data_[j] * q) < right_subpile.front().second) {
            if (found_up) {
              if (j - last_up > 1) {
                dst.emplace_back(first_up << 1 | 1, last_up);
                first_up = j;
              }
            } else {
              found_up = true;
              first_up = j;
            }
            last_up = j;
          }
        }
        if (found_up) {
          dst.emplace_back(first_up << 1 | 1, last_up);
        }
        dst[i].first = subpile_end << 1 | 1;

      } else {
        if (dst[i].second == (dst[i + 1].first >> 1)) {
          continue;
        }

        left_subpile.clear();
        found_down = false;

        std::uint32_t subpile_begin =
            std::max(dst[i].first >> 1, dst[i + 1].first >> 1);
        std::uint32_t subpile_end = dst[i].second;

        for (std::uint32_t j = subpile_begin; j < subpile_end + 1; ++j) {
          if (left_subpile.empty() == false &&
              Clamp(data_[j] * q) < left_subpile.front().second) {
            if (found_down) {
              if (j - last_down > 1) {
                dst.emplace_back(first_down << 1, last_down);
                first_down = j;
              }
            } else {
              found_down = true;
              first_down = j;
            }
            last_down = j;
          }
          subpile_add(left_subpile, data_[j], j);
        }
        if (found_down) {
          dst.emplace_back(first_down << 1, last_down);
        }
        dst[i].second = subpile_begin;
      }

      is_changed = true;
      break;
    }

    if (!is_changed) {
      break;
    }
  }

  // narrow slopes
  for (std::uint32_t i = 0; i < dst.size() - 1; ++i) {
    if ((dst[i].first & 1) && !(dst[i + 1].first & 1)) {
      std::uint32_t subpile_begin = dst[i].second;
      std::uint32_t subpile_end = dst[i + 1].first >> 1;

      if (subpile_end - subpile_begin > static_cast<std::uint32_t>(w)) {
        continue;
      }

      std::uint16_t max_coverage = 0;
      for (std::uint32_t j = subpile_begin + 1; j < subpile_end; ++j) {
        max_coverage = std::max(max_coverage, data_[j]);
      }

      std::uint32_t valid_point = dst[i].first >> 1;
      for (std::uint32_t j = dst[i].first >> 1; j <= subpile_begin; ++j) {
        if (max_coverage > Clamp(data_[j] * q)) {
          valid_point = j;
        }
      }
      dst[i].second = valid_point;

      valid_point = dst[i + 1].second;
      for (uint32_t j = subpile_end; j <= dst[i + 1].second; ++j) {
        if (max_coverage > Clamp(data_[j] * q)) {
          valid_point = j;
          break;
        }
      }
      dst[i + 1].first = valid_point << 1 | 0;
    }
  }

  return dst;
}

std::vector<Pile::Region> Pile::MergeRegions(const std::vector<Pile::Region>& src) {  // NOLINT
  std::vector<Region> dst;
  std::vector<bool> is_merged(src.size(), 0);
  for (std::uint32_t i = 0; i < src.size(); ++i) {
    if (is_merged[i]) {
      continue;
    }
    Region r = src[i];
    while (true) {
      is_merged[i] = false;
      for (std::uint32_t j = i + 1; j < src.size(); ++j) {
        if (is_merged[j]) {
          continue;
        }
        if (r.first < src[j].second && r.second > src[j].first) {
          is_merged[i] = true;
          is_merged[j] = true;
          r.first = std::min(r.first, src[j].first);
          r.second = std::max(r.second, src[j].second);
        }
      }
      if (!is_merged[i]) {
        break;
      }
    }
    dst.emplace_back(r);
  }
  return dst;
}

}  // namespace merlion
