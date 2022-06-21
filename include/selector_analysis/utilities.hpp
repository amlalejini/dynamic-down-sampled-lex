#pragma once

#include "emp/base/vector.hpp"
#include "emp/math/math.hpp"

namespace selector_analysis {

double VecDistDiffLengths(
  const emp::vector<double>& a,
  const emp::vector<double>& b
) {
  auto& larger = (a.size() >= b.size()) ? a : b;
  auto& smaller = (a.size() >= b.size()) ? b : a;
  double total = 0.0;
  for (size_t i = 0; i < smaller.size(); ++i) {
    total += emp::Abs(smaller[i] - larger[i]);
  }
  for (size_t i = smaller.size(); i < larger.size(); ++i) {
    total += emp::Abs(larger[i]);
  }
  return total;
}

double VecDist(
  const emp::vector<double>& a,
  const emp::vector<double>& b
) {
  if (a.size() != b.size()) {
    return VecDistDiffLengths(a, b);
  }
  emp_assert(a.size() == b.size());
  double total = 0.0;
  for (size_t i = 0; i < a.size(); ++i) {
    total += emp::Abs(a[i] - b[i]);
  }
  return total;
}

}