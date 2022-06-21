#pragma once

#include "emp/base/vector.hpp"
#include "emp/math/math.hpp"
#include "emp/math/Random.hpp"
#include "emp/math/random_utils.hpp"

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

emp::vector<size_t> MaxMinSample(
  emp::Random& random,
  size_t sample_size,
  const emp::vector< emp::vector<double> >& profiles,
  int start_profile=-1
) {
  const size_t total_tests = profiles.size();
  // Initialize a cache for distances between test profiles
  emp::vector< emp::vector<double> > dist_cache(
    total_tests,
    emp::vector<double>(total_tests, -1)
  );
  emp::vector<size_t> available_ids(total_tests, 0);
  std::iota(
    available_ids.begin(),
    available_ids.end(),
    0
  );

  emp::vector<size_t> sampled_test_ids;
  if (start_profile >= 0) {
    emp_assert(start_profile < profiles.size());
    // Seed sample with specified test case.
    sampled_test_ids.emplace_back(start_profile);
    emp_assert((size_t)start_profile == available_ids[start_profile]);
    std::swap(available_ids[start_profile], available_ids[available_ids.size()-1]);
    available_ids.pop_back();
    emp::Shuffle(random, available_ids);
  } else {
    // Seed sample with random test case.
    emp::Shuffle(random, available_ids);
    sampled_test_ids.emplace_back(available_ids.back());
    available_ids.pop_back();
  }
  // emp::vector<size_t> sampled_test_ids({available_ids.back()});

  while (sampled_test_ids.size() < sample_size) {
    double maxmin_dist = -1;
    size_t maxmin_id = 0;
    size_t maxmin_avail_idx = 0;

    // For each possible test to be sampled, find its minimum distance to an already sampled test.
    for (size_t avail_idx = 0; avail_idx < available_ids.size(); ++avail_idx) {
      size_t cur_test_id = available_ids[avail_idx];
      double cur_min_dist = 0;
      const auto& cur_test_profile = profiles[cur_test_id];

      for (size_t sampled_idx = 0; sampled_idx < sampled_test_ids.size(); ++sampled_idx) {
        const size_t sampled_id = sampled_test_ids[sampled_idx];
        emp_assert(sampled_id != cur_test_id);
        double dist = 0.0;
        if (dist_cache[cur_test_id][sampled_id] != -1) {
          dist = dist_cache[cur_test_id][sampled_id];
        } else {
          const auto& sampled_test_profile = profiles[sampled_id];
          dist = VecDist(sampled_test_profile, cur_test_profile);
          dist_cache[cur_test_id][sampled_id] = dist;
          dist_cache[sampled_id][cur_test_id] = dist;
        }
        if (sampled_idx == 0 || dist < cur_min_dist) {
          cur_min_dist = dist;
        }
        if (dist == 0) break; // Not going to get bigger than zero.
      }

      if ( avail_idx==0 || cur_min_dist > maxmin_dist ) {
        maxmin_dist = cur_min_dist;
        maxmin_id = cur_test_id;
        maxmin_avail_idx = avail_idx;
      }
    }
    // move sampled id to back of available
    std::swap(available_ids[maxmin_avail_idx], available_ids[available_ids.size()-1]);
    emp_assert(available_ids.back() == maxmin_id);
    available_ids.pop_back();
    sampled_test_ids.emplace_back(maxmin_id);
  }

  return sampled_test_ids;
}

}