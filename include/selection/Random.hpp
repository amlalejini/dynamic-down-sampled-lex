#pragma once

#include <algorithm>

#include "emp/math/Random.hpp"

#include "BaseSelect.hpp"

namespace selection {

class RandomSelect : public BaseSelect {
protected:
  emp::Random& random;
  size_t num_candidates;

public:
  RandomSelect(
    emp::Random& a_random,
    size_t a_num_candidates
  ) :
    random(a_random),
    num_candidates(a_num_candidates)
  { }

  void SetNumCandidates(size_t n) {
    emp_assert(n > 0);
    num_candidates = n;
  }

  size_t GetNumCandidates() const {
    return num_candidates;
  }

  emp::vector<size_t>& operator()(size_t n) override {
    emp_assert(num_candidates > 0);
    selected.resize(n, 0);
    std::generate(
      selected.begin(),
      selected.end(),
      [this]() { return random.GetUInt(num_candidates); }
    );
    return selected;
  }

  emp::vector<size_t>& operator()(
    size_t n,
    const emp::vector<size_t>& cand_ids
  ) {
    emp_assert(cand_ids.size() > 0);
    selected.resize(n, 0);
    std::generate(
      selected.begin(),
      selected.end(),
      [this, &cand_ids]() { return cand_ids[random.GetUInt(cand_ids.size())]; }
    );
    return selected;
  }

};

}