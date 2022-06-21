#pragma once

#include <map>
#include <functional>
#include <algorithm>

#include "emp/base/vector.hpp"

#include "BaseSelect.hpp"

namespace selection {

class EliteSelect : public BaseSelect {
public:
  using score_fun_t = std::function<double(void)>;

protected:
  emp::vector< emp::vector<score_fun_t> >& score_fun_sets;
  emp::vector<double> agg_scores;
  size_t elite_count;                  ///< How many distinct candidates should be chosen (in rank order by score)

  emp::vector<size_t>& Select(size_t n) {
    const size_t num_candidates = agg_scores.size();
    emp_assert(elite_count <= num_candidates, elite_count, num_candidates);

    selected.resize(n, 0);
    std::multimap<double, size_t> fit_map;
    for (size_t id = 0; id < num_candidates; ++id) {
      fit_map.insert(
        std::make_pair(agg_scores[id], id)
      );
    }

    // Identify the elites
    emp::vector<size_t> elites(elite_count, 0);
    auto m = fit_map.rbegin();
    for (size_t i = 0; i < elite_count; ++i) {
      elites[i] = m->second;
      ++m;
    }
    // Fill selected with elites
    for (size_t i = 0; i < n; ++i) {
      selected[i] = elites[i % elite_count];
    }

    return selected;
  }

public:
  EliteSelect(
    emp::vector< emp::vector<score_fun_t> >& a_score_fun_sets,
    size_t a_elite_count=1
  ) :
    score_fun_sets(a_score_fun_sets),
    elite_count(a_elite_count)
  { }

  emp::vector<size_t>& operator()(
    size_t n,
    const emp::vector<size_t>& fun_ids,
    const emp::vector<size_t>& cand_ids
  ) {
    const size_t num_candidates = cand_ids.size();
    const size_t fun_cnt = fun_ids.size();

    // Update agg scores to hold sum of relevant functions for each candidate.
    agg_scores.resize(num_candidates, 0.0);
    for (size_t cand_i=0; cand_i < num_candidates; ++cand_i) {
      const size_t cand_id = cand_ids[cand_i];
      agg_scores[cand_i] = 0.0;
      for (size_t fun_i=0; fun_i < fun_cnt; ++fun_i) {
        const size_t fun_id = fun_ids[fun_i];
        emp_assert(cand_id < score_fun_sets.size());
        emp_assert(fun_id < score_fun_sets[cand_id].size());
        agg_scores[cand_i] += score_fun_sets[cand_id][fun_id]();
      }
    }

    // Run tournament over aggregate scores
    Select(n);

    // Transform selected ids into candidate ids
    for (size_t sel_i = 0; sel_i < selected.size(); ++sel_i) {
      selected[sel_i] = cand_ids[selected[sel_i]];
    }
    return selected;
  }

  emp::vector<size_t>& operator()(size_t n) override {
    const size_t num_candidates = score_fun_sets.size();

    // Update agg_scores table
    agg_scores.resize(num_candidates, 0);
    for (size_t cand_id=0; cand_id < num_candidates; ++cand_id) {
      agg_scores[cand_id] = 0.0;
      for (size_t fun_id=0; fun_id < score_fun_sets[cand_id].size(); ++fun_id) {
        agg_scores[cand_id] += score_fun_sets[cand_id][fun_id]();
      }
    }

    // Run tournament selection on updated aggregate scores
    return Select(n);
  }

};

} // namespace selection