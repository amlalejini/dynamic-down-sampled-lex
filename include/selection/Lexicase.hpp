#pragma once

#include <functional>
#include <algorithm>

#include "emp/base/vector.hpp"
#include "emp/datastructs/vector_utils.hpp"
#include "emp/math/Random.hpp"
#include "emp/math/random_utils.hpp"

#include "BaseSelect.hpp"

namespace selection {

class LexicaseSelect : public BaseSelect {

public:
  using score_fun_t = std::function<double(void)>;

protected:
  emp::vector< emp::vector<score_fun_t> >& score_fun_sets;
  emp::Random& random;

  emp::vector< emp::vector<double> > score_table;  ///< Holds results of calling score functions

  emp::vector<size_t> fun_ordering;                 ///< Used internally to track function ordering. WARNING - Don't change values in this!
  emp::vector<size_t> all_candidates;               ///< Used internally to track complete set of candidate ids.

  emp::vector<size_t>& Select(
    size_t n
  ) {

    const size_t fun_cnt = score_table.size();
    emp_assert(fun_cnt > 0);
    const size_t num_candidates = score_table[0].size();
    emp_assert(num_candidates > 0);

    selected.resize(n, 0);

    // Update function ordering
    fun_ordering.resize(fun_cnt);
    std::iota(
      fun_ordering.begin(),
      fun_ordering.end(),
      0
    );

    all_candidates.resize(num_candidates);
    std::iota(
      all_candidates.begin(),
      all_candidates.end(),
      0
    );

    emp::vector<size_t> cur_pool, next_pool;
    for (size_t sel_i = 0; sel_i < n; ++sel_i) {
      // Randomize the function ordering
      emp::Shuffle(random, fun_ordering);
      // Step through each function
      cur_pool = all_candidates;
      int depth = 0;
      // For each function, filter the population down to only the best performers.
      for (size_t fun_id : fun_ordering) {
        ++depth;
        double max_score = score_table[fun_id][cur_pool[0]]; // Max score starts as the first candidate's score on this function.
        next_pool.push_back(cur_pool[0]);                    // Seed the keeper pool with the first candidate.

        for (size_t i = 1; i < cur_pool.size(); ++i) {
          const size_t cand_id = cur_pool[i];
          const double cur_score = score_table[fun_id][cand_id];
          if (cur_score > max_score) {
            max_score = cur_score;              // This is the new max score for this function
            next_pool.resize(1);                // Clear out candidates with the former max score for this function.
            next_pool[0] = cand_id;             // Add this candidate as the only one with the new max score.
          } else if (cur_score == max_score) {
            next_pool.emplace_back(cand_id);    // Same as current max score. Save this candidate too.
          }
        }
        // Make next_pool into new cur_pool; make cur_pool allocated space for next_pool
        std::swap(cur_pool, next_pool);
        next_pool.resize(0);
        if (cur_pool.size() == 1) break; // Step if we're down to just one candidate.
      }
      // Select a random survivor (all equal at this point)
      emp_assert(cur_pool.size() > 0);
      const size_t win_id = cur_pool[random.GetUInt(cur_pool.size())];
      selected[sel_i] = win_id;
    }
    return selected;
  }

public:
  // todo - max depth
  LexicaseSelect(
    emp::vector< emp::vector<score_fun_t> >& a_score_fun_sets,
    emp::Random& a_random
  ) :
    score_fun_sets(a_score_fun_sets),
    random(a_random)
  {
    name="LexicaseSelect";
  }

  /// Run lexicase selection, but use only specified function ids and select from only the specified set of candidates.
  /// n: number of candidates to select (size of the returned vector)
  /// fun_ids: functions that can be used for selection (must be a subset of available functions)
  /// cand_ids: candidate ids that can be selected (must be a subset of available candidates)
  emp::vector<size_t>& operator()(
    size_t n,
    const emp::vector<size_t>& fun_ids,
    const emp::vector<size_t>& cand_ids
  ) {
    emp_assert(fun_ids.size() > 0);
    emp_assert(cand_ids.size() > 0);

    const size_t num_candidates = cand_ids.size();
    const size_t fun_cnt = fun_ids.size();

    // Update the score table
    score_table.resize(fun_cnt);
    for (size_t fun_i = 0; fun_i < fun_cnt; ++fun_i) {
      score_table[fun_i].resize(num_candidates);
      const size_t fun_id = fun_ids[fun_i];
      for (size_t cand_i = 0; cand_i < num_candidates; ++cand_i) {
        const size_t cand_id = cand_ids[cand_i];
        emp_assert(cand_id < score_fun_sets.size());
        emp_assert(fun_id < score_fun_sets[cand_id].size());
        score_table[fun_i][cand_i] = score_fun_sets[cand_id][fun_id]();
      }
    }

    // Run lexicase selection using constructed score table.
    Select(n);

    // Populate selected with actual ids of selected candidates
    for (size_t sel_i = 0; sel_i < selected.size(); ++sel_i) {
      selected[sel_i] = cand_ids[selected[sel_i]];
    }

    return selected;
  }

  // Run lexicase selection w/all candidates and all functions (as determined by the function set).
  emp::vector<size_t>& operator()(size_t n) override {
    const size_t num_candidates = score_fun_sets.size(); // How many candidates are there to select from?
    emp_assert(num_candidates > 0);
    const size_t fun_cnt = score_fun_sets[0].size();

    // Update the score table.
    score_table.resize(fun_cnt);

    for (size_t fun_i = 0; fun_i < fun_cnt; ++fun_i) {
      score_table[fun_i].resize(num_candidates);
      for (size_t cand_i = 0; cand_i < num_candidates; ++cand_i) {
        emp_assert(fun_cnt == score_fun_sets[cand_i].size());
        score_table[fun_i][cand_i] = score_fun_sets[cand_i][fun_i]();
      }
    }

    return Select(n);
  }
};

} // End selection namespace