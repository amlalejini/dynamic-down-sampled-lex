#pragma once

#include <functional>
#include <algorithm>

#include "emp/base/vector.hpp"
#include "emp/math/Random.hpp"

#include "BaseSelect.hpp"

namespace selection {


class TournamentSelect : public BaseSelect {
public:
  using score_fun_t = std::function<double(void)>;

protected:
  emp::vector< emp::vector<score_fun_t> >& score_fun_sets;
  emp::Random& random;
  size_t tournament_size;

  emp::vector<double> agg_scores;

  emp::vector<size_t>& Select(size_t n) {
    const size_t num_candidates = agg_scores.size();
    emp_assert(num_candidates > 0);

    selected.resize(n, 0);
    emp::vector<size_t> entries(tournament_size, 0); // Tracks entries into each tournament.
    for (size_t t=0; t < n; ++t) {
      // form a tournament
      std::generate(
        entries.begin(),
        entries.end(),
        [this, num_candidates](){ return random.GetUInt(num_candidates); }
      );
      // pick an arbitrary winner
      size_t winner_id = entries[0];
      double winner_fit = agg_scores[entries[0]];
      // update arbitrary winner to real winner
      for (size_t i = 1; i < entries.size(); ++i) {
        const size_t entry_id = entries[i];
        const double entry_fit = agg_scores[entry_id];
        if (entry_fit > winner_fit) {
          winner_id = entry_id;
        }
      }
      // save winner of tournament
      selected[t] = winner_id;
    }
    return selected;
  }

public:
  TournamentSelect(
    emp::vector< emp::vector<score_fun_t> >& a_score_fun_sets,
    emp::Random& a_random,
    size_t a_tournament_size=4
  ) :
    score_fun_sets(a_score_fun_sets),
    random(a_random),
    tournament_size(a_tournament_size)
  { }

  emp::vector<size_t>& operator()(
    size_t n,
    const emp::vector<size_t>& fun_ids,
    const emp::vector<size_t>& cand_ids
  ) {
    emp_assert(tournament_size > 0, "Tournament size must be greater than 0.", tournament_size);
    emp_assert(tournament_size <= cand_ids.size(), "Tournament size should not exceed number of candidates for selection.", tournament_size, cand_ids.size());
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

    // Run tournament selection over aggregate scores
    Select(n);

    // Transform selected ids into candidate ids
    for (size_t sel_i = 0; sel_i < selected.size(); ++sel_i) {
      selected[sel_i] = cand_ids[selected[sel_i]];
    }
    return selected;
  }

  emp::vector<size_t>& operator()(size_t n) override {
    emp_assert(tournament_size > 0, "Tournament size must be greater than 0.", tournament_size);

    const size_t num_candidates = score_fun_sets.size();
    emp_assert(tournament_size <= num_candidates, "Tournament size should not exceed number of individuals that we can select from.", tournament_size, num_candidates);

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

}