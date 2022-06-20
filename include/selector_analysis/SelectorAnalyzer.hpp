#pragma once

// Standard includes
#include <unordered_set>
#include <functional>

// Empirical includes
#include "emp/math/Random.hpp"

// Local includes
#include "selection/SelectionSchemes.hpp"
#include "selector_analysis/config.hpp"
#include "selector_analysis/PopulationSet.hpp"

namespace selector_analysis {

const std::unordered_set<std::string> valid_selection_methods = {
  "lexicase",
  "random"
  // "lexicase-even-leads",
  // ...
  // "elite",
  // "tournament",
  // "non-dominated-elite",
  // "non-dominated-tournament",
  // "random",
  // "none"
};

const std::unordered_set<std::string> valid_test_case_sampling_methods = {
  "none",
  "random-sample"
  // "cohorts",
  // "maxmin-sample",
  // "maxmin-fullinfo-sample"
  // ....
};

// TODO - Partitioning methods?

// Run one selector over all loaded populations N number of times.
// Output statistics about which candidates were selected.
class SelectorAnalyzer {
public:
  using config_t = SelectorAnalysisConfig;


protected:
  const config_t& config;
  emp::Random random;

  emp::Ptr<selection::BaseSelect> selector=nullptr;
  std::function<emp::vector<size_t>&(void)> do_selection_fun;

  emp::vector<size_t> sel_valid_candidates; ///< IDs of candidates that are valid to be selected (useful for partitioning)
  emp::vector<size_t> sel_valid_funs;       ///<

  emp::vector< std::function<double(void)> > aggregate_score_funs;                ///< One function for each individual in the population.
  emp::vector< emp::vector< std::function<double(void)> > > test_score_fun_sets;  ///< One set of functions for each individual in the population.

  PopulationSet pop_set;
  size_t cur_selection_round=0;
  size_t cur_pop_idx=0;
  size_t cur_pop_size=0;

  void Setup();

  void SetupPopulations();
  void SetupDataTracking();

  void SetupSelection();
  void SetupLexicaseSelection();
  void SetupRandomSelection();

  // Called at beginning of analyzing population
  void SetupPop(size_t pop_id);

public:

  SelectorAnalyzer(
    const config_t& in_config
  ) :
    config(in_config),
    random(config.SEED())
  {
    Setup();
  }

  ~SelectorAnalyzer() {
    if (selector!=nullptr) selector.Delete();
  }

  void Run();

};

void SelectorAnalyzer::Setup() {
  std::cout << "Setting up SelectorAnalyzer" << std::endl;
  // Configure populations
  SetupPopulations();
  // Configure selection
  // TODO

  // TODO - finish setup
  std::cout << "Done setting up SelectorAnalyzer" << std::endl;
}

void SelectorAnalyzer::SetupPopulations() {
  // Load populations from file
  pop_set.LoadFromCSV(config.POPULATIONS_FILE());
}

void SelectorAnalyzer::SetupDataTracking() {
  // TODO

  // Info to track:
  // - Evaluations? --> How many test scores would have need to be evaluated by selection scheme?
}

// Needs to be run for each population to be analyzed
void SelectorAnalyzer::SetupSelection() {
  // TODO -- setup actual selection scheme + sampling + partitioning?
  // TODO - any general-purpose wiring

  if (config.SELECTION_METHOD() == "lexicase") {
    SetupLexicaseSelection();
  } else if (config.SELECTION_METHOD() == "random") {
    SetupRandomSelection();
  } else {
    // code should never reach this else (unless I forget to add a selection scheme here that is in the valid selection method set)
    emp_assert(false, "Unimplemented selection scheme.", config.SELECTION_METHOD());
  }
}

void SelectorAnalyzer::SetupLexicaseSelection() {
  selector = emp::NewPtr<selection::LexicaseSelect>(
    test_score_fun_sets,
    random
  );
  do_selection_fun = [this]() -> emp::vector<size_t>& {
    emp::Ptr<selection::LexicaseSelect> sel = selector.Cast<selection::LexicaseSelect>();
    return (*sel)(
      cur_pop_size,
      sel_valid_funs,
      sel_valid_candidates
    );
  };
}

void SelectorAnalyzer::SetupRandomSelection() {
  // TODO
}

void SelectorAnalyzer::SetupPop(size_t pop_id) {
  emp_assert(pop_id < pop_set.GetSize());
  cur_pop_idx = pop_id;
  auto& cur_pop = pop_set.GetPop(cur_pop_idx);
  cur_pop_size = cur_pop.GetSize();
  // Wire up aggregate funs
  // TODO - test this!
  aggregate_score_funs.clear();
  for (size_t org_id = 0; org_id < cur_pop_size; ++org_id) {
    aggregate_score_funs.emplace_back(
      [this, org_id]() {
        return pop_set.GetPop(cur_pop_idx).GetOrg(org_id).agg_score;
      }
    );
  }

  // Wire up test score function sets
  // TODO - test this!
  test_score_fun_sets.clear();
  for (size_t org_id = 0; org_id < cur_pop_size; ++org_id) {
    test_score_fun_sets.emplace_back();
    const size_t fun_set_size = cur_pop.GetOrg(org_id).test_case_scores.size();
    emp_assert(cur_pop.GetNumTestCases() == fun_set_size);
    for (size_t fun_i = 0; fun_i < fun_set_size; ++fun_i) {
      test_score_fun_sets[org_id].emplace_back(
        [this, org_id, fun_i]() -> double {
          return pop_set.GetPop(cur_pop_idx).GetOrg(org_id).test_case_scores[fun_i];
        }
      );
    }
  }
}

void SelectorAnalyzer::Run() {
  // For each population, run selection scheme for configured number of replicates.
  for (size_t pop_i=0; pop_i < pop_set.GetSize(); ++pop_i) {
    // (1) Update pop info
    SetupPop(pop_i);
    // TODO - Setup sampling/partitioning?
    for (cur_selection_round=0; cur_selection_round < config.SELECTION_ROUNDS(); ++cur_selection_round) {
      // todo - run Selection on current population
      // todo - output data
    }

  }


}

} // namespace selector_analysis