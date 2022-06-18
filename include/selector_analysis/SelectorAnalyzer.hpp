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
  // "lexicase",
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
  // "none",
  // "down-sample",
  // "cohorts",
  // "maxmin-sample",
  // "maxmin-fullinfo-sample"
  // ....
};

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

  PopulationSet pop_set;
  size_t cur_selection_round=0;
  size_t cur_pop_idx=0;

  void Setup();

  void SetupPopulations();
  void SetupDataTracking();
  void SetupSelection();


public:

  SelectorAnalyzer(
    const config_t& in_config
  ) :
    config(in_config),
    random(config.SEED())
  {
    Setup();
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

void SelectorAnalyzer::SetupSelection() {
  // TODO -- setup actual selection scheme + sampling + partitioning?
}

void SelectorAnalyzer::SetupDataTracking() {
  // TODO
}

void SelectorAnalyzer::Run() {
  // For each population, run selection scheme for configured number of replicates.
  for (cur_pop_idx=0; cur_pop_idx < pop_set.GetSize(); ++cur_pop_idx) {

    for (cur_selection_round=0; cur_selection_round < config.SELECTION_ROUNDS(); ++cur_selection_round) {
      // todo - run Selection on current population
      // todo - output data
    }

  }


}

} // namespace selector_analysis