#ifndef SELECTOR_ANALYSIS_ANALYZER_HPP
#define SELECTOR_ANALYSIS_ANALYZER_HPP

// Standard includes

// Empirical includes
#include "emp/math/Random.hpp"

// Local includes
#include "selection/SelectionSchemes.hpp"
#include "selector_analysis/config.hpp"
#include "selector_analysis/PopulationSet.hpp"

namespace selector_analysis {

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

  void Setup();

  void SetupPopulations();


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
  // Configure selection(?)

  // TODO - finish setup
  std::cout << "Done setting up SelectorAnalyzer" << std::endl;
}

void SelectorAnalyzer::SetupPopulations() {
  // Load populations from file
  pop_set.LoadFromCSV(config.POPULATIONS_FILE());
}

} // namespace selector_analysis


#endif