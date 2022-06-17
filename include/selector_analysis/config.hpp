#ifndef SELECTOR_ANALYSIS_CONFIG_HPP
#define SELECTOR_ANALYSIS_CONFIG_HPP

#include "emp/config/config.hpp"

namespace selector_analysis {

  EMP_BUILD_CONFIG(SelectorAnalysisConfig,
    GROUP(GENERAL, "General settings"),
    VALUE(SEED, int, 0, "Random number seed."),
    VALUE(OUTPUT_DIR, std::string, "./output/", "What directory are we dumping all this data")
  )

} // namespace diag

#endif