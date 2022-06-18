#pragma once

#include "emp/config/config.hpp"

namespace selector_analysis {

  EMP_BUILD_CONFIG(SelectorAnalysisConfig,
    GROUP(GENERAL, "General settings"),
    VALUE(SEED, int, 0, "Random number seed."),
    VALUE(POPULATIONS_FILE, std::string, "pops.csv", "Path to the file with populations to load into the analyzer."),
    VALUE(SELECTION_ROUNDS, size_t, 10, "How many times to run selection scheme on each loaded population?"),
    VALUE(OUTPUT_DIR, std::string, "./output/", "What directory are we dumping all this data"),
  )

} // namespace diag
