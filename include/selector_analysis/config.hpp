#pragma once

#include "emp/config/config.hpp"

namespace selector_analysis {

  EMP_BUILD_CONFIG(SelectorAnalysisConfig,
    GROUP(GENERAL, "General settings"),
    VALUE(SEED, int, 0, "Random number seed."),
    VALUE(POPULATIONS_FILE, std::string, "pops.csv", "Path to the file with populations to load into the analyzer."),
    VALUE(OUTPUT_DIR, std::string, "./output/", "What directory are we dumping all this data"),
    VALUE(GENS, size_t, 0, "Should we run selection forward in an 'ecology' mode?"),

    GROUP(SELECTION, "Selection settings"),
    VALUE(SELECTION_METHOD, std::string, "random", "Which selection method should be used?"),
    VALUE(SELECTION_ROUNDS, size_t, 1, "How many times to run selection scheme on each loaded population?"),
    VALUE(TOURNAMENT_SIZE, size_t, 4, "Size of tournaments for tournament-based selection schemes"),
    VALUE(ELITE_COUNT, size_t, 1, "For elite selection, top N candidates to select as parents."),

    GROUP(TEST_SAMPLING, "Sampling settings"),
    VALUE(TEST_SAMPLING_METHOD, std::string, "none", "Method of sampling test cases to use."),
    VALUE(TEST_SAMPLING_PROP, double, 0.1, "Proportion of tests to sample for selection." ),

    GROUP(PARTITIONING, "Partitioning settings"),
    VALUE(PARTITIONING_METHOD, std::string, "none", "Method of partitioning population and/or test cases"),
    VALUE(COHORT_PARTITIONING_PROP, double, 0.1, "Proportion of test set and population sizes to use for test and population partitions, respectively."),

  )

} // namespace diag
