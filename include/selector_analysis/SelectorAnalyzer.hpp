#pragma once

// Standard includes
#include <unordered_set>
#include <functional>
#include <sys/stat.h>
#include <filesystem>

// Empirical includes
#include "emp/data/DataFile.hpp"
#include "emp/datastructs/set_utils.hpp"
#include "emp/math/Random.hpp"
#include "emp/bits/BitVector.hpp"
#include "emp/math/stats.hpp"

// Local includes
#include "selection/SelectionSchemes.hpp"
#include "selector_analysis/config.hpp"
#include "selector_analysis/PopulationSet.hpp"
#include "selector_analysis/utilities.hpp"

namespace selector_analysis {

const std::unordered_set<std::string> valid_selection_methods = {
  "lexicase",
  "tournament",
  "random",
  "elite",
  // "lexicase-even-leads",
  // ...
  // "non-dominated-elite",
  // "non-dominated-tournament",
  // "none"
};

const std::unordered_set<std::string> valid_test_case_sampling_methods = {
  "none",
  "random-sample",
  "maxmin-sample"
  // "maxmin-fullinfo-sample"
  // ....
};

const std::unordered_set<std::string> valid_partitioning_methods = {
  "none",
  "random-cohort"
};

template<typename CONFIG_T>
void SnapshotConfig(
  const CONFIG_T& cfg,
  const std::string& file_path
) {
  emp::DataFile snapshot_file(file_path);
  std::function<std::string(void)> get_param;
  std::function<std::string(void)> get_value;
  snapshot_file.AddFun<std::string>(
    [&get_param]() { return get_param(); },
    "parameter"
  );
  snapshot_file.AddFun<std::string>(
    [&get_value]() { return get_value(); },
    "value"
  );
  snapshot_file.PrintHeaderKeys();

  for (const auto& entry : cfg) {
    get_param = [&entry]() { return entry.first; };
    get_value = [&entry]() { return emp::to_string(entry.second->GetValue()); };
    snapshot_file.Update();
  }
}

// Run one selector over all loaded populations N number of times.
// Output statistics about which candidates were selected.
class SelectorAnalyzer {
public:
  using config_t = SelectorAnalysisConfig;

protected:

  const config_t& config;
  emp::Random random;

  emp::Ptr<selection::BaseSelect> selector=nullptr;
  std::function<emp::vector<size_t>&(size_t /* n */, emp::vector<size_t>& /* functions */, emp::vector<size_t>& /* candidates */)> do_selection_fun;
  std::function<void(void)> run_selection_routine;

  std::function<void(void)> do_sample_tests_fun;
  std::function<void(void)> do_partition_fun;

  emp::vector<size_t> sel_valid_candidates; ///< IDs of candidates that are valid to be selected (useful for partitioning)
  emp::vector<size_t> sel_valid_tests;      ///<

  emp::vector< emp::vector<size_t> > sel_candidate_partitions;
  emp::vector< emp::vector<size_t> > sel_test_partitions;

  emp::vector< std::function<double(void)> > aggregate_score_funs;                ///< One function for each individual in the population.
  emp::vector< emp::vector< std::function<double(void)> > > test_score_fun_sets;  ///< One set of functions for each individual in the population.

  std::string output_dir;
  emp::Ptr<emp::DataFile> summary_file=nullptr;
  // TODO - config snapshot output

  PopulationSet pop_set;
  Population cur_pop;
  Population next_pop;
  int cur_gen=0;
  size_t cur_selection_round=0;
  size_t cur_pop_idx=0;
  size_t cur_pop_size=0;
  size_t cur_total_tests=0;
  emp::vector<size_t> cur_selected;

  // Containers for organizing statistics -- note, may want to organize this a bit better/more generically later
  struct PopStatistics {
    std::unordered_set<size_t> tests_covered;
    size_t test_coverage;

    // set of phenotypes (as bitstrings?)
    double phenotype_entropy;
    double phenotype_richness;
    double max_aggregate_score;

    // TODO - single-test specialist count?
    // TODO - max test coverage?
    // TODO - calc pareto front size?

    void Reset() {
      tests_covered.clear();
      test_coverage=0;
      phenotype_richness=0.0;
      phenotype_richness=0;
      max_aggregate_score=0.0;
    }

    void Calculate(const Population& pop) {
      Reset();
      const size_t num_tests = pop.GetNumTestCases();
      const size_t pop_size = pop.GetSize();
      // Collect covered tests
      for (size_t test_i = 0; test_i < num_tests; ++test_i) {
        for (size_t cand_i = 0; cand_i < pop_size; ++cand_i) {
          const bool pass = pop.GetOrg(cand_i).test_case_scores[test_i] >= 1.0;
          if (pass) {
            tests_covered.emplace(test_i);
            break; // Don't need to look at anymore candidates, already covered.
          }
        }
      }
      test_coverage = tests_covered.size();

      // Collect binary pass/fail phenotypes
      emp::vector<emp::BitVector> phenotypes;
      for (size_t cand_i = 0; cand_i < pop_size; ++cand_i) {
        phenotypes.emplace_back(num_tests, false);
        auto& cur_phen = phenotypes.back();
        for (size_t test_i = 0; test_i < num_tests; ++test_i) {
          cur_phen[test_i] = pop.GetOrg(cand_i).test_case_scores[test_i] >= 1.0;
        }
      }

      phenotype_richness = emp::UniqueCount(phenotypes);
      phenotype_entropy = emp::ShannonEntropy(phenotypes);

      // Collect aggregate scores
      emp::vector<double> agg_scores(pop_size, 0.0);
      for (size_t cand_i = 0; cand_i < pop_size; ++cand_i) {
        agg_scores[cand_i] = pop.GetOrg(cand_i).agg_score;
      }
      max_aggregate_score = emp::FindMax(agg_scores);

    }
  };

  struct SelectedStatistics {
    emp::vector<size_t> selected_uids;
    size_t num_unique_cand_selected;
    double entropy_cand_selected;

    size_t num_unique_uids_selected;
    double entropy_uids_selected;

    void Reset() {
      selected_uids.clear();
      num_unique_cand_selected=0;
      entropy_cand_selected=0;
      num_unique_uids_selected=0;
      entropy_uids_selected=0;
    }

    void Calculate(const emp::vector<size_t>& selected, const Population& pop) {
      Reset();
      num_unique_cand_selected = emp::UniqueCount(selected);
      entropy_cand_selected = emp::ShannonEntropy(selected);
      selected_uids.resize(selected.size(), 0);
      for (size_t i = 0; i < selected_uids.size(); ++i) {
        selected_uids[i] = pop.GetOrg(selected[i]).uid;
      }
      num_unique_uids_selected = emp::UniqueCount(selected_uids);
      entropy_uids_selected = emp::ShannonEntropy(selected_uids);
    }

  };

  PopStatistics init_pop_stats;
  PopStatistics cur_pop_stats;
  SelectedStatistics cur_selected_stats;

  // struct SelectionSchemeStatistics {
  //   size_t
  // };

  void Setup();

  void SetupPopulations();
  void SetupDataCollection();

  // note - might be able to make sampling/partitioning more flexible by using structs/classes for each sampling/partitioning method (where base version does nothing)
  void SetupTestSampling();
  void SetupTestSamplingNone();
  void SetupTestSamplingRandom();
  void SetupTestSamplingMaxMin();

  void SetupPartitioning();
  void SetupPartitioningNone();
  void SetupPartitioningRandomCohort();

  void SetupSelection();
  void SetupLexicaseSelection();
  void SetupRandomSelection();
  void SetupTournamentSelection();
  void SetupEliteSelection();

  // Called at beginning of analyzing population
  void SetupPop(size_t pop_id);

  void DoReproduction();

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
    if (summary_file!=nullptr) summary_file.Delete();
  }

  void Run();

};

void SelectorAnalyzer::Setup() {
  std::cout << "Setting up SelectorAnalyzer" << std::endl;
  // Configure populations
  SetupPopulations();
  // Configure selection
  SetupSelection();
  // Configure testcase sampling
  SetupTestSampling();
  // Configure partitoning method
  SetupPartitioning();
  // Configure data tracking
  SetupDataCollection();

  // TODO - finish setup
  std::cout << "Done setting up SelectorAnalyzer" << std::endl;
}

void SelectorAnalyzer::SetupPopulations() {
  std::cout << "Loading populations from " << config.POPULATIONS_FILE() << std::endl;
  // Load populations from file
  pop_set.LoadFromCSV(config.POPULATIONS_FILE());
  std::cout << "  ...successfully loaded " << pop_set.GetSize() << " population(s)." << std::endl;
}

void SelectorAnalyzer::SetupDataCollection() {
  // TODO

  // Prepare output directory
  output_dir = config.OUTPUT_DIR();
  mkdir(output_dir.c_str(), ACCESSPERMS);
  if (output_dir.back() != '/') output_dir += '/';

  // Setup summary file
  summary_file = emp::NewPtr<emp::DataFile>(output_dir + "summary.csv");
  // population id
  summary_file->AddFun<size_t>(
    [this]() { return cur_pop_idx; },
    "pop_id"
  );
  // generation
  summary_file->AddFun<int>(
    [this]() { return cur_gen; },
    "generation"
  );
  // sel_replicate
  summary_file->AddFun<size_t>(
    [this]() { return cur_selection_round; },
    "sel_replicate"
  );
  // cur_pop_size
  summary_file->AddFun<size_t>(
    [this]() { return cur_pop_size; },
    "pop_size"
  );
  // cur_total_tests
  summary_file->AddFun<size_t>(
    [this]() { return cur_total_tests; },
    "total_test_cases"
  );
  // population_info (as json? "{}")
  summary_file->AddFun<std::string>(
    [this]() {
      std::ostringstream stream;
      stream << "\"{";
      bool first=true;
      for (auto& entry : cur_pop.GetPopInfo()) {
        if (!first) {
          stream << ",";
        }
        first=false;
        stream << entry.first << ":" << entry.second;
      }
      stream << "}\"";
      return stream.str();
    },
    "pop_info"
  );
  // coverage profile
  summary_file->AddFun<std::string>(
    [this]() {
      std::ostringstream stream;
      stream << "\"[";
      // cur_pop_stats
      for (size_t i = 0; i < cur_total_tests; ++i) {
        if (i) stream << ",";
        const bool covered_i = emp::Has(cur_pop_stats.tests_covered, i);
        stream << (size_t)covered_i;
      }
      stream << "]\"";
      return stream.str();
    },
    "pop_test_coverage_profile"
  );
  // test_coverage
  summary_file->AddFun<double>(
    [this]() {
      return cur_pop_stats.test_coverage;
    },
    "pop_test_coverage"
  );
  summary_file->AddFun<double>(
    [this]() {
      return cur_pop_stats.max_aggregate_score;
    },
    "pop_max_agg_score"
  );
  // phenotype_entropy
  summary_file->AddFun<double>(
    [this]() {
      return cur_pop_stats.phenotype_entropy;
    },
    "pop_phenotype_entropy"
  );
  // phenotype_richness
  summary_file->AddFun<double>(
    [this]() {
      return cur_pop_stats.phenotype_richness;
    },
    "pop_phenotype_richness"
  );

  // test_coverage_loss
  summary_file->AddFun<double>(
    [this]() {
      return (init_pop_stats.test_coverage - cur_pop_stats.test_coverage);
    },
    "pop_test_coverage_loss"
  );
  // phenotype_entropy_delta
  summary_file->AddFun<double>(
    [this]() {
      return (init_pop_stats.phenotype_entropy - cur_pop_stats.phenotype_entropy);
    },
    "pop_phenotype_entropy_diff"
  );
  // phenotype_richness_loss
  summary_file->AddFun<double>(
    [this]() {
      return (init_pop_stats.phenotype_richness - cur_pop_stats.phenotype_richness);
    },
    "pop_phenotype_entropy_loss"
  );
  // pop_max_agg_score_diff
  summary_file->AddFun<double>(
    [this]() {
      return (init_pop_stats.max_aggregate_score - cur_pop_stats.max_aggregate_score);
    },
    "pop_max_agg_score_diff"
  );

  summary_file->AddFun<size_t>(
    [this]() { return cur_selected_stats.num_unique_cand_selected; },
    "num_unique_cand_selected"
  );

  summary_file->AddFun<double>(
    [this]() { return cur_selected_stats.entropy_cand_selected; },
    "entropy_cand_selected"
  );

  summary_file->AddFun<double>(
    [this]() { return cur_selected_stats.num_unique_uids_selected; },
    "num_unique_uids_selected"
  );

  summary_file->AddFun<double>(
    [this]() { return cur_selected_stats.entropy_uids_selected; },
    "entropy_uids_selected"
  );


  // TODO - add additional statistics output

  summary_file->PrintHeaderKeys();

  // Info to track:
  // - Evaluations? --> How many test scores would have need to be evaluated by selection scheme?
  SnapshotConfig(config, output_dir+"run_config.csv");
}

// Needs to be run for each population to be analyzed
void SelectorAnalyzer::SetupSelection() {
  emp_assert(emp::Has(valid_selection_methods, config.SELECTION_METHOD()));
  std::cout << "Setting up selection method: " << config.SELECTION_METHOD() << std::endl;
  if (config.SELECTION_METHOD() == "lexicase") {
    SetupLexicaseSelection();
  } else if (config.SELECTION_METHOD() == "tournament") {
    SetupTournamentSelection();
  } else if (config.SELECTION_METHOD() == "elite") {
    SetupEliteSelection();
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
  // TODO - make this work with partitioning! (cur_pop_size not what we want for # of selected things)
  do_selection_fun = [this](size_t n, emp::vector<size_t>& tests, emp::vector<size_t>& candidates ) -> emp::vector<size_t>& {
    emp::Ptr<selection::LexicaseSelect> sel = selector.Cast<selection::LexicaseSelect>();
    return (*sel)(
      n,
      tests,
      candidates
    );
  };
}

void SelectorAnalyzer::SetupTournamentSelection() {
  selector = emp::NewPtr<selection::TournamentSelect>(
    test_score_fun_sets,
    random,
    config.TOURNAMENT_SIZE()
  );
  do_selection_fun = [this](
    size_t n,
    emp::vector<size_t>& tests,
    emp::vector<size_t>& candidates
  ) -> emp::vector<size_t>& {
    emp::Ptr<selection::TournamentSelect> sel = selector.Cast<selection::TournamentSelect>();
    return (*sel)(
      n,
      tests,
      candidates
    );
  };
}

void SelectorAnalyzer::SetupEliteSelection() {
  selector = emp::NewPtr<selection::EliteSelect>(
    test_score_fun_sets,
    config.ELITE_COUNT()
  );
  do_selection_fun = [this](
    size_t n,
    emp::vector<size_t>& tests,
    emp::vector<size_t>& candidates
  ) -> emp::vector<size_t>& {
    emp::Ptr<selection::EliteSelect> sel = selector.Cast<selection::EliteSelect>();
    return (*sel)(
      n,
      tests,
      candidates
    );
  };
}

void SelectorAnalyzer::SetupRandomSelection() {
  selector = emp::NewPtr<selection::RandomSelect>(
    random,
    cur_pop_size
  );
  do_selection_fun = [this](size_t n, emp::vector<size_t>& tests, emp::vector<size_t>& candidates ) -> emp::vector<size_t>& {
    emp::Ptr<selection::RandomSelect> sel = selector.Cast<selection::RandomSelect>();
    return (*sel)(n, candidates);
  };
}

void SelectorAnalyzer::SetupTestSampling() {
  emp_assert(emp::Has(valid_test_case_sampling_methods, config.TEST_SAMPLING_METHOD()));
  std::cout << "Setting up test sampling method: " << config.TEST_SAMPLING_METHOD() << std::endl;
  if (config.TEST_SAMPLING_METHOD() == "none") {
    SetupTestSamplingNone();
  } else if (config.TEST_SAMPLING_METHOD() == "random-sample") {
    SetupTestSamplingRandom();
  } else if (config.TEST_SAMPLING_METHOD() == "maxmin-sample") {
    SetupTestSamplingMaxMin();
  } else {
    emp_assert(false, "Unimplemented test sampling method.", config.TEST_SAMPLING_METHOD());
  }
}

void SelectorAnalyzer::SetupTestSamplingNone() {
  // Configure do sample function
  do_sample_tests_fun = [this]() {
    // Use all test cases for selection
    sel_valid_tests.resize(cur_pop.GetNumTestCases());
    std::iota(
      sel_valid_tests.begin(),
      sel_valid_tests.end(),
      0
    );
  };
}

void SelectorAnalyzer::SetupTestSamplingRandom() {
  // Configure do sample function
  emp_assert(config.TEST_SAMPLING_PROP() > 0);
  do_sample_tests_fun = [this]() {
    // Use a random subset of tests for selection
    emp_assert((int)(config.TEST_SAMPLING_PROP() * (double)cur_total_tests) > 0);
    const size_t sample_size = (size_t)(config.TEST_SAMPLING_PROP() * (double)cur_total_tests);
    // Begin with all test cases
    sel_valid_tests.resize(cur_total_tests, 0);
    std::iota(
      sel_valid_tests.begin(),
      sel_valid_tests.end(),
      0
    );
    // Shuffle and chop down to sample size
    emp::Shuffle(random, sel_valid_tests);
    sel_valid_tests.resize(sample_size);
  };
}

void SelectorAnalyzer::SetupTestSamplingMaxMin() {
  emp_assert(config.TEST_SAMPLING_PROP() > 0);
  emp_assert(config.MAXMIN_POP_PROP() > 0);

  do_sample_tests_fun = [this]() {
    // Calculate the test sample size
    emp_assert((int)(config.TEST_SAMPLING_PROP() * (double)cur_total_tests) > 0);
    const size_t test_sample_size = (size_t)(config.TEST_SAMPLING_PROP() * (double)cur_total_tests);
    // Calculate number of candidates to use in order to perform max min sampling (need at least 2)
    emp_assert((int)(config.MAXMIN_POP_PROP() * (double)cur_pop_size) > 1);
    const size_t cand_sample_size = (config.MAXMIN_POP_PROP() * (double)cur_pop_size);
    emp_assert(cand_sample_size > 1 && cand_sample_size <= cur_pop_size);

    // Pick a subset of possible parents to base maxmin sampling on.
    emp::vector<size_t> sampled_cand_ids(cur_pop_size, 0);
    std::iota(
      sampled_cand_ids.begin(),
      sampled_cand_ids.end(),
      0
    );
    emp::Shuffle(random, sampled_cand_ids);
    sampled_cand_ids.resize(cand_sample_size);
    // Build a matrix of sampled test profiles
    emp::vector< emp::vector<double> > sampled_test_profiles(
      cur_total_tests,
      emp::vector<double>(cand_sample_size, -1)
    );
    for (size_t cand_i = 0; cand_i < sampled_cand_ids.size(); ++cand_i) {
      const size_t cand_id = sampled_cand_ids[cand_i];
      const auto& cand_scores = cur_pop.GetOrg(cand_id).test_case_scores;
      for (size_t test_id = 0; test_id < cand_scores.size(); ++test_id) {
        sampled_test_profiles[test_id][cand_i] = cand_scores[test_id];
      }
    }
    sel_valid_tests = MaxMinSample(random, test_sample_size, sampled_test_profiles);
  };

}

void SelectorAnalyzer::SetupPartitioning() {
  emp_assert(emp::Has(valid_partitioning_methods, config.PARTITIONING_METHOD()));
  std::cout << "Setting up partitioning method: " << config.PARTITIONING_METHOD() << std::endl;
  if (config.PARTITIONING_METHOD() == "none") {
    SetupPartitioningNone();
  } else if (config.PARTITIONING_METHOD() == "random-cohort") {
    SetupPartitioningRandomCohort();
  } else {
    emp_assert(false, "Unimplemented partitioning method.", config.PARTITIONING_METHOD());
  }
}

void SelectorAnalyzer::SetupPartitioningRandomCohort() {
  // Configure do partitioning function
  do_partition_fun = [this]() {
    // Partition population and tests into an equal number of cohorts
    emp_assert(config.COHORT_PARTITIONING_PROP() > 0 and config.COHORT_PARTITIONING_PROP() <= 1.0);
    const size_t num_valid_tests = sel_valid_tests.size();
    const size_t num_valid_candidates = sel_valid_candidates.size();
    emp_assert(num_valid_tests > 0);
    // Calculate population cohort partition sizes
    emp_assert((size_t)(config.COHORT_PARTITIONING_PROP() * (double)num_valid_candidates) > 0, "Too small of a population to ensure one individual per cohort.");
    const size_t base_pop_cohort_size = (size_t)(config.COHORT_PARTITIONING_PROP() * (double)num_valid_candidates);
    const size_t num_pop_cohorts = num_valid_candidates / base_pop_cohort_size;

    // Calculate test case cohort partition sizes
    emp_assert((size_t)(config.COHORT_PARTITIONING_PROP() * (double)num_valid_tests) > 0, "Too few tests (post-sampling) to ensure at least one test per cohort.");
    const size_t base_test_cohort_size = (size_t)(config.COHORT_PARTITIONING_PROP() * (double)num_valid_tests);

    // Shuffle candidates (to randomly assign into cohorts)
    emp::vector<size_t> candidates(num_valid_candidates, 0);
    std::copy(
      sel_valid_candidates.begin(),
      sel_valid_candidates.end(),
      candidates.begin()
    );
    emp::Shuffle(random, candidates);

    // Shuffle test ids (to randomly assign into cohorts)
    emp::vector<size_t> tests(num_valid_tests, 0);
    std::copy(
      sel_valid_tests.begin(),
      sel_valid_tests.end(),
      tests.begin()
    );
    emp::Shuffle(random, tests);

    // Size the partitions
    sel_candidate_partitions.resize(num_pop_cohorts);
    int leftover_pop = num_valid_candidates - (num_pop_cohorts*base_pop_cohort_size);
    size_t cur_cand_i=0;

    sel_test_partitions.resize(num_pop_cohorts); // same number of population and test cohorts
    int leftover_tests = num_valid_tests - (num_pop_cohorts*base_test_cohort_size);
    size_t cur_test_i = 0;

    for (size_t cohort_i = 0; cohort_i < sel_candidate_partitions.size(); ++cohort_i) {
      // -- Fill population cohort --
      auto& pop_cohort = sel_candidate_partitions[cohort_i];
      const size_t pop_cohort_size = base_pop_cohort_size + ((size_t)(leftover_pop>0));
      leftover_pop -= (int)(leftover_pop > 0);
      pop_cohort.resize(pop_cohort_size, 0);
      for (size_t i = 0; i < pop_cohort_size; ++i) {
        emp_assert(cur_cand_i < candidates.size());
        pop_cohort[i] = candidates[cur_cand_i];
        ++cur_cand_i;
      }
      // -- Fill test cohort --
      auto& test_cohort = sel_test_partitions[cohort_i];
      const size_t test_cohort_size = base_test_cohort_size + ((size_t)leftover_tests>0);
      leftover_tests -= (int)(leftover_tests > 0);
      test_cohort.resize(test_cohort_size, 0);
      for (size_t i = 0; i < test_cohort_size; ++i) {
        emp_assert(cur_test_i < tests.size());
        test_cohort[i] = tests[cur_test_i];
        ++cur_test_i;
      }
    }
    // Should not be any leftover candidates because we based cohort sizing off of population size.
    emp_assert(leftover_pop == 0);
    // But, there might be leftover tests, so we need to distribute them evenly across cohorts.
    size_t cohort_i = 0;
    for (; cur_test_i < tests.size(); ++cur_test_i) {
      sel_test_partitions[cohort_i%num_pop_cohorts].emplace_back(tests[cur_test_i]);
      ++cohort_i;
    }

    emp_assert(cur_cand_i == candidates.size(), cur_cand_i, candidates.size()); // Ensure that we used all of the candidates.
    emp_assert(cur_test_i == num_valid_tests, cur_test_i, num_valid_tests, tests.size()); // Ensure that we used all of the tests.
  };

  // Setup run selection routine
  run_selection_routine = [this]() {
    do_sample_tests_fun(); // Sample tests (if necessary)
    do_partition_fun(); // Partition (if necessary)
    // std::cout << "  - used tests: " << sel_valid_tests << std::endl;
    // std::cout << "  - selectable candidates: " << sel_valid_candidates << std::endl;
    // std::cout << "  - selectable uids:";
    // for (size_t i = 0; i < cur_pop.GetSize(); ++i) {
    //   std::cout << " " << cur_pop.GetOrg(i).uid;
    // }
    // std::cout << std::endl;

    // Run selection for each cohort
    cur_selected.clear();
    emp_assert(sel_candidate_partitions.size() == sel_test_partitions.size());
    for (size_t cohort_i = 0; cohort_i < sel_candidate_partitions.size(); ++cohort_i) {
      const size_t pop_cohort_size = sel_candidate_partitions[cohort_i].size();
      const auto& selected = do_selection_fun(pop_cohort_size, sel_test_partitions[cohort_i], sel_candidate_partitions[cohort_i]);
      std::copy(
        selected.begin(),
        selected.end(),
        std::back_inserter(cur_selected)
      );
    }
    // std::cout << "  - selected candidates: " << cur_selected << std::endl;
    // std::cout << "  - selected uids: ";
    // for (size_t i = 0; i < cur_selected.size(); ++i) {
    //   std::cout << " " << cur_pop.GetOrg(cur_selected[i]).uid;
    // }
    // std::cout << std::endl;

    emp_assert(cur_selected.size() == cur_pop_size, cur_selected.size(), cur_pop_size);
  };
}

void SelectorAnalyzer::SetupPartitioningNone() {
  // do_partition_fun does nothing.
  do_partition_fun = [this]() { };
  // Setup default run selection routine
  run_selection_routine = [this]() {
    do_sample_tests_fun(); // Sample tests (if necessary)
    // std::cout << "  - used tests: " << sel_valid_tests << std::endl;
    // std::cout << "  - selectable candidates: " << sel_valid_candidates << std::endl;
    // std::cout << "  - selectable uids:";
    // for (size_t i = 0; i < cur_pop.GetSize(); ++i) {
    //   std::cout << " " << cur_pop.GetOrg(i).uid;
    // }
    // std::cout << std::endl;

    const auto& selected = do_selection_fun(cur_pop_size, sel_valid_tests, sel_valid_candidates);
    cur_selected.resize(selected.size(), 0);
    std::copy(
      selected.begin(),
      selected.end(),
      cur_selected.begin()
    );
    emp_assert(cur_selected.size() == cur_pop_size);

    // std::cout << "  - selected candidates: " << cur_selected << std::endl;
    // std::cout << "  - selected uids: ";
    // for (size_t i = 0; i < cur_selected.size(); ++i) {
    //   std::cout << " " << cur_pop.GetOrg(cur_selected[i]).uid;
    // }
    // std::cout << std::endl;
  };
}

void SelectorAnalyzer::SetupPop(size_t pop_id) {
  emp_assert(pop_id < pop_set.GetSize());
  cur_pop_idx = pop_id;
  cur_gen = 0;

  // Initialize cur_pop with specified loaded population.
  auto& loaded_pop = pop_set.GetPop(cur_pop_idx);
  cur_pop.SetNumTestCases(loaded_pop.GetNumTestCases());  // Update number of test cases
  cur_pop.UpdatePopSize(loaded_pop.GetSize());            // Update population size
  cur_pop.SetPopInfo(loaded_pop.GetPopInfo());            // Update population information
  next_pop.SetNumTestCases(loaded_pop.GetNumTestCases());
  next_pop.UpdatePopSize(loaded_pop.GetSize());
  next_pop.SetPopInfo(loaded_pop.GetPopInfo());
  // Update each organism
  for (size_t org_id = 0; org_id < cur_pop.GetSize(); ++org_id) {
    cur_pop.SetOrg(org_id, loaded_pop.GetOrg(org_id));
  }

  // auto& cur_pop = pop_set.GetPop(cur_pop_idx);
  cur_pop_size = cur_pop.GetSize();
  cur_total_tests = cur_pop.GetNumTestCases();
  // Wire up aggregate score and test score functions
  // TODO - test this! ==> just print out the results of running
  aggregate_score_funs.clear();
  test_score_fun_sets.clear();
  for (size_t org_id = 0; org_id < cur_pop_size; ++org_id) {
    // Wire aggregate score
    aggregate_score_funs.emplace_back(
      [this, org_id]() {
        return cur_pop.GetOrg(org_id).agg_score;
      }
    );
    // Wire test score function set
    test_score_fun_sets.emplace_back();
    const size_t fun_set_size = cur_pop.GetOrg(org_id).test_case_scores.size();
    emp_assert(cur_total_tests == fun_set_size);
    for (size_t fun_i = 0; fun_i < fun_set_size; ++fun_i) {
      test_score_fun_sets[org_id].emplace_back(
        [this, org_id, fun_i]() -> double {
          return cur_pop.GetOrg(org_id).test_case_scores[fun_i];
        }
      );
    }
  }

  sel_valid_candidates.resize(cur_pop_size, 0);
  std::iota(
    sel_valid_candidates.begin(),
    sel_valid_candidates.end(),
    0
  );
  sel_valid_tests.resize(cur_total_tests, 0);
  std::iota(
    sel_valid_tests.begin(),
    sel_valid_tests.end(),
    0
  );

  // Calculate statistics for initial loaded population.
  init_pop_stats.Calculate(cur_pop);
  cur_pop_stats.Calculate(cur_pop);


  // for (size_t org_i = 0; org_i < cur_pop.GetSize(); ++org_i) {
  //   auto& org = cur_pop.GetOrg(org_i);
  //   std::cout << org.test_case_scores << std::endl;
  // }
  // std::cout << "\"[";
  // // cur_pop_stats
  // for (size_t i = 0; i < cur_total_tests; ++i) {
  //   if (i) std::cout << ",";
  //   const bool covered_i = emp::Has(cur_pop_stats.tests_covered, i);
  //   std::cout << (size_t)covered_i;
  // }
  // std::cout << "]\"";

}

void SelectorAnalyzer::DoReproduction() {
  emp_assert(cur_selected.size() == next_pop.GetSize());
  emp_assert(cur_selected.size() == cur_pop.GetSize());
  // Reproduce selected organisms into next_pop
  for (size_t sel_i = 0; sel_i < cur_selected.size(); ++sel_i) {
    const size_t sel_id = cur_selected[sel_i];
    next_pop.SetOrg(sel_i, cur_pop.GetOrg(sel_id));
  }
  // Move next pop into cur pop
  // std::swap(cur_pop, next_pop); // TODO - test!
  for (size_t org_id = 0; org_id < cur_pop_size; ++org_id) {
    cur_pop.SetOrg(org_id, next_pop.GetOrg(org_id));
  }
}

void SelectorAnalyzer::Run() {
  // For each population, run selection scheme for configured number of replicates.
  for (size_t pop_i=0; pop_i < pop_set.GetSize(); ++pop_i) {
    std::cout << "Analzying pop " << pop_i << std::endl;
    for (cur_selection_round=0; cur_selection_round < config.SELECTION_ROUNDS(); ++cur_selection_round) {
      // Re-init population.

      SetupPop(pop_i); // TODO - turn this into a signal?
      // std::cout << "\nPop size: " <<  cur_pop.GetSize() << std::endl;
      cur_gen = 0;
      cur_selected_stats.Reset();
      summary_file->Update(); // Collect info about original population, mark as generation 0.
      // Increment generations, run first parent selection.
      ++cur_gen;
      run_selection_routine();
      cur_selected_stats.Calculate(cur_selected, cur_pop); // Run selected stats *before* doing reproduction
      // Update current population with selected individuals
      DoReproduction();
      // Calculate statistics on population post-selection
      cur_pop_stats.Calculate(cur_pop);
      // Output population statistics
      summary_file->Update();
      // todo - output data
      if (config.GENS() > 0) {
        ++cur_gen;
        for (; cur_gen <= (int)config.GENS(); ++cur_gen) {
          std::cout << " Gen=" << cur_gen << std::endl;
          run_selection_routine();
          cur_selected_stats.Calculate(cur_selected, cur_pop);
          DoReproduction();
          cur_pop_stats.Calculate(cur_pop);
          summary_file->Update();
        }
      }
    }

  }

}

} // namespace selector_analysis