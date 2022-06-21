#pragma once

// Standard includes
#include <unordered_set>
#include <functional>

// Empirical includes
#include "emp/math/Random.hpp"
#include "emp/datastructs/set_utils.hpp"

// Local includes
#include "selection/SelectionSchemes.hpp"
#include "selector_analysis/config.hpp"
#include "selector_analysis/PopulationSet.hpp"

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
const std::unordered_set<std::string> valid_partitioning_methods = {
  "none",
  "random-cohort"
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

  PopulationSet pop_set;
  size_t cur_selection_round=0;
  size_t cur_pop_idx=0;
  size_t cur_pop_size=0;
  size_t cur_total_tests=0;
  emp::vector<size_t> cur_selected;

  void Setup();

  void SetupPopulations();
  void SetupDataTracking();

  // note - might be able to make sampling/partitioning more flexible by using structs/classes for each sampling/partitioning method (where base version does nothing)
  void SetupTestSampling();
  void SetupTestSamplingNone();
  void SetupTestSamplingRandom();

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
  SetupSelection();
  // Configure testcase sampling
  SetupTestSampling();
  // Configure partitoning method
  SetupPartitioning();

  // TODO - finish setup
  std::cout << "Done setting up SelectorAnalyzer" << std::endl;
}

void SelectorAnalyzer::SetupPopulations() {
  std::cout << "Loading populations from " << config.POPULATIONS_FILE() << std::endl;
  // Load populations from file
  pop_set.LoadFromCSV(config.POPULATIONS_FILE());
  std::cout << "  ...successfully loaded " << pop_set.GetSize() << " population(s)." << std::endl;
}

void SelectorAnalyzer::SetupDataTracking() {
  // TODO

  // Info to track:
  // - Evaluations? --> How many test scores would have need to be evaluated by selection scheme?
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
  } else {
    emp_assert(false, "Unimplemented test sampling method.", config.TEST_SAMPLING_METHOD());
  }
}

void SelectorAnalyzer::SetupTestSamplingNone() {
  // Configure do sample function
  do_sample_tests_fun = [this]() {
    // Use all test cases for selection
    sel_valid_tests.resize(pop_set.GetPop(cur_pop_idx).GetNumTestCases());
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
    emp_assert((config.TEST_SAMPLING_PROP() * (double)cur_total_tests) > 0);
    const size_t sample_size = (size_t)(config.TEST_SAMPLING_PROP() * (double)cur_total_tests);

    // Use all test cases for selection
    // emp::vector<size_t> test_ids(cur_total_tests, 0);
    sel_valid_tests.resize(cur_total_tests, 0);
    std::iota(
      sel_valid_tests.begin(),
      sel_valid_tests.end(),
      0
    );
    emp::Shuffle(random, sel_valid_tests);
    sel_valid_tests.resize(sample_size);
    // TODO - test sample!
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
    // std::cout << "  [do partitioning]" << std::endl;
    // Partition population and tests into an equal number of cohorts
    emp_assert(config.COHORT_PARTITIONING_PROP() > 0 and config.COHORT_PARTITIONING_PROP() <= 1.0);
    const size_t num_valid_tests = sel_valid_tests.size();
    const size_t num_valid_candidates = sel_valid_candidates.size();
    emp_assert(num_valid_tests > 0);
    // std::cout << "    num valid tests: " << num_valid_tests << std::endl;
    // Calculate population cohort partition sizes
    emp_assert((size_t)(config.COHORT_PARTITIONING_PROP() * (double)num_valid_candidates) > 0, "Too small of a population to ensure one individual per cohort.");
    const size_t base_pop_cohort_size = (size_t)(config.COHORT_PARTITIONING_PROP() * (double)num_valid_candidates);
    const size_t num_pop_cohorts = num_valid_candidates / base_pop_cohort_size;

    // Calculate test case cohort partition sizes
    emp_assert((size_t)(config.COHORT_PARTITIONING_PROP() * (double)num_valid_tests) > 0, "Too few tests (post-sampling) to ensure at least one test per cohort.");
    const size_t base_test_cohort_size = (size_t)(config.COHORT_PARTITIONING_PROP() * (double)num_valid_tests);

    // std::cout << "    base pop cohort size: " << base_pop_cohort_size << std::endl;
    // std::cout << "    base test cohort size: " << base_test_cohort_size << std::endl;
    // std::cout << "    num cohorts: " << num_pop_cohorts << std::endl;

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

    // std::cout << "    leftover pop: " << leftover_pop << std::endl;
    // std::cout << "    leftover tests: " << leftover_tests << std::endl;

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
      // std::cout << "    Pop cohort " << cohort_i << " ("<<pop_cohort.size()<<"): " << pop_cohort << std::endl;
    }
    // Should not be any leftover candidates because we based cohort sizing off of population size.
    emp_assert(leftover_pop == 0);
    // But, there might be leftover tests, so we need to distribute them evenly across cohorts.
    size_t cohort_i = 0;
    for (; cur_test_i < tests.size(); ++cur_test_i) {
      sel_test_partitions[cohort_i%num_pop_cohorts].emplace_back(tests[cur_test_i]);
      ++cohort_i;
    }

    // for (size_t cohort_i = 0; cohort_i < num_pop_cohorts; ++cohort_i) {
    //   std::cout << "    Test cohort " << cohort_i << " (" << sel_test_partitions[cohort_i].size() << "): " << sel_test_partitions[cohort_i] << std::endl;
    // }
    emp_assert(cur_cand_i == candidates.size(), cur_cand_i, candidates.size()); // Ensure that we used all of the candidates.
    emp_assert(cur_test_i == num_valid_tests, cur_test_i, num_valid_tests, tests.size()); // Ensure that we used all of the tests.
  };

  // Setup run selection routine
  run_selection_routine = [this]() {
    do_sample_tests_fun(); // Sample tests (if necessary)
    do_partition_fun(); // Partition (if necessary)
    std::cout << "  - used tests: " << sel_valid_tests << std::endl;
    std::cout << "  - selectable candidates: " << sel_valid_candidates << std::endl;
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
    std::cout << "  - selected candidates: " << cur_selected << std::endl;
    emp_assert(cur_selected.size() == cur_pop_size, cur_selected.size(), cur_pop_size);
  };
}

void SelectorAnalyzer::SetupPartitioningNone() {
  // do_partition_fun does nothing.
  do_partition_fun = [this]() { };
  // Setup default run selection routine
  run_selection_routine = [this]() {
    do_sample_tests_fun(); // Sample tests (if necessary)
    std::cout << "  - used tests: " << sel_valid_tests << std::endl;
    std::cout << "  - selectable candidates: " << sel_valid_candidates << std::endl;
    const auto& selected = do_selection_fun(cur_pop_size, sel_valid_tests, sel_valid_candidates);
    cur_selected.resize(selected.size(), 0);
    std::copy(
      selected.begin(),
      selected.end(),
      cur_selected.begin()
    );
    emp_assert(cur_selected.size() == cur_pop_size);
    std::cout << "  - selected candidates: " << cur_selected << std::endl;
  };
}

void SelectorAnalyzer::SetupPop(size_t pop_id) {
  emp_assert(pop_id < pop_set.GetSize());
  cur_pop_idx = pop_id;
  auto& cur_pop = pop_set.GetPop(cur_pop_idx);
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
        return pop_set.GetPop(cur_pop_idx).GetOrg(org_id).agg_score;
      }
    );
    // Wire test score function set
    test_score_fun_sets.emplace_back();
    const size_t fun_set_size = cur_pop.GetOrg(org_id).test_case_scores.size();
    emp_assert(cur_total_tests == fun_set_size);
    for (size_t fun_i = 0; fun_i < fun_set_size; ++fun_i) {
      test_score_fun_sets[org_id].emplace_back(
        [this, org_id, fun_i]() -> double {
          return pop_set.GetPop(cur_pop_idx).GetOrg(org_id).test_case_scores[fun_i];
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

}

void SelectorAnalyzer::Run() {
  // For each population, run selection scheme for configured number of replicates.
  for (size_t pop_i=0; pop_i < pop_set.GetSize(); ++pop_i) {
    std::cout << "Analzying pop " << pop_i << std::endl;
    // Update pop info
    SetupPop(pop_i); // TODO - turn this into a signal?
    for (cur_selection_round=0; cur_selection_round < config.SELECTION_ROUNDS(); ++cur_selection_round) {
      // std::cout << "  selection round " << cur_selection_round << std::endl;
      run_selection_routine();
      // std::cout << "  - selected (" << cur_selected.size() << "): " << cur_selected << std::endl;
      // todo - output data

      // TODO - support running for N "generations?"
    }

  }

}

} // namespace selector_analysis