/// World that will manage solutions during the evolutionary run

#ifndef DIA_WORLD_HPP
#define DIA_WORLD_HPP

///< standard headers
#include <functional>
#include <map>
#include <set>
#include <fstream>
#include <string.h>
#include <set>
#include <string>
#include <algorithm>
#include <cmath>
#include <deque>
#include <filesystem>
#include <sys/stat.h>

///< empirical headers
#include "emp/Evolve/World.hpp"

///< experiment headers
#include "config.hpp"
#include "org.hpp"
#include "problem.hpp"
#include "selection.hpp"

namespace diag {

class DiagWorld : public emp::World<Org>
{
  // object types for consistency between working class
  public:
    ///< Org related

    // solution genome + diagnotic problem types
    using genome_t = emp::vector<double>;
    // score vector for a solution
    using score_t = emp::vector<double>;
    // boolean optimal vector per objective
    using optimal_t = emp::vector<bool>;
    // target vector type
    using target_t = emp::vector<double>;

    ///< selection related types

    // vector of position ids
    using ids_t = emp::vector<size_t>;
    // matrix of population score vectors
    using fmatrix_t = emp::vector<score_t>;
    // matrix of population genomes
    using gmatrix_t = emp::vector<genome_t>;
    // map holding population id groupings by fitness (keys in decending order)
    using fitgp_t = std::map<double, ids_t, std::greater<double>>;
    // vector of double vectors for K neighborhoods
    using neigh_t = emp::vector<score_t>;
    // vector of vector size_t for Pareto grouping
    using pareto_t = emp::vector<emp::vector<size_t>>;

    ///< world related types

    // evaluation function type
    using eval_t = std::function<double(Org &)>;
    // selection function type
    using sele_t = std::function<ids_t()>;

    ///< data tracking stuff (ask about)
    using nodef_t = emp::Ptr<emp::DataMonitor<double>>;
    using nodeo_t = emp::Ptr<emp::DataMonitor<size_t>>;
    using como_t = std::map<size_t, ids_t>;

  private:

    // experiment configurations
    DiaConfig & config;
    // enum class Scheme {
    //   TRUNCATION=0,
    //   TOURNAMENT=1,
    //   FITNESS_SHARING=2,
    //   LEXICASE=3,
    //   NONDOMINATED=4,
    //   NOVELTY=5
    // };

    // target vector
    target_t target;
    // vector holding population aggregate scores (by position id)
    score_t fit_vec;
    // vector holding parent solutions selected by selection scheme
    ids_t parent_vec;
    // novelty minimum
    double pmin = 0.0;
    // generations since solution added to archive
    size_t archive_gens = 0;


    // evaluation lambda we set
    eval_t evaluate;
    // selection lambda we set
    sele_t select;

    // select.h var
    emp::Ptr<Selection> selection;
    // problem.h var
    emp::Ptr<Diagnostic> diagnostic;

    ///< data file & node related variables
    // file we are working with
    emp::Ptr<emp::DataFile> data_file;       ///< Population-level data(?)
    emp::Ptr<emp::DataFile> elite_data_file; ///< Stores information about the elite organism
    std::string output_dir;
    // systematics tracking
    // emp::Ptr<systematics_t> sys_ptr;
    // node to track population fitnesses
    nodef_t pop_fit;
    // node to track population opitmized count
    nodeo_t pop_opti;
    // node to track parent fitnesses
    nodef_t pnt_fit;
    // node to track parent optimized count
    nodeo_t pnt_opti;
    // node to track streak counts
    nodeo_t pop_str;
    // // csv file to track best performing solutions
    // std::ofstream elite_csv;

    ///< data we are tracking during an evolutionary run
    // elite solution position
    size_t elite_pos;
    // common solution position
    size_t comm_pos;
    // optimal solution position
    size_t opti_pos;
    // streak solution position
    size_t strk_pos;
    // population activation gene vector
    optimal_t pop_acti_gene;

    // Pareto group count
    size_t pareto_cnt = 0;

    // novelty search archive
    std::deque<score_t> archive;
    // elite solution position
    double arc_elite = 0.0;
    // archive optimal trait vector
    optimal_t arc_opti_trt;
    // archive activation gene vector
    optimal_t arc_acti_gene;

    // lexicase - sampled test cases
    emp::vector<size_t> sampled_trait_ids;
    emp::vector<size_t> possible_trait_ids;
    std::function<void()> do_sample_traits;
    emp::vector< emp::vector<double> > test_score_matrix; // [test_id][org_id]

    void SetupTraitSampling();
    double PhenDist(
      const emp::vector<double>& a,
      const emp::vector<double>& b
    ) {
      if (a.size() != b.size()) {
        return PhenDistDiffLengths(a, b);
      }
      emp_assert(a.size() == b.size());
      double total = 0.0;
      for (size_t i = 0; i < a.size(); ++i) {
        total += emp::Abs(a[i] - b[i]);
      }
      return total;
    }
    double PhenDistDiffLengths(
      const emp::vector<double>& a,
      const emp::vector<double>& b
    ) {
      auto& larger = (a.size() >= b.size()) ? a : b;
      auto& smaller = (a.size() >= b.size()) ? b : a;
      double total = 0.0;
      for (size_t i = 0; i < smaller.size(); ++i) {
        total += emp::Abs(smaller[i] - larger[i]);
      }
      for (size_t i = smaller.size(); i < larger.size(); ++i) {
        total += emp::Abs(larger[i]);
      }
      return total;
    }

    // ---- call all functions to initiallize the world ----
    void Initialize();
    void SetOnUpdate();
    void SetMutation();
    void SetSelection();
    void SetOnOffspringReady();
    void SetEvaluation();
    void SetDataTracking();
    void PopulateWorld();

    // ---- principle steps during an evolutionary run ----
    void ResetData();
    void DoEvaluationStep();
    void DoSelectionStep();
    void DoReproductionStep();
    void RecordData();

    // ---- selection scheme implementations ----
    void SetupSelection_Truncation();
    void SetupSelection_Tournament();
    void SetupSelection_FitnessSharing();
    void SetupSelection_EpsilonLexicase();
    void SetupSelection_NonDominatedSorting();
    void SetupSelection_NoveltySearch();
    void SetupSelection_Lexicase();
    void SetupSelection_LexicaseEvenLead();
    // void DownSampledLexicase();

    // ---- evaluation function implementations ----
    void SetupDiagnostic_Exploitation();
    void SetupDiagnostic_StructuredExploitation();
    void SetupDiagnostic_StrongEcology();
    void SetupDiagnostic_Exploration();
    void SetupDiagnostic_WeakEcology();

    // ---- data tracking ----
    size_t UniqueObjective();
    size_t FindUniqueStart();
    void FindEverything();
    size_t ActivationGeneOverlap();

    // ---- helper functions ----
    // create a matrix of popultion score vectors
    fmatrix_t PopFitMat();

    // create matrix of population genomes
    gmatrix_t PopGenomes();

    // update archive
    bool ArchiveUpdate(const score_t & score, const fmatrix_t & dmat);
    void ArchiveDataUpdate(const size_t org_id);

  public:

    DiagWorld(
      DiaConfig & _config
    ) :
      config(_config)
    {
      // set random pointer seed
      this->NewRandom(config.SEED());
      // initialize the world
      Initialize();
    }

    ~DiagWorld() {
      selection.Delete();
      diagnostic.Delete();
      pop_fit.Delete();
      pop_opti.Delete();
      pnt_fit.Delete();
      pnt_opti.Delete();
      data_file.Delete();
      elite_data_file.Delete();
    }

    void RunStep();
    void Run();
};

void DiagWorld::RunStep() {
  this->Update();
}

void DiagWorld::Run() {
  for (size_t ud = 0; ud <= config.MAX_GENS(); ud++)
  {
    RunStep();
  }
}


///< functions called to setup the world

void DiagWorld::Initialize()
{
  std::cout << "==========================================" << std::endl;
  std::cout << "BEGINNING INITIAL SETUP" << std::endl;
  std::cout << "==========================================" << std::endl;

  // reset the world upon start
  Reset();
  // set world to well mixed so we don't over populate
  SetPopStruct_Mixed(true);

  // stuff we need to initialize for the experiment
  SetEvaluation();
  SetMutation();
  SetOnUpdate();
  SetDataTracking();
  SetSelection();
  SetOnOffspringReady();
  PopulateWorld();

  std::cout << "==========================================" << std::endl;
  std::cout << "FINISHED INITIAL SETUP" << std::endl;
  std::cout << "==========================================" << std::endl;
}

void DiagWorld::SetOnUpdate()
{
  std::cout << "------------------------------------------------" << std::endl;
  std::cout << "Setting OnUpdate function..." << std::endl;

  // set up the evolutionary algorithm
  OnUpdate([this](size_t gen)
  {
    // step 0: reset all data collection variables
    ResetData();

    // step 1: evaluate all solutions on diagnostic
    DoEvaluationStep();

    // step 2: select parent solutions for
    DoSelectionStep();

    // step 3: gather and record data
    RecordData();

    // step 4: reproduce and create new solutions
    DoReproductionStep();
  });

  std::cout << "Finished setting the OnUpdate function! \n" << std::endl;
}

void DiagWorld::SetMutation()
{
  std::cout << "------------------------------------------------" << std::endl;
  std::cout << "Setting mutation function..." << std::endl;

  // set the mutation function
  SetMutFun(
    [this](Org & org, emp::Random & random) {
      // number of mutations and solution genome
      size_t mcnt = 0;
      genome_t & genome = org.GetGenome();

      // quick checks
      emp_assert(genome.size() == config.OBJECTIVE_CNT());
      emp_assert(target.size() == config.OBJECTIVE_CNT());

      for (size_t i = 0; i < genome.size(); ++i) {
        // if we do a mutation at this objective
        if (random_ptr->P(config.MUTATE_PER())) {
          const double mut = random_ptr->GetRandNormal(config.MEAN(), config.STD());

          // mutation puts objective above target
          if (config.TARGET() < genome[i] + mut) {
            // we wrap it back around target value
            genome[i] = target[i] - (genome[i] + mut - target[i]);

          // mutation puts objective in the negatives
          } else if (genome[i] + mut < config.LOWER_BND()) {
            // genome[i] = std::abs(genome[i] + mut) + config.LOWER_BND();
            genome[i] = config.LOWER_BND();

          // else we can simply add mutation
          } else {
            genome[i] = genome[i] + mut;
          }

          ++mcnt;
        }
      }
      return mcnt;
    }
  );

  std::cout << "Mutation function set!\n" << std::endl;
}

void DiagWorld::SetSelection()
{
  std::cout << "------------------------------------------------" << std::endl;
  std::cout << "Setting Selection function..." << std::endl;

  selection = emp::NewPtr<Selection>(*random_ptr);
  std::cout << "Created selection" << std::endl;

  do_sample_traits = []() { /*do nothing by default*/; };

  if (config.SELECTION() == "truncation") {
    SetupSelection_Truncation();
  } else if (config.SELECTION() == "tournament") {
    SetupSelection_Tournament();
  } else if (config.SELECTION() == "fitness-sharing") {
    SetupSelection_FitnessSharing();
  } else if (config.SELECTION() == "lexicase") {
    SetupTraitSampling();
    SetupSelection_Lexicase();
  } else if (config.SELECTION() == "lexicase-eps") {
    SetupTraitSampling();
    SetupSelection_EpsilonLexicase();
  } else if (config.SELECTION() == "lexicase-even-lead") {
    SetupTraitSampling();
    SetupSelection_LexicaseEvenLead();
  } else if (config.SELECTION() == "nondominated-sorting") {
    SetupSelection_NonDominatedSorting();
  } else if (config.SELECTION() == "novelty") {
    SetupSelection_NoveltySearch();
  } else {
    std::cout << "ERROR UNKNOWN SELECTION CALL, " << config.SELECTION() << std::endl;
    emp_assert(false);
  }

  std::cout << "Finished setting the Selection function! \n" << std::endl;
}

void DiagWorld::SetupTraitSampling() {
  std::cout << "Setting up trait sampling." << std::endl;
  // Setup trait sampling
  sampled_trait_ids.resize(config.OBJECTIVE_CNT(), 0);
  std::iota(
    sampled_trait_ids.begin(),
    sampled_trait_ids.end(),
    0
  );
  possible_trait_ids.resize(config.OBJECTIVE_CNT(), 0);
  std::iota(
    possible_trait_ids.begin(),
    possible_trait_ids.end(),
    0
  );
  emp_assert(config.LEX_DS_RATE() <= 1.0 && config.LEX_DS_RATE() > 0.0);
  emp_assert((config.LEX_DS_RATE() * (double)config.OBJECTIVE_CNT()) > 0);

  if (config.LEX_DS_MODE() == "random") {
    std::cout << "  Sample type: random" << std::endl;
    do_sample_traits = [this]() {
      // calculate the sample size
      const size_t sample_size = (size_t)(config.LEX_DS_RATE() * (double)config.OBJECTIVE_CNT());
      sampled_trait_ids.resize(sample_size, 0);
      emp_assert(sampled_trait_ids.size() <= possible_trait_ids.size(), "The number of sampled traits should be <= than the total number of traits");
      emp::Shuffle(GetRandom(), possible_trait_ids);
      for (size_t i = 0; i < sampled_trait_ids.size(); ++i) {
        sampled_trait_ids[i] = possible_trait_ids[i];
      }
    };

  } else if (config.LEX_DS_MODE() == "maxmin-full") {
    // full info maxmin
    std::cout << "  Sample type: maxmin-full" << std::endl;
    // TODO - check that this works as expected
    do_sample_traits = [this]() {
      // std::cout << "--------SAMPLING---------" << std::endl;
      // std::cout << "Test profiles: " << std::endl;
      // for (size_t test_id = 0; test_id < config.OBJECTIVE_CNT(); ++test_id) {
      //   std::cout << "["<<test_id<<"]: " << test_score_matrix[test_id] << std::endl;
      // }

      // calculate the sample size
      const size_t sample_size = (size_t)(config.LEX_DS_RATE() * (double)config.OBJECTIVE_CNT());
      emp_assert(sampled_trait_ids.size() <= possible_trait_ids.size(), "The number of sampled traits should be <= than the total number of traits");

      emp::vector<size_t> available_ids(possible_trait_ids.size(), 0);
      std::copy(
        possible_trait_ids.begin(),
        possible_trait_ids.end(),
        available_ids.begin()
      );
      emp::Shuffle(GetRandom(), available_ids);
      std::unordered_set<size_t> included_ids;
      sampled_trait_ids.clear();
      sampled_trait_ids.emplace_back(available_ids.back());
      available_ids.pop_back();
      // std::cout << "First sampled trait: " << sampled_trait_ids[0] << std::endl;

      emp::vector< emp::vector<double> > dist_cache(possible_trait_ids.size(), emp::vector<double>(possible_trait_ids.size(), -1));

      while (sampled_trait_ids.size() < sample_size) {
        // std::cout << "Sampling test " << sampled_trait_ids.size() << std::endl;
        // Find available id with the maximum min distance to a chosen id
        double maxmin_dist = -1;
        size_t maxmin_id = 0;
        size_t maxmin_avail_idx = 0;

        // For each possible test to be sampled, find its minimum distance to an already sampled test.
        for (size_t avail_idx = 0; avail_idx < available_ids.size(); ++avail_idx) {
          size_t cur_id = available_ids[avail_idx];
          // std::cout << "  Analyzing test " << cur_id << std::endl;
          double cur_min_dist = 0; // find minimum distance between current availble id and all already sampled ids
          const auto& cur_test_profile = test_score_matrix[cur_id];
          // std::cout << "    Test profile" << test_score_matrix[cur_id] << std::endl;

          for (size_t sampled_idx = 0; sampled_idx < sampled_trait_ids.size(); ++sampled_idx) {
            const size_t sampled_id = sampled_trait_ids[sampled_idx];
            // std::cout << "    Comparing " << cur_id << " vs sampled " << sampled_id << std::endl;
            emp_assert(sampled_id != cur_id); // An available ID should never already be sampled.
            double dist=0.0;
            if (dist_cache[cur_id][sampled_id] != -1) {
              dist = dist_cache[cur_id][sampled_id];
            } else {
              const auto& sampled_test_profile = test_score_matrix[sampled_id];
              dist = PhenDist(sampled_test_profile, cur_test_profile);
              dist_cache[cur_id][sampled_id] = dist;
              dist_cache[sampled_id][cur_id] = dist;
            }
            // std::cout << "      dist = " << dist << std::endl;
            if (sampled_idx == 0 || dist < cur_min_dist) {
              cur_min_dist = dist;
            }
            if (dist == 0) break; // Not going to get bigger than zero.
          }

          // std::cout << "    min dist for this test " << cur_min_dist << std::endl;

          if ( avail_idx==0 || cur_min_dist > maxmin_dist ) {
            maxmin_dist = cur_min_dist;
            maxmin_id = cur_id;
            maxmin_avail_idx = avail_idx;
          }
        }

        // std::cout << "  found maxmin_dist ("<<maxmin_id<<") = " << maxmin_dist << std::endl;
        // move sampled id to back of available
        std::swap(available_ids[maxmin_avail_idx], available_ids[available_ids.size()-1]);
        emp_assert(available_ids.back() == maxmin_id);
        available_ids.pop_back();
        sampled_trait_ids.emplace_back(maxmin_id);
      }

      /////////


      // std::cout << "Cached distances: " << std::endl;
      // for (size_t test_id = 0; test_id < config.OBJECTIVE_CNT(); ++test_id) {
      //   std::cout << "["<<test_id<<"]: " << dist_cache[test_id] << std::endl;
      // }

      // std::cout << "Sampled test ids: " << sampled_trait_ids << std::endl;
      /////////

    };
  } else if (config.LEX_DS_MODE() == "maxmin-pop-sample") {
    // full info maxmin
    std::cout << "  Sample type: maxmin-pop-sample" << std::endl;
    // TODO - check that this works as expected
    do_sample_traits = [this]() {
      // std::cout << "--------SAMPLING---------" << std::endl;
      // std::cout << "Test profiles: " << std::endl;
      // for (size_t test_id = 0; test_id < config.OBJECTIVE_CNT(); ++test_id) {
      //   std::cout << "["<<test_id<<"]: " << test_score_matrix[test_id] << std::endl;
      // }

      // calculate the sample size
      const size_t sample_size = (size_t)(config.LEX_DS_RATE() * (double)config.OBJECTIVE_CNT());
      emp_assert(sampled_trait_ids.size() <= possible_trait_ids.size(), "The number of sampled traits should be <= than the total number of traits");

      emp::vector<size_t> available_ids(possible_trait_ids.size(), 0);
      std::copy(
        possible_trait_ids.begin(),
        possible_trait_ids.end(),
        available_ids.begin()
      );
      emp::Shuffle(GetRandom(), available_ids);
      std::unordered_set<size_t> included_ids;
      sampled_trait_ids.clear();
      sampled_trait_ids.emplace_back(available_ids.back());
      available_ids.pop_back();
      // std::cout << "First sampled trait: " << sampled_trait_ids[0] << std::endl;

      // Only base sample calculation on a sample of possible parents
      emp::vector<size_t> sampled_org_ids(GetSize(), 0);
      std::iota(
        sampled_org_ids.begin(),
        sampled_org_ids.end(),
        0
      );
      emp::Shuffle(GetRandom(), sampled_org_ids);
      const size_t org_sample_size = (size_t)(config.LEX_DS_POP_RATE() * (double)GetSize());
      sampled_org_ids.resize(org_sample_size);
      // Build matrix of sampled test profiles
      emp::vector< emp::vector<double> > sampled_test_profiles(possible_trait_ids.size(), emp::vector<double>(org_sample_size, -1));
      for (size_t org_idx = 0; org_idx < sampled_org_ids.size(); ++org_idx) {
        const size_t org_id = sampled_org_ids[org_idx];
        auto& org_scores = GetOrg(org_id).GetScore();
        for (size_t test_id = 0; test_id < org_scores.size(); ++test_id) {
          sampled_test_profiles[test_id][org_idx] = test_score_matrix[test_id][org_id];
        }
      }
      // Initialize a cache for distances between test profiles
      emp::vector< emp::vector<double> > dist_cache(possible_trait_ids.size(), emp::vector<double>(possible_trait_ids.size(), -1));

      while (sampled_trait_ids.size() < sample_size) {
        // std::cout << "Sampling test " << sampled_trait_ids.size() << std::endl;
        // Find available id with the maximum min distance to a chosen id
        double maxmin_dist = -1;
        size_t maxmin_id = 0;
        size_t maxmin_avail_idx = 0;

        // For each possible test to be sampled, find its minimum distance to an already sampled test.
        for (size_t avail_idx = 0; avail_idx < available_ids.size(); ++avail_idx) {
          size_t cur_id = available_ids[avail_idx];
          // std::cout << "  Analyzing test " << cur_id << std::endl;
          double cur_min_dist = 0; // find minimum distance between current availble id and all already sampled ids
          const auto& cur_test_profile = sampled_test_profiles[cur_id];
          // std::cout << "    Test profile" << test_score_matrix[cur_id] << std::endl;

          for (size_t sampled_idx = 0; sampled_idx < sampled_trait_ids.size(); ++sampled_idx) {
            const size_t sampled_id = sampled_trait_ids[sampled_idx];
            // std::cout << "    Comparing " << cur_id << " vs sampled " << sampled_id << std::endl;
            emp_assert(sampled_id != cur_id); // An available ID should never already be sampled.
            double dist=0.0;
            if (dist_cache[cur_id][sampled_id] != -1) {
              dist = dist_cache[cur_id][sampled_id];
            } else {
              const auto& sampled_test_profile = sampled_test_profiles[sampled_id];
              dist = PhenDist(sampled_test_profile, cur_test_profile);
              dist_cache[cur_id][sampled_id] = dist;
              dist_cache[sampled_id][cur_id] = dist;
            }
            // std::cout << "      dist = " << dist << std::endl;
            if (sampled_idx == 0 || dist < cur_min_dist) {
              cur_min_dist = dist;
            }
            if (dist == 0) break; // Not going to get bigger than zero.
          }

          // std::cout << "    min dist for this test " << cur_min_dist << std::endl;

          if ( avail_idx==0 || cur_min_dist > maxmin_dist ) {
            maxmin_dist = cur_min_dist;
            maxmin_id = cur_id;
            maxmin_avail_idx = avail_idx;
          }
        }

        // std::cout << "  found maxmin_dist ("<<maxmin_id<<") = " << maxmin_dist << std::endl;
        // move sampled id to back of available
        std::swap(available_ids[maxmin_avail_idx], available_ids[available_ids.size()-1]);
        emp_assert(available_ids.back() == maxmin_id);
        available_ids.pop_back();
        sampled_trait_ids.emplace_back(maxmin_id);
      }

      /////////


      // std::cout << "Cached distances: " << std::endl;
      // for (size_t test_id = 0; test_id < config.OBJECTIVE_CNT(); ++test_id) {
      //   std::cout << "["<<test_id<<"]: " << dist_cache[test_id] << std::endl;
      // }

      // std::cout << "Sampled test ids: " << sampled_trait_ids << std::endl;
      /////////

    };

  } else if (config.LEX_DS_MODE() == "none") {
    std::cout << "  Sample type: none" << std::endl;
    do_sample_traits = []() { /*do nothing by default*/; };

  } else {
    std::cout << "  Unknown sample type." << std::endl;
    exit(-1);
  }

}

void DiagWorld::SetOnOffspringReady()
{
  std::cout << "------------------------------------------------" << std::endl;
  std::cout << "Setting OnOffspringReady function..." << std::endl;

  OnOffspringReady(
    [this](Org & org, size_t parent_pos) {
      // quick checks
      emp_assert(fun_do_mutations);
      emp_assert(random_ptr);
      emp_assert(org.GetGenome().size() == config.OBJECTIVE_CNT());
      emp_assert(org.GetM() == config.OBJECTIVE_CNT());

      // do mutations on offspring
      size_t mcnt = fun_do_mutations(org, *random_ptr);

      // no mutations were applied to offspring
      if(mcnt == 0) {
        Org & parent = *pop[parent_pos];

        // quick checks
        emp_assert(parent.GetGenome().size() == config.OBJECTIVE_CNT());
        emp_assert(parent.GetM() == config.OBJECTIVE_CNT());

        // give everything to offspring from parent
        org.MeClone();
        org.Inherit(parent.GetScore(), parent.GetOptimal(), parent.GetCount(), parent.GetAggregate(), parent.GetStart(), parent.GetStreak());

      } else{
        org.Reset();
      }
    }
  );
  std::cout << "Finished setting OnOffspringReady function!\n" << std::endl;
}

void DiagWorld::SetEvaluation()
{
  std::cout << "------------------------------------------------" << std::endl;
  std::cout << "Setting Evaluation function..." << std::endl;

  target_t tar(config.OBJECTIVE_CNT(), config.TARGET());
  target.clear(); target.resize(config.OBJECTIVE_CNT());
  std::copy(tar.begin(), tar.end(), target.begin());

  diagnostic = emp::NewPtr<Diagnostic>(target, config.CREDIT());
  std::cout << "Created diagnostic emp::Ptr" << std::endl;

  test_score_matrix.resize(config.OBJECTIVE_CNT(), emp::vector<double>(config.POP_SIZE(), 0));

  if (config.DIAGNOSTIC() == "exploitation") {
    SetupDiagnostic_Exploitation();
  } else if (config.DIAGNOSTIC() == "struct-exploitation") {
    SetupDiagnostic_StructuredExploitation();
  } else if (config.DIAGNOSTIC() == "strong-ecology") {
    SetupDiagnostic_StrongEcology();
  } else if (config.DIAGNOSTIC() == "weak-ecology") {
    SetupDiagnostic_WeakEcology();
  } else if (config.DIAGNOSTIC() == "exploration") {
    SetupDiagnostic_Exploration();
  } else {
    std::cout << "ERROR: UNKNOWN DIAGNOSTIC, " << config.DIAGNOSTIC() << std::endl;
    exit(-1);
  }

  std::cout << "Evaluation function set!\n" <<std::endl;
}

void DiagWorld::SetDataTracking()
{
  std::cout << "------------------------------------------------" << std::endl;
  std::cout << "Setting up data tracking..." << std::endl;
  if (data_file != nullptr) {
    data_file.Delete();
  }
  if (elite_data_file != nullptr) {
    elite_data_file.Delete();
  }

  output_dir = config.OUTPUT_DIR();
  if (output_dir.back() != '/') {
    output_dir += "/";
  }
  mkdir(output_dir.c_str(), ACCESSPERMS);
  data_file = emp::NewPtr<emp::DataFile>(output_dir + "data.csv");

  // initialize all nodes
  std::cout << "Initializing nodes..." << std::endl;
  pop_fit.New();
  pop_opti.New();
  pnt_fit.New();
  pnt_opti.New();
  pop_str.New();
  std::cout << "Nodes initialized!" << std::endl;

  // track population aggregate score stats: average, variance, min, max
  data_file->AddMean(*pop_fit, "pop_fit_avg", "Population average aggregate performance.");
  data_file->AddVariance(*pop_fit, "pop_fit_var", "Population variance aggregate performance.");
  data_file->AddMax(*pop_fit, "pop_fit_max", "Population maximum aggregate performance.");
  data_file->AddMin(*pop_fit, "pop_fit_min", "Population minimum aggregate performance.");

  // track population optimized objective count stats: average, variance, min, max
  data_file->AddMean(*pop_opti, "pop_opt_avg", "Population average objective optimization count.");
  data_file->AddVariance(*pop_opti, "pop_opt_var", "Population variance objective optimization count.");
  data_file->AddMax(*pop_opti, "pop_opt_max", "Population maximum objective optimization count.");
  data_file->AddMin(*pop_opti, "pop_opt_min", "Population minimum objective optimization count.");

  // track parent aggregate score stats: average, variance, min, max
  data_file->AddMean(*pnt_fit, "pnt_fit_avg", "Parent average aggregate performance.");
  data_file->AddVariance(*pnt_fit, "pnt_fit_var", "Parent variance aggregate performance.");
  data_file->AddMax(*pnt_fit, "pnt_fit_max", "Parent maximum aggregate performance.");
  data_file->AddMin(*pnt_fit, "pnt_fit_min", "Parent minimum aggregate performance.");

  // track parent optimized objective count stats: average, variance, min, max
  data_file->AddMean(*pnt_opti, "pnt_opt_avg", "Parent average objective optimization count.");
  data_file->AddVariance(*pnt_opti, "pnt_opt_var", "Parent variance objective optimization count.");
  data_file->AddMax(*pnt_opti, "pnt_opt_max", "Parent maximum objective optimization count.");
  data_file->AddMin(*pnt_opti, "pnt_opt_min", "Parent minimum objective optimization count.");

  // track parent optimized objective count stats: average, variance, min, max
  data_file->AddMean(*pop_str, "pop_str_avg", "Population average streak count.");
  data_file->AddVariance(*pop_str, "pop_str_var", "Population variance streak count.");
  data_file->AddMax(*pop_str, "pop_str_max", "Population maximum streak count.");
  data_file->AddMin(*pop_str, "pop_str_min", "Population minimum streak count.");

  std::cout << "Added all data nodes to data file!" << std::endl;

  // update we are at
  data_file->AddFun<size_t>(
    [this]() {
      return update;
    },
    "gen",
    "Current generation at!"
  );

  // unique optimized objectives count
  data_file->AddFun<size_t>([this]()
  {
    return UniqueObjective();
  }, "pop_uni_obj", "Number of unique optimized traits per generation!");

  // elite solution aggregate performance
  data_file->AddFun<double>([this]()
  {
    // quick checks
    emp_assert(elite_pos != config.POP_SIZE());
    emp_assert(pop.size() == config.POP_SIZE());

    Org & org = *pop[elite_pos];

    return org.GetAggregate();
  }, "ele_agg_per", "Elite solution aggregate performance!");

  // elite solution optimized objectives count
  data_file->AddFun<size_t>([this]()
  {
    // quick checks
    emp_assert(elite_pos != config.POP_SIZE());
    emp_assert(pop.size() == config.POP_SIZE());

    Org & org = *pop[elite_pos];

    return org.GetCount();
  }, "ele_opt_cnt", "Elite solution optimized objective count!");

  // optimized solution aggregate performance
  data_file->AddFun<double>([this]()
  {
    // quick checks
    emp_assert(opti_pos != config.POP_SIZE());
    emp_assert(pop.size() == config.POP_SIZE());

    Org & org = *pop[opti_pos];

    return org.GetAggregate();
  }, "opt_agg_per", "Otpimal solution aggregate performance");

  // optimized solution optimized objectives count
  data_file->AddFun<size_t>([this]()
  {
    // quick checks
    emp_assert(opti_pos != config.POP_SIZE());
    emp_assert(pop.size() == config.POP_SIZE());

    Org & org = *pop[opti_pos];

    return org.GetCount();
  }, "opt_obj_cnt", "Otpimal solution aggregate performance");

  // streak solution aggregate performance
  data_file->AddFun<double>([this]()
  {
    // quick checks
    emp_assert(strk_pos != config.POP_SIZE());
    emp_assert(pop.size() == config.POP_SIZE());

    Org & org = *pop[strk_pos];

    return org.GetAggregate();
  }, "str_agg_per", "Otpimal solution aggregate performance");

  // streak solution optimized objectives count
  data_file->AddFun<size_t>([this]()
  {
    // quick checks
    emp_assert(strk_pos != config.POP_SIZE());
    emp_assert(pop.size() == config.POP_SIZE());

    Org & org = *pop[strk_pos];

    return org.GetCount();
  }, "str_obj_cnt", "Otpimal solution aggregate performance");

  // loss of diversity
  data_file->AddFun<double>([this]()
  {
    // quick checks
    emp_assert(parent_vec.size() == config.POP_SIZE());

    std::set<size_t> unique;
    for(auto & id : parent_vec) {unique.insert(id);}

    // ask Charles
    const double num = static_cast<double>(unique.size());
    const double dem = static_cast<double>(config.POP_SIZE());

    return num / dem;
  }, "los_div", "Loss in diversity generated by the selection scheme!");

  // selection pressure
  data_file->AddFun<double>([this]()
  {
    // quick checks
    emp_assert(pop_fit->GetCount() == config.POP_SIZE());
    emp_assert(pnt_fit->GetCount() == config.POP_SIZE());

    const double pop = pop_fit->GetMean();
    const double pnt = pnt_fit->GetMean();
    const double var = pop_fit->GetVariance();

    if(var == 0.0) {return 0.0;}

    return (pop - pnt) / var;
  }, "sel_pre", "Selection pressure applied by selection scheme!");

  // selection variance
  data_file->AddFun<double>([this]()
  {
    // quick checks
    emp_assert(pop_fit->GetCount() == config.POP_SIZE());
    emp_assert(pnt_fit->GetCount() == config.POP_SIZE());

    const double pop = pop_fit->GetVariance();
    const double pnt = pnt_fit->GetVariance();

    if(pnt == 0.0) {return 0.0;}

    return pop / pnt;
  }, "sel_var", "Selection pressure applied by selection scheme!");

  // unique starting positions
  data_file->AddFun<size_t>([this]()
  {
    emp_assert(pop_acti_gene.size()  == config.OBJECTIVE_CNT());
    return FindUniqueStart();
  }, "uni_str_pos", "Number of unique starting positions in the population!");

  // Pareto group count
  data_file->AddFun<size_t>([this]()
  {
    if(config.SELECTION() == "nondominated-sorting")
    {
      emp_assert(pareto_cnt != 0);
      return pareto_cnt;
    }
    else
    {
      emp_assert(pareto_cnt == 0);
      return pareto_cnt;
    }
  }, "pareto_cnt", "Number of Pareto groups generated!");

  // archive group count
  data_file->AddFun<size_t>([this]()
  {
    return archive.size();
  }, "archive_cnt", "Number phenotypes in the archive!");

  // archive group count
  data_file->AddFun<double>([this]()
  {
    return pmin;
  }, "pmin", "pmin used for archive approval!");

  // archive group count
  data_file->AddFun<double>([this]()
  {
    return arc_elite;
  }, "arc_elite", "archive best fitness found so far!");

  // archive unique optimal traits
  data_file->AddFun<size_t>([this]()
  {
    return std::accumulate(arc_opti_trt.begin(), arc_opti_trt.end(), 0);
  }, "arc_opti_trt", "Unique optimal traits found in the archive!");

  // archive unique activation genes
  data_file->AddFun<size_t>([this]()
  {
    return std::accumulate(arc_acti_gene.begin(), arc_acti_gene.end(), 0);
  }, "arc_acti_gene", "Unique activation genes found in the archive!");

  // archive unique activation genes
  data_file->AddFun<size_t>([this]()
  {
    return ActivationGeneOverlap();
  }, "overlap", "Unique activation genes found in the archive!");

  data_file->PrintHeaderKeys();

  // create elite csv plus headers
  elite_data_file = emp::NewPtr<emp::DataFile>(output_dir + "elite.csv");
  elite_data_file->AddFun<size_t>(
    [this]() {
      return update;
    },
    "gen",
    "Current generation"
  );
  elite_data_file->AddFun<std::string>(
    [this]() {
      std::ostringstream stream;
      const auto& elite_org = GetOrg(elite_pos);
      const auto& genome = elite_org.GetGenome();
      stream << "\"[";
      for (size_t i = 0; i < genome.size(); ++i) {
        if (i) stream << ",";
        stream << genome[i];
      }
      stream << "]\"";
      return stream.str();
    },
    "genome",
    "Genome of the current elite organism"
  );
  elite_data_file->AddFun<std::string>(
    [this]() {
      std::ostringstream stream;
      const auto& elite_org = GetOrg(elite_pos);
      const auto& phenotype = elite_org.GetScore();
      stream << "\"[";
      for (size_t i = 0; i < phenotype.size(); ++i) {
        if (i) stream << ",";
        stream << phenotype[i];
      }
      stream << "]\"";
      return stream.str();
    },
    "phenotype",
    "Phenotype of the current elite organism"
  );
  elite_data_file->AddFun<double>(
    [this]() {
      return GetOrg(elite_pos).GetAggregate();
    },
    "aggregate_score",
    "Aggregate score of elite organism"
  );
  elite_data_file->AddFun<size_t>(
    [this]() {
      return GetOrg(elite_pos).GetStart();
    },
    "start_pos",
    "Activation gene for elite organism"
  );
  elite_data_file->PrintHeaderKeys();

  std::cout << "Finished setting data tracking!\n" << std::endl;
}

void DiagWorld::PopulateWorld()
{
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "Populating world with initial solutions..." << std::endl;

  // random starting organisms
  if(config.START())
  {
    for(size_t i = 0; i < config.POP_SIZE(); ++i)
    {
      genome_t g = emp::RandomDoubleVector(*random_ptr, config.OBJECTIVE_CNT(), config.LOWER_BND(), config.UPPER_BND());
      Inject(g,1);
    }
  }
  // same starting organisms
  else
  {
    Org org(config.OBJECTIVE_CNT());
    Inject(org.GetGenome(), config.POP_SIZE());
  }

  std::cout << "Initialing world complete!" << std::endl;
}


///< principle steps during an evolutionary run

void DiagWorld::ResetData()
{
  // reset all data nodes
  pop_fit->Reset();
  pop_opti->Reset();
  pnt_fit->Reset();
  pnt_opti->Reset();
  pop_str->Reset();

  // reset all positon ids
  elite_pos = config.POP_SIZE();
  comm_pos = config.POP_SIZE();
  opti_pos = config.POP_SIZE();
  strk_pos = config.POP_SIZE();
  // unique_starts = config.POP_SIZE();
  pareto_cnt = 0;

  // reset all vectors/maps holding current gen data
  fit_vec.clear();
  parent_vec.clear();
  pop_acti_gene.clear();
}

void DiagWorld::DoEvaluationStep()
{
  // quick checks
  emp_assert(fit_vec.size() == 0);
  emp_assert(0 < pop.size());
  emp_assert(pop.size() == config.POP_SIZE());

  // iterate through the world and populate fitness vector
  fit_vec.resize(config.POP_SIZE());
  for(size_t org_id = 0; org_id < pop.size(); ++org_id) {
    Org& org = GetOrg(org_id);

    // no evaluate needed if offspring is a clone
    fit_vec[org_id] = (org.GetClone()) ? org.GetAggregate() : evaluate(org);

    // Fill out test matrix
    for (size_t test_id = 0; test_id < test_score_matrix.size(); ++test_id) {
      const emp::vector<double>& org_test_scores = org.GetScore();
      emp_assert(test_id < org_test_scores.size());
      test_score_matrix[test_id][org_id] = org_test_scores[test_id];
    }
  }
}

void DiagWorld::DoSelectionStep()
{
  // quick checks
  emp_assert(parent_vec.size() == 0);
  emp_assert(0 < pop.size());
  emp_assert(pop.size() == config.POP_SIZE());

  // store parents
  auto parents = select();
  emp_assert(parents.size() == config.POP_SIZE());

  parent_vec.resize(config.POP_SIZE());
  std::copy(parents.begin(), parents.end(), parent_vec.begin());
}

void DiagWorld::RecordData()
{
  /// Add data to all nodes

  // get pop data
  emp_assert(pop.size() == config.POP_SIZE());
  for(size_t i = 0; i < pop.size(); ++i)
  {
    const Org & org = *pop[i];
    pop_fit->Add(org.GetAggregate());
    pop_opti->Add(org.GetCount());
    pop_str->Add(org.GetStreak());

    const Org & par = *pop[parent_vec[i]];
    pnt_fit->Add(par.GetAggregate());
    pnt_opti->Add(par.GetCount());
  }
  emp_assert(pop_fit->GetCount() == config.POP_SIZE());
  emp_assert(pop_opti->GetCount() == config.POP_SIZE());
  emp_assert(pop_str->GetCount() == config.POP_SIZE());

  // get parent data
  emp_assert(parent_vec.size() == config.POP_SIZE());
  emp_assert(pnt_fit->GetCount() == config.POP_SIZE());
  emp_assert(pnt_opti->GetCount() == config.POP_SIZE());

  /// get all position ids

  FindEverything();

  /// fill vectors & map
  emp_assert(fit_vec.size() == config.POP_SIZE()); // should be set already
  emp_assert(parent_vec.size() == config.POP_SIZE()); // should be set already

  /// update the file
  data_file->Update();
  elite_data_file->Update();

  // output this so we know where we are in terms of generations and fitness
  Org & org = GetOrg(elite_pos);
  Org & opt = GetOrg(opti_pos);
  std::cout << "gen=" << GetUpdate() << ", max_fit=" << org.GetAggregate()  << ", max_opt=" << opt.GetCount() << std::endl;
}

void DiagWorld::DoReproductionStep()
{
  // quick checks
  emp_assert(parent_vec.size() == config.POP_SIZE());
  emp_assert(pop.size() == config.POP_SIZE());

  // go through parent ids and do births
  for(auto & id : parent_vec){
    DoBirth(GetGenomeAt(id), id);
  }
}


///< selection scheme implementations

void DiagWorld::SetupSelection_Truncation()
{
  std::cout << "Setting selection scheme: Truncation" << std::endl;

  // set select lambda to mu lambda selection
  select = [this]()
  {
    // quick checks
    emp_assert(selection); emp_assert(pop.size() == config.POP_SIZE());
    emp_assert(0 < pop.size()); emp_assert(fit_vec.size() == config.POP_SIZE());

    // group population by fitness
    fitgp_t group = selection->FitnessGroup(fit_vec);

    return selection->MLSelect(config.TRUNC_SIZE(), config.POP_SIZE(), group);
  };

  std::cout << "Truncation selection scheme set!" << std::endl;
}

void DiagWorld::SetupSelection_Tournament()
{
  std::cout << "Setting selection scheme: Tournament" << std::endl;

  select = [this]()
  {
    // quick checks
    emp_assert(selection); emp_assert(pop.size() == config.POP_SIZE());
    emp_assert(0 < pop.size()); emp_assert(fit_vec.size() == config.POP_SIZE());

    // will hold parent ids + get pop agg score values
    ids_t parent(pop.size());

    // get pop size amount of parents
    for(size_t i = 0; i < parent.size(); ++i)
    {
      parent[i] = selection->Tournament(config.TOUR_SIZE(), fit_vec);
    }

    return parent;
  };

  std::cout << "Tournament selection scheme set!" << std::endl;
}

void DiagWorld::SetupSelection_FitnessSharing()
{
  std::cout << "Setting selection scheme: FitnessSharing" << std::endl;
  std::cout << "Fitness Sharing applied on: ";
  if(!config.FIT_APPLI()){
    std::cout << "Genome" << std::endl;
  }
  else
  {
    std::cout << "Phenotype" << std::endl;
  }

  select = [this]()
  {
    // quick checks
    emp_assert(selection); emp_assert(pop.size() == config.POP_SIZE());
    emp_assert(0 < pop.size()); emp_assert(fit_vec.size() == config.POP_SIZE());

    // If we get asked to do standard tournament selection with fitness sharing enabled
    if(config.FIT_SIGMA() == 0.0)
    {
      // do stochastic remainder selection with unmodified fitness
      return selection->StochasticRemainder(fit_vec);
    }

    // are we using genomes or phenotypes for similarity comparison?
    gmatrix_t comps = (config.FIT_APPLI()) ? PopFitMat() : PopGenomes();

    // generate distance matrix + fitness transformation
    fmatrix_t dist_mat = selection->SimilarityMatrix(comps, config.PNORM_EXP());
    score_t tscore = selection->FitnessSharing(dist_mat, fit_vec, config.FIT_ALPHA(), config.FIT_SIGMA());

    return selection->StochasticRemainder(tscore);
  };

  std::cout << "Fitness sharing selection scheme set!" << std::endl;
}

void DiagWorld::SetupSelection_Lexicase() {
  std::cout << "Setting selection scheme: Lexicase" << std::endl;

  // note - these select functions could also be made faster by not returning a copy of the selected parents
  select = [this]() {
    emp_assert(selection);
    emp_assert(pop.size() == config.POP_SIZE());
    emp_assert(0 < pop.size());

    fmatrix_t matrix = PopFitMat(); // NOTE - this could be made faster

    do_sample_traits();
    return selection->Lexicase(
      matrix,
      pop.size(),
      sampled_trait_ids
    );

  };

}

void DiagWorld::SetupSelection_LexicaseEvenLead() {
  std::cout << "Setting selection scheme: LexicaseEvenLead" << std::endl;

  // note - these select functions could also be made faster by not returning a copy of the selected parents
  select = [this]() {
    emp_assert(selection);
    emp_assert(pop.size() == config.POP_SIZE());
    emp_assert(0 < pop.size());

    fmatrix_t matrix = PopFitMat(); // NOTE - this could be made faster

    do_sample_traits();
    return selection->LexicaseEvenLead(
      matrix,
      pop.size(),
      sampled_trait_ids
    );

  };

}

void DiagWorld::SetupSelection_EpsilonLexicase()
{
  std::cout << "Setting selection scheme: EpsilonLexicase" << std::endl;
  std::cout << "Epsilon: " << config.LEX_EPS() << std::endl;

  select = [this]()
  {
    // quick checks
    emp_assert(selection);
    emp_assert(pop.size() == config.POP_SIZE());
    emp_assert(0 < pop.size());

    fmatrix_t matrix = PopFitMat(); // NOTE - this could be made faster if we didn't copy this matrix and instead used a reference

    do_sample_traits();
    return selection->EpsilonLexicase(
      matrix,
      config.LEX_EPS(),
      pop.size(),
      sampled_trait_ids
    );
  };

  std::cout << "Epsilon Lexicase selection scheme set!" << std::endl;
}

void DiagWorld::SetupSelection_NonDominatedSorting()
{
  std::cout << "Setting selection scheme: NonDominatedSorting" << std::endl;

  select = [this]()
  {
    // quick checks
    emp_assert(selection); emp_assert(pop.size() == config.POP_SIZE());
    emp_assert(0 < pop.size()); emp_assert(pareto_cnt == 0);

    // select parent ids
    ids_t parent(pop.size());

    // get Pareto groups with ids
    fmatrix_t matrix = PopFitMat();
    pareto_t pgroups = selection->ParetoFrontGroups(matrix);

    // construct fitnesses
    // ParetoFitness
    score_t fitess = selection->ParetoFitness(pgroups, matrix, config.NDS_ALP(), config.NDS_SIG(), config.NDS_RED(), config.NDS_MAX());

    // track data
    pareto_cnt = pgroups.size();

    return selection->StochasticRemainder(fitess);
  };

  std::cout << "NonDominated Sorting selection scheme set!" << std::endl;
}

void DiagWorld::SetupSelection_NoveltySearch()
{
  std::cout << "Setting selection scheme: NoveltySearch" << std::endl;
  std::cout << "Starting PMIN: " << config.NOVEL_PMIN() << std::endl;
  // save starting pmin
  pmin = config.NOVEL_PMIN();

  // initialize vectors for achive
  arc_opti_trt = optimal_t(config.OBJECTIVE_CNT(), false);
  arc_acti_gene = optimal_t(config.OBJECTIVE_CNT(), false);
  std::cout << "Created vectors for tracking archive data" << std::endl;

  select = [this]()
  {
    // quick checks
    emp_assert(selection); emp_assert(pop.size() == config.POP_SIZE());
    emp_assert(0 < pop.size()); emp_assert(fit_vec.size() == config.POP_SIZE());

    // If we get asked to do standard tournament selection with novelty euclidean selected
    if(config.NOVEL_K() == 0)
    {
      // select parent ids
      ids_t parent(pop.size());

      for(size_t i = 0; i < parent.size(); ++i)
      {
        parent[i] = selection->Tournament(config.TOUR_SIZE(), fit_vec);
      }

      return parent;
    }

    // find nearest neighbors for each solution
    // NOVEL_K ammount of summation of squares score (x-y)^2
    fmatrix_t fit_mat = PopFitMat();
    fmatrix_t neighborhood = selection->NoveltySearchNearSum(fit_mat, config.NOVEL_K(), config.POP_SIZE(), config.PNORM_EXP());

    // calculate novelty scores
    score_t tscore = selection->NoveltySOS(neighborhood, config.NOVEL_K(), config.PNORM_EXP());

    // check if we need to reduce pmin
    if(archive_gens == config.NOVEL_GEN())
    {
      archive_gens = 0;
      pmin -= pmin * config.NOVEL_DOWN();
    }

    // archive tracking and updated gen count
    if(!ArchiveUpdate(tscore, fit_mat))
    {
      archive_gens++;
    }

    // select parent ids
    ids_t parent(pop.size());

    for(size_t i = 0; i < parent.size(); ++i)
    {
      parent[i] = selection->Tournament(config.TOUR_SIZE(), tscore);
    }

    return parent;
  };

  std::cout << "Novelty Search - Archive selection scheme set!" << std::endl;
}


///< evaluation function implementations

void DiagWorld::SetupDiagnostic_Exploitation()
{
  std::cout << "Setting exploitation diagnostic..." << std::endl;

  evaluate = [this](Org & org)
  {
    // set score & aggregate
    score_t score = diagnostic->Exploitation(org.GetGenome());
    org.SetScore(score);
    org.AggregateScore();

    // set the starting position
    org.StartPosition();

    // set optimal vector and count
    optimal_t opti = diagnostic->OptimizedVector(org.GetGenome(), config.ACCURACY());
    org.SetOptimal(opti);
    org.CountOptimized();

    // set streak
    org.CalcStreak();

    return org.GetAggregate();
  };

  std::cout << "Exploitation diagnotic set!" << std::endl;
}

void DiagWorld::SetupDiagnostic_StructuredExploitation()
{
  std::cout << "Setting structured exploitation diagnostic..." << std::endl;

  evaluate = [this](Org & org)
  {
    // set score & aggregate
    score_t score = diagnostic->StructExploitation(org.GetGenome());
    org.SetScore(score);
    org.AggregateScore();

    // set the starting position
    org.StartPosition();

    // set optimal vector and count
    optimal_t opti = diagnostic->OptimizedVector(org.GetGenome(), config.ACCURACY());
    org.SetOptimal(opti);
    org.CountOptimized();

    // set streak
    org.CalcStreak();

    return org.GetAggregate();
  };

  std::cout << "Structured exploitation diagnotic set!" << std::endl;
}

void DiagWorld::SetupDiagnostic_StrongEcology()
{
  std::cout << "Setting strong ecology diagnostic..." << std::endl;

  evaluate = [this](Org & org)
  {
    // set score & aggregate
    score_t score = diagnostic->StrongEcology(org.GetGenome());
    org.SetScore(score);
    org.AggregateScore();

    // set the starting position
    org.StartPosition();

    // set optimal vector and count
    optimal_t opti = diagnostic->OptimizedVector(org.GetGenome(), config.ACCURACY());
    org.SetOptimal(opti);
    org.CountOptimized();

    // set streak
    org.CalcStreak();

    return org.GetAggregate();
  };

  std::cout << "Strong ecology diagnotic set!" << std::endl;
}

void DiagWorld::SetupDiagnostic_Exploration()
{
  std::cout << "Setting exploration diagnostic..." << std::endl;

  evaluate = [this](Org & org)
  {
    // set score & aggregate
    score_t score = diagnostic->Exploration(org.GetGenome());
    org.SetScore(score);
    org.AggregateScore();

    // set the starting position
    org.StartPosition();

    // set optimal vector and count
    optimal_t opti = diagnostic->OptimizedVector(org.GetGenome(), config.ACCURACY());
    org.SetOptimal(opti);
    org.CountOptimized();

    // set streak
    org.CalcStreak();

    return org.GetAggregate();
  };

  std::cout << "Exploration diagnotic set!" << std::endl;
}

void DiagWorld::SetupDiagnostic_WeakEcology()
{
  std::cout << "Setting weak ecology diagnostic..." << std::endl;

  evaluate = [this](Org & org)
  {
    // set score & aggregate
    score_t score = diagnostic->WeakEcology(org.GetGenome());
    org.SetScore(score);
    org.AggregateScore();

    // set the starting position
    org.StartPosition();

    // set optimal vector and count
    optimal_t opti = diagnostic->OptimizedVector(org.GetGenome(), config.ACCURACY());
    org.SetOptimal(opti);
    org.CountOptimized();

    // set streak
    org.CalcStreak();

    return org.GetAggregate();
  };

  std::cout << "Weak ecology diagnotic set!" << std::endl;
}

///< data tracking

size_t DiagWorld::UniqueObjective()
{
  // quick checks
  emp_assert(0 < pop.size()); emp_assert(pop.size() == config.POP_SIZE());

  optimal_t unique;

  // novelty search unqiue objective trait count
  if(config.SELECTION() == "novelty")
  {
    emp_assert(arc_opti_trt.size() == config.OBJECTIVE_CNT());
    unique = arc_opti_trt;
  }
  else{ unique = optimal_t(config.OBJECTIVE_CNT(), false); }


  for(size_t o = 0; o < config.OBJECTIVE_CNT(); ++o)
  {
    // iterate pop to check is a solution has the objective optimized
    for(size_t p = 0; p < pop.size(); ++p)
    {
      Org & org = *pop[p];

      // quick checks
      emp_assert(org.GetOptimal().size() == config.OBJECTIVE_CNT());

      if(org.OptimizedAt(o))
      {
        unique[o] = true;
        break;
      }
    }
  }

  return std::accumulate(unique.begin(), unique.end(), 0);
}

size_t DiagWorld::FindUniqueStart()
{
  // quick checks
  emp_assert(0 < pop.size()); emp_assert(pop.size() == config.POP_SIZE());\
  emp_assert(pop_acti_gene.size() == config.OBJECTIVE_CNT());
  return std::accumulate(pop_acti_gene.begin(), pop_acti_gene.end(), 0);
}

void DiagWorld::FindEverything()
{
  // quick checks
  emp_assert(fit_vec.size() == config.POP_SIZE());
  emp_assert(pop_acti_gene.size() == 0);
  emp_assert(elite_pos == config.POP_SIZE());

  // bools to make sure got everything
  bool elite_b = false,  opti_b = false, strk_b = false;

  // collect number of unique starting positions
  pop_acti_gene = optimal_t(config.OBJECTIVE_CNT(), false);

  // loop and get data
  for(size_t i = 0; i < pop.size(); ++i)
  {
    const Org & org = *pop[i];

    // update the population unique activation genes
    pop_acti_gene[org.GetStart()] = true;

    // check if we need to do anything below
    if(elite_b && opti_b && strk_b) {continue;}

    // find max fit solution
    if(org.GetAggregate() == pop_fit->GetMax() && !elite_b) {elite_b = true; elite_pos = i;}
    // find max optimal count solution
    if(org.GetCount() == pop_opti->GetMax() && !opti_b) {opti_b = true; opti_pos = i;}
    // find max fit solution
    if(org.GetStreak() == pop_str->GetMax() && !strk_b) {strk_b = true; strk_pos = i;}
  }
}

size_t DiagWorld::ActivationGeneOverlap()
{
  // quick checks
  emp_assert(pop_acti_gene.size() == config.OBJECTIVE_CNT());
  size_t count = 0;

  // if novelty selection is running calculate overlap
  if(config.SELECTION() == "novelty")
  {
    emp_assert(arc_acti_gene.size() == config.OBJECTIVE_CNT());
    for(size_t i = 0; i < config.OBJECTIVE_CNT(); ++i)
    {
      if(pop_acti_gene[i] && arc_acti_gene[i]) {count++;}
    }
    return count;
  }
  else
  {
    emp_assert(arc_acti_gene.size() == 0);
    return count;
  }
}


///< helper functions

DiagWorld::fmatrix_t DiagWorld::PopFitMat()
{
  // quick checks
  emp_assert(pop.size() == config.POP_SIZE());

  // create matrix of population score vectors
  fmatrix_t matrix(pop.size());

  for(size_t i = 0; i < pop.size(); ++i)
  {
    Org & org = *pop[i];
    emp_assert(org.GetScore().size() == config.OBJECTIVE_CNT());

    // charles ask if this is the actual org score vector or a deep copy made
    matrix[i] = org.GetScore();
  }

  return matrix;
}

DiagWorld::gmatrix_t DiagWorld::PopGenomes()
{
  // quick checks
  emp_assert(pop.size() == config.POP_SIZE());

  gmatrix_t matrix(pop.size());

  for(size_t i = 0; i < pop.size(); ++i)
  {
    Org & org = *pop[i];
    emp_assert(org.GetGenome().size() == config.OBJECTIVE_CNT());

    // charles ask if this is the actual org genome or a deep copy made
    matrix[i] = org.GetGenome();
  }

  return matrix;
}

bool DiagWorld::ArchiveUpdate(const score_t & score, const fmatrix_t & dmat)
{
  //quick checks
  emp_assert(0 < score.size()); emp_assert(0 < dmat.size());
  emp_assert(score.size() == dmat.size());

  // archive insertion count
  size_t insert = 0;

  // check each solution novelty score
  for(size_t i = 0; i < score.size(); ++i)
  {
    // insert solution if lucky
    if(random_ptr->P(config.NOVEL_RI()))
    {
      // add phenotype to archive
      archive.push_back(dmat[i]);

      // update archive stats with solution data (if possible)
      ArchiveDataUpdate(i);

      if(config.NOVEL_CAP() < archive.size() && config.NOVEL_CQS()) {archive.pop_front();}
    }
    else if(score[i] > pmin)
    {
      // increment insertion counts for update
      ++insert;

      // add phenotype to the archive
      archive.push_back(dmat[i]);

      // update archive stats with solution data (if possible)
      ArchiveDataUpdate(i);

      if(config.NOVEL_CAP() < archive.size() && config.NOVEL_CQS()) {archive.pop_front();}
    }
  }

  if(4 < insert)
  {
    pmin += pmin * config.NOVEL_UP();
  }

  return 0 < insert;
}

void DiagWorld::ArchiveDataUpdate(const size_t org_id)
{
  // quick checks
  emp_assert(pop.size() == config.POP_SIZE());
  emp_assert(0 <= org_id); emp_assert(org_id < pop.size());
  emp_assert(config.SELECTION() == "novelty");

  // get org from
  Org & org = *pop[org_id];

  // update archive data if needed

  // archive current trait aggregate maximum
  if(arc_elite < org.GetAggregate()) {arc_elite = org.GetAggregate();}
  // update the archive activation gene vector
  arc_acti_gene[org.GetStart()] = true;
  // update the archive optimal trait vector
  for(size_t o = 0; o < config.OBJECTIVE_CNT(); ++o)
  {
    if(org.OptimizedAt(o))
      {
        arc_opti_trt[o] = true;
      }
  }
}

} // namespace diag

#endif