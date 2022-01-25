#ifndef DIA_CONFIG_H
#define DIA_CONFIG_H

#include "config/config.h"

EMP_BUILD_CONFIG(DiaConfig,
  GROUP(WORLD, "How should the world be setup?"),
  VALUE(POP_SIZE,     size_t,      512,    "Population size."),
  VALUE(MAX_GENS,     size_t,    50000,    "Maximum number of generations."),
  VALUE(SEED,            int,        1,    "Random number seed."),
  VALUE(START,          bool,     true,    "Do we start randomly (true) or from the lowest point (false)"),

  GROUP(DIAGNOSTICS, "How are the diagnostics setup?"),
  VALUE(LOWER_BND,           double,       0.0,      "Lower bound for random starts."),
  VALUE(UPPER_BND,           double,       1.0,      "Upper bound for random starts."),
  VALUE(TARGET,              double,     100.0,      "Upper bound for random starts."),
  VALUE(ACCURACY,            double,      0.99,      "Accuracy percentage needed to be considered an optimal trait"),
  VALUE(CREDIT,              double,      0.00,      "Maximum credit a solution can get on an objective if applicable"),
  VALUE(OBJECTIVE_CNT,       size_t,       100,      "Number of traits an organism has"),
  VALUE(SELECTION,           size_t,         6,      "Which selection are we doing? \n0: Truncation\n1: Tournament\n2: Fitness Sharing\n"
                                                     "3: Novelty Aggregate\n4: Espilon Lexicase\n5: Novelty Euclidean\n6: Nondominated Sorting"
                                                     "\n7: Novelty Search"),
  VALUE(DIAGNOSTIC,          size_t,         0,      "Which diagnostic are we doing? \n0: Exploitation\n1: Structured Exploitation\n"
                                                     "2: Strong Ecology \n3: Exploration \n4: Weak Ecology"),

  GROUP(MUTATIONS, "Mutation rates for organisms."),
  VALUE(MUTATE_PER,       double,     0.007,        "Probability of instructions being mutated"),
  VALUE(MEAN,             double,     0.0,          "Mean of Gaussian Distribution for mutations"),
  VALUE(STD,              double,     1.0,          "Standard Deviation of Gaussian Distribution for mutations"),

  GROUP(PARAMETERS, "Parameter estimations all selection schemes."),
  VALUE(TRUNC_SIZE,       size_t,             8,       "Parameter estimate for μ."),
  VALUE(TOUR_SIZE,        size_t,             8,       "Parameter estimate for tournament size."),

  VALUE(FIT_SIGMA,        double,           0.0,       "Parameter estimate for proportion of similarity threshold sigma (based on maximum distance between solutions)."),
  VALUE(FIT_ALPHA,        double,           1.0,       "Parameter estimate for penalty function shape alpha."),
  VALUE(FIT_APPLI,        bool,           false,       "Fitness sharing applied: 0->Genome, 1->Phenotype"),

  VALUE(PNORM_EXP,        double,            2.0,       "Paramter we are using for the p-norm function."),
  VALUE(NOVEL_K,          size_t,             15,       "Parameter estimate k-nearest neighbors."),
  VALUE(NOVEL_PMIN,       double,            7.0,       "Minimum novelty score needed to enter archive."),
  VALUE(NOVEL_UP,         double,           0.25,       "Increase pmin by this much."),
  VALUE(NOVEL_DOWN,       double,           0.05,       "Decrease pmin by this much."),
  VALUE(NOVEL_RI,         double,        0.00001,       "Probability of random solution inserted in archive."),
  VALUE(NOVEL_GEN,        size_t,            500,       "Number of generations to lower pmin."),
  VALUE(NOVEL_CAP,        size_t,            512,       "Cap on number of solutions allowed in the archive"),
  VALUE(NOVEL_CQS,          bool,          false,       "Do we cap solutions?"),

  VALUE(LEX_EPS,          double,            0.0,       "Parameter estimate for lexicase epsilon."),

  VALUE(PAT_MAX,          double,   9000000000000000.0,       "Large dummy number ."),
  VALUE(PAT_DFT,          double,            1.0,       "Default large value for different Pareto groups."),
  VALUE(PAT_ALP,          double,            2.0,       "Alpha value for Pareto fitness sharing."),
  VALUE(PAT_SIG,          double,            0.1,       "Sigma value for Pareto fitness sharing."),

  GROUP(SYSTEMATICS, "Output rates for OpenWorld"),
  VALUE(SNAP_INTERVAL,             size_t,             10000,          "How many updates between prints?"),
  VALUE(PRINT_INTERVAL,            size_t,              1000,          "How many updates between prints?"),
  VALUE(OUTPUT_DIR,           std::string,              "./",          "What directory are we dumping all this data")
)

#endif