#define CATCH_CONFIG_MAIN

#include "Catch/single_include/catch2/catch.hpp"

#include <unordered_set>

#include "emp/datastructs/vector_utils.hpp"

#include "selection/Lexicase.hpp"

TEST_CASE("Test selection::LexicaseSelect", "[selection][lexicase]") {
  using score_fun_t = std::function<double(void)>;

  constexpr size_t seed = 2;
  emp::Random random(seed);

  emp::vector< emp::vector<score_fun_t> > score_fun_sets;
  emp::vector< emp::vector<double> > scores{
    /* 0= */ {1.0, 1.0, 1.0},
    /* 1= */ {1.0, 0.0, 0.0},
    /* 2= */ {0.0, 1.0, 0.0},
    /* 3= */ {0.0, 0.0, 1.0}
  };

  for (size_t i = 0; i < scores.size(); ++i) {
    score_fun_sets.emplace_back();
    for (size_t j = 0; j < scores[i].size(); ++j) {
      score_fun_sets[i].emplace_back(
        [i,j,&scores](){return scores[i][j];}
      );
    }
  }

  selection::LexicaseSelect lex_selector(
    score_fun_sets,
    random
  );

  // With regular lexicase, should only ever select ID 0.
  for (size_t i = 0; i < 10; ++i) {
    lex_selector(100);
    const auto& selected = lex_selector.GetSelected();
    REQUIRE(selected.size() == 100);
    for (size_t id : selected) {
      REQUIRE(id == 0);
    }
  }

  // If we do not include ID 0 in candidates that can be considered, that ID should never show up
  for (size_t i = 0; i < 10; ++i) {
    lex_selector(
      100,
      {0, 1, 2},
      {1, 2, 3}
    );
    const auto& selected = lex_selector.GetSelected();
    REQUIRE(selected.size() == 100);
    for (size_t id : selected) {
      REQUIRE(id != 0);
    }
  }

  // Filtering down to just one allowed candidate should only return that candidate.
  for (size_t i = 0; i < 10; ++i) {
    lex_selector(
      100,
      {0, 1, 2},
      {2}
    );
    const auto& selected = lex_selector.GetSelected();
    REQUIRE(selected.size() == 100);
    for (size_t id : selected) {
      REQUIRE(id == 2);
    }
  }

  // Down-sampling to test 2 should always result in selection of candidate 3 (if we don't allow 0).
  for (size_t i = 0; i < 10; ++i) {
    lex_selector(
      100,
      {2},
      {1, 2, 3}
    );
    const auto& selected = lex_selector.GetSelected();
    REQUIRE(selected.size() == 100);
    for (size_t id : selected) {
      REQUIRE(id == 3);
    }
  }

}