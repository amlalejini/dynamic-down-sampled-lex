#define CATCH_CONFIG_MAIN

// testing files
#include "/mnt/c/Users/josex/Desktop/Research/Repos/Catch/catch.hpp"
#include "../source/selection.h"

// empirical headers
#include "base/vector.h"
#include "tools/Random.h"
#include "tools/random_utils.h"

// library includes
#include <algorithm>
#include <string>
#include <cmath>

// const vars for test
constexpr size_t SEED = 17;

constexpr double PNORM_EXP = 2.0;
constexpr double POW_EXP = 0.5;

constexpr double SIMP_EXP = 2.0;
constexpr double SIMP_ERR = -1.0;

constexpr size_t NEAR_K = 3;


// In Tests directory, to run:
// clang++ -std=c++17 -I ../../../Empirical/source/ selection-test.cpp -o selection-test; ./selection-test

TEST_CASE("Selection class Pnorm function", "[initialization]")
{
  // all the vars we will be altering
  const size_t size = 5; double res; double exp;
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(SEED);
  emp::vector<double> x(size); emp::vector<double> y(size);
  Selection select(random);


  // first argument of std::pow calculated through wolfram alpha

  // just checking simple x vector computation
  x = {1.0,1.0,1.0,1.0,1.0};
  y = {0.0,0.0,0.0,0.0,0.0};
  res = select.Pnorm(x,y,PNORM_EXP);
  exp = std::pow(5.0, POW_EXP);
  REQUIRE(res == Approx(exp));
  // check that the reveserse also holds
  res = select.Pnorm(y,x,PNORM_EXP);
  REQUIRE(res == Approx(exp));


  // another bigger simple x vector computation
  x = {10.0,10.0,10.0,10.0,10.0};
  y = {0.0,0.0,0.0,0.0,0.0};
  res = select.Pnorm(x,y,PNORM_EXP);
  exp = std::pow(500.0, POW_EXP);
  REQUIRE(res == Approx(exp));
  // check that the reveserse also holds
  res = select.Pnorm(y,x,PNORM_EXP);
  REQUIRE(res == Approx(exp));


  // compute difference accuracy between vector x & y
  x = {10.0,10.0,10.0,10.0,10.0};
  y = {5.0,5.0,5.0,5.0,5.0};
  res = select.Pnorm(x,y,PNORM_EXP);
  exp = std::pow(125.0, POW_EXP);
  REQUIRE(res == Approx(exp));
  // check that the reveserse also holds
  res = select.Pnorm(y,x,PNORM_EXP);
  REQUIRE(res == Approx(exp));


  // check with small precision abnormalities
  x = {10.0,10.0,10.0,10.0,10.0};
  y = {5.5,5.5,5.5,5.5,5.5};
  res = select.Pnorm(x,y,PNORM_EXP);
  exp = std::pow(101.25, POW_EXP);
  REQUIRE(res == Approx(exp));
  // check that the reveserse also holds
  res = select.Pnorm(y,x,PNORM_EXP);
  REQUIRE(res == Approx(exp));

  random.Delete();
}

TEST_CASE("Selection class similarity matrix function","[similarity]")
{
  // all the vars we will be altering
  const size_t size = 5; const size_t pop = 4; double diff;
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(SEED);
  emp::vector<double> x(size); emp::vector<double> y(size);
  emp::vector<double> u(size); emp::vector<double> v(size);
  emp::vector<emp::vector<double>> tmat(size-1);
  emp::vector<emp::vector<double>> mat(size-1);
  emp::vector<emp::vector<double>> exp(size-1);
  Selection select(random);


  // create genomes
  x = {1,1,1,1,1}; y = {1,1,1,1,1}; u = {1,1,1,1,1}; v = {1,1,1,1,1};
  // create genome matrix
  mat = {x,y,u,v};
  // test function
  tmat = select.SimilarityMatrix(mat, SIMP_EXP);
  // create expected return (0's gained from premature calculation)
  exp = {{SIMP_ERR,SIMP_ERR,SIMP_ERR,SIMP_ERR},{0,SIMP_ERR,SIMP_ERR,SIMP_ERR},{0,0,SIMP_ERR,SIMP_ERR},{0,0,0,SIMP_ERR}};
  // level one checks
  REQUIRE(tmat.size() == exp.size());
  // level two checks
  for(size_t i = 0; i < pop; ++i){REQUIRE(tmat[i].size() == exp[i].size());}
  // level three checks
  for(size_t i = 0; i < pop; ++i){REQUIRE_THAT(tmat[i], Catch::Matchers::Equals(exp[i]));}


  // alter second genome
  // create genomes
  x = {1,1,1,1,1}; y = {0,0,0,0,0}; u = {1,1,1,1,1}; v = {1,1,1,1,1};
  // create genome matrix
  mat = {x,y,u,v};
  // test function
  tmat = select.SimilarityMatrix(mat, SIMP_EXP);
  // create expected return (0's gained from premature calculation)
  diff = std::pow(5.0, 0.5);
  exp = {{SIMP_ERR,SIMP_ERR,SIMP_ERR,SIMP_ERR},{diff,SIMP_ERR,SIMP_ERR,SIMP_ERR},{0,diff,SIMP_ERR,SIMP_ERR},{0,diff,0,SIMP_ERR}};
  // level one checks
  REQUIRE(tmat.size() == exp.size());
  // level two checks
  for(size_t i = 0; i < pop; ++i){REQUIRE(tmat[i].size() == exp[i].size());}
  // level three checks
  for(size_t i = 0; i < pop; ++i){REQUIRE_THAT(tmat[i], Catch::Matchers::Equals(exp[i]));}


  // alter third genome
  // create genomes
  x = {1,1,1,1,1}; y = {1,1,1,1,1}; u = {0,0,0,0,0}; v = {1,1,1,1,1};
  // create genome matrix
  mat = {x,y,u,v};
  // test function
  tmat = select.SimilarityMatrix(mat, SIMP_EXP);
  // create expected return (0's gained from premature calculation)
  diff = std::pow(5.0, 0.5);
  exp = {{SIMP_ERR,SIMP_ERR,SIMP_ERR,SIMP_ERR},{0,SIMP_ERR,SIMP_ERR,SIMP_ERR},{diff,diff,SIMP_ERR,SIMP_ERR},{0,0,diff,SIMP_ERR}};
  // level one checks
  REQUIRE(tmat.size() == exp.size());
  // level two checks
  for(size_t i = 0; i < pop; ++i){REQUIRE(tmat[i].size() == exp[i].size());}
  // level three checks
  for(size_t i = 0; i < pop; ++i){REQUIRE_THAT(tmat[i], Catch::Matchers::Equals(exp[i]));}


  // alter last genome
  // create genomes
  x = {1,1,1,1,1}; y = {1,1,1,1,1}; u = {1,1,1,1,1}; v = {0,0,0,0,0};
  // create genome matrix
  mat = {x,y,u,v};
  // test function
  tmat = select.SimilarityMatrix(mat, SIMP_EXP);
  // create expected return (0's gained from premature calculation)
  diff = std::pow(5.0, 0.5);
  exp = {{SIMP_ERR,SIMP_ERR,SIMP_ERR,SIMP_ERR},{0,SIMP_ERR,SIMP_ERR,SIMP_ERR},{0,0,SIMP_ERR,SIMP_ERR},{diff,diff,diff,SIMP_ERR}};
  // level one checks
  REQUIRE(tmat.size() == exp.size());
  // level two checks
  for(size_t i = 0; i < pop; ++i){REQUIRE(tmat[i].size() == exp[i].size());}
  // level three checks
  for(size_t i = 0; i < pop; ++i){REQUIRE_THAT(tmat[i], Catch::Matchers::Equals(exp[i]));}

  random.Delete();
}

TEST_CASE("Selection class fitness nearest neigbor function", "[neighbor]")
{
  // all the vars we will be altering
  const size_t pop = 10;
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(SEED);
  emp::vector<emp::vector<double>> neigh(pop);
  emp::vector<double> score(pop);
  Selection select(random);


  // strict ascending order
  score = {1,2,3,4,5,6,7,8,9,10};
  // get neighborhood scores
  neigh = select.FitNearestN(score, NEAR_K);
  // build correct vector of vecotr fitnesses
  




  random.Delete();
}



