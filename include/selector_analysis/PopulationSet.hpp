#pragma once

// Standard includes
#include <map>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <algorithm>

// Empirical includes
#include "emp/base/vector.hpp"
#include "emp/datastructs/set_utils.hpp"
#include "emp/tools/string_utils.hpp"

// Local includes
#include "csv-parser/parser.hpp"

namespace selector_analysis {

class Population {
public:
  using pop_info_t = std::unordered_map<std::string, std::string>;
protected:

  std::string name;
  pop_info_t info; ///< Info about how this population was generated
  emp::vector< emp::vector<double> > pop_test_case_scores;
  size_t num_test_cases;

public:

  Population(
    const std::string& pop_name,
    const pop_info_t& pop_info
  ) :
    name(pop_name),
    info(pop_info)
  { }

  size_t GetSize() const { return pop_test_case_scores.size(); }

  void SetNumTestCases(size_t n)  { num_test_cases = n; }
  size_t GetNumTestCases() const { return num_test_cases; }

  void UpdatePopSize(size_t n) {
    pop_test_case_scores.resize(n, emp::vector<double>(num_test_cases, 0.0));
  }

  const pop_info_t& GetPopInfo() const { return info; }

  void SetScores(size_t cand_id, const emp::vector<double>& test_scores) {
    if (cand_id >= pop_test_case_scores.size()) {
      pop_test_case_scores.resize(cand_id + 1);
    }
    pop_test_case_scores[cand_id] = test_scores;
  }
};

/// Loads population test score profiles from file.
/// Maps population id to population object, which stores loaded info about the particular population.
class PopulationSet {
public:
  using csv_parser_t = aria::csv::CsvParser;

protected:
  emp::vector<Population> pops;
  std::map<std::string, size_t> pop_id_map; ///< Maps loaded population name to population object in pops vector

  emp::vector<std::string> pop_desc_fields = emp::vector<std::string>({
    "ko_prob",
    "num_specialist_cases",
    "num_test_cases",
    "pass_prob",
    "pop_id",
    "pop_size",
    "pop_type",
    "specialists_per_case"
  });

public:

  size_t GetSize() const { return pops.size(); }

  /// Given path to file containing population test case scores (e.g., those generated by the gen_pop_profiles.py script),
  /// load all populations in file into this PopulationSet.
  /// If overlapping/conflicting pop ids (those being loaded vs those already loaded), skip conflicts (but complain).
  /// Required fields:
  /// TODO
  void LoadFromCSV(const std::string& path);

};

void PopulationSet::LoadFromCSV(const std::string& path) {
  std::ifstream pops_fstream(path);
  if (!pops_fstream.is_open()) {
    std::cout << "Failed to open population file (" << path << "). Exiting." << std::endl;
    exit(-1);
  }
  // if file is empty, failure!
  if (pops_fstream.eof()) return;

  std::string cur_line;
  emp::vector<std::string> line_components;

  // Collect header
  std::getline(pops_fstream, cur_line);
  emp::slice(cur_line, line_components, ',');
  std::unordered_map<std::string, size_t> header_lu;
  for (size_t i = 0; i < line_components.size(); ++i) {
    header_lu[line_components[i]] = i;
  }
  // Verify required header info
  // TODO
  // emp_assert(emp::Has(header_lu, *todo*));

  // Extract populations (individual-by-individual)
  csv_parser_t parser(pops_fstream);
  for (auto& row : parser) {
    auto& pop_name = row[header_lu["pop_id"]];
    // Extract population description
    std::unordered_map<std::string, std::string> pop_desc;
    for (auto& field : pop_desc_fields) {
      pop_desc[field] = row[header_lu[field]];
    }
    // Extract test case scores
    std::string test_scores_str = row[header_lu["test_scores"]];
    emp::remove_chars(test_scores_str, "[]");
    line_components.clear();
    emp::slice(test_scores_str, line_components, ',');
    emp::vector<double> test_scores(line_components.size(), 0.0);
    for (size_t i = 0; i < test_scores.size(); i++) {
      test_scores[i] = emp::from_string<double>(line_components[i]);
    }

    // If this is the first time we've seen this population, create a new population w/appropriate info.
    size_t pop_id = 0;
    if (!emp::Has(pop_id_map, pop_name)) {
      // Add pop name to id map.
      pop_id = pops.size();
      pop_id_map[pop_name] = pop_id;
      // Add pop to pops
      pops.emplace_back(
        pop_name,
        pop_desc
      );
      // Configure number of test cases
      pops.back().SetNumTestCases(emp::from_string<size_t>(pop_desc["num_test_cases"]));
      pops.back().UpdatePopSize(emp::from_string<size_t>(pop_desc["pop_size"]));
    }

    // Get population to add to
    emp_assert(emp::Has(pop_id_map, pop_name));
    pop_id = pop_id_map[pop_name];
    emp_assert(pop_id < pops.size());
    Population& population = pops[pop_id];
    emp_assert(population.GetPopInfo() == pop_desc);
    const size_t cand_id = emp::from_string<size_t>(row[header_lu["cand_id"]]);
    emp_assert(cand_id < population.GetSize(), cand_id, population.GetSize());
    emp_assert(test_scores.size() == population.GetNumTestCases());
    // Add test case profile to population
    population.SetScores(cand_id, test_scores);
  }
}

}