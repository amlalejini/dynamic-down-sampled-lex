// This is the main function for the NATIVE version of this project.

#include <iostream>
#include <limits>

#include "emp/base/vector.hpp"
#include "emp/config/command_line.hpp"
#include "emp/config/ArgManager.hpp"

#include "selector_analysis/config.hpp"
#include "selector_analysis/PopulationSet.hpp"
#include "selector_analysis/SelectorAnalyzer.hpp"
// #include "diagnostics/world.hpp"
// #include "diagnostics/org.hpp"

int main(int argc, char* argv[])
{
  selector_analysis::SelectorAnalysisConfig config;
  config.Read("selector_analysis.cfg", false);
  auto args = emp::cl::ArgManager(argc, argv);
  if (args.ProcessConfigOptions(config, std::cout, "selector_analysis.cfg", "selector_analysis.hpp") == false) exit(0);
  if (args.TestUnknown() == false) exit(0);  // If there are leftover args, throw an error.

  std::cout << "==============================" << std::endl;
  std::cout << "|    How am I configured?    |" << std::endl;
  std::cout << "==============================" << std::endl;
  config.Write(std::cout);
  std::cout << "==============================\n"
            << std::endl;

  // todo run analysis
  selector_analysis::SelectorAnalyzer analyzer(config);
  analyzer.Run();

  // selector_analysis::PopulationSet pop_set;
  // pop_set.LoadFromCSV("pop.csv");
  // std::cout << "Pop set size: " << pop_set.GetSize() << std::endl;
}