// This is the main function for the NATIVE version of this project.

#include <iostream>
#include <limits>

#include "emp/base/vector.hpp"
#include "emp/config/command_line.hpp"
#include "emp/config/ArgManager.hpp"

#include "diagnostics/config.hpp"
#include "diagnostics/world.hpp"
#include "diagnostics/org.hpp"

// Hello world

int main(int argc, char* argv[])
{
  diag::DiaConfig config;
  config.Read("diagnostics.cfg", false);
  auto args = emp::cl::ArgManager(argc, argv);
  if (args.ProcessConfigOptions(config, std::cout, "diagnostics.cfg", "diagnostics-macros.h") == false) exit(0);
  if (args.TestUnknown() == false) exit(0);  // If there are leftover args, throw an error.

  std::cout << "==============================" << std::endl;
  std::cout << "|    How am I configured?    |" << std::endl;
  std::cout << "==============================" << std::endl;
  config.Write(std::cout);
  std::cout << "==============================\n"
            << std::endl;


  diag::DiagWorld world(config);

  // for (size_t ud = 0; ud <= config.MAX_GENS(); ud++)
  // {
  //   world.Update();
  // }
}