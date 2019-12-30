//
// Created by giacomo on 30/12/19.
//

#include "ExternalTesting.h"

#include <fstream>
#include <sstream>
#include <string>

ExternalTesting::ExternalTesting(size_t maximumBranchingFactor, size_t maximumHeight, const std::string &filename)
: Testing(maximumBranchingFactor, maximumHeight) {
    std::ifstream infile(filename);
}
