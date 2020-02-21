//
// Created by giacomo on 30/12/19.
//

#include "tests/tree/TestingTreeExternal.h"

#include <fstream>
#include <sstream>
#include <string>
#include <cassert>
#include <iterator>

TestingTreeExternal::TestingTreeExternal(size_t maximumBranchingFactor, size_t maximumHeight, const std::string &filename)
: TestingTree(maximumBranchingFactor, maximumHeight) {
    std::ifstream infile(filename);
    std::string line;
    bool isFirstLine = true;
    size_t vectorDimension = 0;
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        if (isFirstLine) {
            iss >> vectorDimension;
            isFirstLine = false;
        } else {
            bool isFirstComponent = true;
            std::vector<double> vector;
            std::string dimensionName;

            std::vector<std::string> ret((std::istream_iterator<std::string>(iss)),
                                         std::istream_iterator<std::string>());

            assert(ret.size()-1 == vectorDimension);

            for (size_t i = 0, n = ret.size(); i<n; i++) {
                if (!i) {
                    dimensionName = ret[i];
                    if (dimensionName == "e")
                        dimensionName = "";
                } else {
                    vector.emplace_back(std::stod(ret[i]));
                }
            }
            memoization_map.insert(std::make_pair(dimensionName, vector));
        }
    }
}
