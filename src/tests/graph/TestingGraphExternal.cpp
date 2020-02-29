//
// Created by giacomo on 16/02/20.
//

#include "tests/graph/TestingGraphExternal.h"

TestingGraphExternal::TestingGraphExternal(Graph &ref,const std::string &filename) : TestingGraph(ref) {
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
            size_t dimensionId;

            std::vector<std::string> ret((std::istream_iterator<std::string>(iss)),
                                         std::istream_iterator<std::string>());

            assert(ret.size()-1 == vectorDimension);

            for (size_t i = 0, n = ret.size(); i<n; i++) {
                if (!i) {
                    dimensionId = ref.getId(ret[i]);
                } else {
                    vector.emplace_back(std::stod(ret[i]));
                }
            }
            memoization_map.insert(std::make_pair(dimensionId, vector));
        }
    }
}
