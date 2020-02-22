//
// Created by giacomo on 09/02/20.
//

#include "tests/graph/TestingGraphProposal.h"

TestingGraphProposal::TestingGraphProposal(double distanceFactor, double decayFactor, Graph &ref) : TestingGraph(ref), distance{distanceFactor}, decay{decayFactor} {}

TestingGraphProposal::~TestingGraphProposal() {
    delete prop;
}

std::vector<std::vector<size_t>> TestingGraphProposal::getVectorRepresentation(const size_t &current) const {
    return memoization_map.at(current); // concurrent access to the map
}

double TestingGraphProposal::similarity(const std::vector<std::vector<size_t>> &lhs,
                                        const std::vector<std::vector<size_t>> &rhs) const {
    double similarity = 0;
    for (const auto x : lhs) {
        for (const auto y : rhs) {
            double distanceMetric =  prop->rankingMetric(x, y);
            double normalizedDistance = (distanceMetric == std::numeric_limits<size_t>::max()) ? 1 : distanceMetric / (1+distanceMetric);
            double local_similairity = 1 - normalizedDistance;
            similarity = std::max(local_similairity, similarity);
        }
    }
    return similarity;
}

void TestingGraphProposal::passGraphDataIfRequired(const Graph &graph) {
    //Now, we can read the data from the graph!
    if (prop) {
        delete prop;
        prop = nullptr;
    }
    prop = new Proposal(graph.maxBranch, distance, decay);
    // memoizing the vector representations
    std::cout << "Taking some time to memoize the computations..." << std::endl;
    for (auto it : this->morphismInv) {
        std::vector<std::vector<size_t>> result;
        for (const size_t& x : it.second) {
            result.emplace_back(
                    string_split_to_sizetvector(this->treeIdToPathString[x]));
        }
        memoization_map.insert(std::make_pair(it.first, result));
    }
    std::cout << "...done! " << std::endl;
}
