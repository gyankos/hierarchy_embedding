//
// Created by giacomo on 09/02/20.
//

#ifndef HIERARCHY_TESTS_TESTINGGRAPHPROPOSAL_H
#define HIERARCHY_TESTS_TESTINGGRAPHPROPOSAL_H


#include <cstdio>
#include <tests/TestingGraph.h>

class TestingGraphProposal : public TestingGraph<std::vector<std::vector<size_t>>> {
    Proposal* prop = nullptr;
    std::unordered_map<size_t, std::vector<std::vector<size_t>>> memoization_map;
    double distance;
    double decay;

public:
    TestingGraphProposal( double distanceFactor, double decayFactor, Graph &ref);
    ~TestingGraphProposal();

protected:
    std::vector<std::vector<size_t>> getVectorRepresentation(const size_t &current) const override;
    double similarity(const std::vector<std::vector<size_t>> &lhs, const std::vector<std::vector<size_t>> &rhs) const override;
    void passGraphDataIfRequired(const Graph &graph) override;
};


#endif //HIERARCHY_TESTS_TESTINGGRAPHPROPOSAL_H
