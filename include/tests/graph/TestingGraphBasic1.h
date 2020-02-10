//
// Created by giacomo on 09/02/20.
//

#ifndef HIERARCHY_TESTS_TESTINGGRAPHBASIC1_H
#define HIERARCHY_TESTS_TESTINGGRAPHBASIC1_H

#include <tests/TestingGraph.h>
#include <concept_vector/ConceptVector.h>
#include <math_utils.h>

class TestingGraphBasic1  : public TestingGraph<std::vector<std::vector<double>>> {
    naryTree t{0};
    std::unordered_map<size_t, std::vector<std::vector<double>>> memoizationMap;
    size_t nodes = 0;


public:
    TestingGraphBasic1(Graph &ref);

protected:
    std::vector<std::vector<double>> getVectorRepresentation(const size_t &current) override;
    double
    similarity(const std::vector<std::vector<double>> &lhs, const std::vector<std::vector<double>> &rhs) override;
    void passGraphDataIfRequired(const Graph &graph) override;
};


#endif //HIERARCHY_TESTS_TESTINGGRAPHBASIC1_H
