//
// Created by giacomo on 09/02/20.
//

#ifndef HIERARCHY_TESTS_TestingGraphBasic2_H
#define HIERARCHY_TESTS_TestingGraphBasic2_H

#include <tests/TestingGraph.h>
#include <concept_vector/ConceptVector.h>
#include <math_utils.h>
#include "TestingGraphBasic1.h"

class TestingGraphBasic2  : public TestingGraph<std::vector<std::vector<double>>> {
    naryTree t{0};
    size_t nodes = 0;
    double beta;
    std::unordered_map<size_t, std::vector<std::vector<double>>> memoizationMap;

public:
    TestingGraphBasic2(Graph &ref, double beta);

protected:
    std::vector<std::vector<double>> getVectorRepresentation(const size_t &current) const override;
    double similarity(const std::vector<std::vector<double>> &lhs, const std::vector<std::vector<double>> &rhs) const override;
    void passGraphDataIfRequired(const Graph &graph) override;
};


#endif //HIERARCHY_TESTS_TestingGraphBasic2_H
