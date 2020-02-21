//
// Created by giacomo on 16/02/20.
//

#ifndef HIERARCHY_TESTS_TESTINGGRAPHEXTERNAL_H
#define HIERARCHY_TESTS_TESTINGGRAPHEXTERNAL_H

#include <vector>
#include <tests/TestingGraph.h>

class TestingGraphExternal : public TestingGraph<std::vector<double>> {

    std::map<size_t, std::vector<double>> memoization_map;

public:
    TestingGraphExternal(Graph &ref, const std::string &filename);

protected:
    std::vector<double> getVectorRepresentation(const size_t &current) const override {
        return memoization_map.at(current);
    }

    double similarity(const std::vector<double> &lhs, const std::vector<double> &rhs) const override {
        double distance = euclideanDistance(lhs, rhs);
        double normalized = distance / (distance+1);
        return 1.0 - normalized;
    }

    void passGraphDataIfRequired(const Graph &graph) override {}
};


#endif //HIERARCHY_TESTS_TESTINGGRAPHEXTERNAL_H
