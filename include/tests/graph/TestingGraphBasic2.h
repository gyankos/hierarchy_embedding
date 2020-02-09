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
    std::unordered_map<size_t, size_t> treeIdToLocalTId;
    std::map<std::string, std::vector<double>> keyValueMap;
    size_t nodes = 0;
    double beta;

public:
    TestingGraphBasic2(Graph &ref, double beta);

protected:
    std::vector<std::vector<double>> getVectorRepresentation(const size_t &current) override {
        std::vector<std::vector<double>> result;
        for (const size_t& x : this->morphismInv[current]) {
            result.emplace_back(keyValueMap[this->treeIdToPathString[x]]);
        }
        return result;
    }

    double
    similarity(const std::vector<std::vector<double>> &lhs, const std::vector<std::vector<double>> &rhs) override {
        double sim = std::numeric_limits<double>::min();
        for (auto x : lhs) {
            for (auto y: lhs) {
                sim = std::max(sim, cosine_similarity(x,y));
            }
        }
        return sim;
    }

    void passGraphDataIfRequired(const Graph &graph) override {
        size_t id = 0;
        nodes = 0;
        keyValueMap.clear();
        for (const auto& cp : this->treeIdToPathString) {
            size_t tmp = t.addChild(string_split_to_sizetvector(cp.second), id);
            treeIdToLocalTId[cp.first] = id;
            nodes += tmp;
        }
        keyValueMap.insert(std::make_pair("", ConceptVector::local_density_problem_vector(t, std::vector<size_t>{}, beta)));
        for (const auto& cp : this->treeIdToPathString) {
            auto x = string_split_to_sizetvector(cp.second);
            keyValueMap.insert(std::make_pair(size_vector_to_string(x), ConceptVector::local_density_problem_vector(t, x, beta)));
        }
    }
};


#endif //HIERARCHY_TESTS_TestingGraphBasic2_H
