//
// Created by giacomo on 09/02/20.
//

#ifndef HIERARCHY_TESTS_TestingGraphBasic3_H
#define HIERARCHY_TESTS_TestingGraphBasic3_H

#include <tests/TestingGraph.h>
#include <concept_vector/ConceptVector.h>
#include <math_utils.h>

class TestingGraphBasic3  : public TestingGraph<std::vector<std::vector<double>>> {
    naryTree t{0};
    std::unordered_map<size_t, std::vector<std::vector<double>>> memoizationMap;
    size_t nodes = 0;
    double beta, alpha;

public:
    TestingGraphBasic3(Graph &ref, double beta, double alpha);

protected:
    std::vector<std::vector<double>> getVectorRepresentation(const size_t &current) const override {
        return memoizationMap.at(current);
    }

    double
    similarity(const std::vector<std::vector<double>> &lhs, const std::vector<std::vector<double>> &rhs) const override {
        return max_cosine_similarity(lhs, rhs);
    }

    void passGraphDataIfRequired(const Graph &graph) override {
        size_t id = 0;
        nodes = 0;
        std::unordered_map<size_t, size_t> treeIdToLocalTId;
        std::map<std::string, std::vector<double>> keyValueMap;
        for (const auto& cp : this->treeIdToPathString) {
            size_t tmp = t.addChild(string_split_to_sizetvector(cp.second), id);
            treeIdToLocalTId[cp.first] = id;
            nodes += tmp;
        }
        keyValueMap.insert(std::make_pair("", ConceptVector::multiple_descent_concept_problem_vector(t, std::vector<size_t>{}, beta, alpha)));
        for (const auto& cp : this->treeIdToPathString) {
            auto x = string_split_to_sizetvector(cp.second);
            keyValueMap.insert(std::make_pair(size_vector_to_string(x), ConceptVector::multiple_descent_concept_problem_vector(t, x, beta, alpha)));
        }
        for (auto it : this->morphismInv) {
            std::vector<std::vector<double>> result;
            for (const size_t& x : it.second) {
                result.emplace_back(keyValueMap[this->treeIdToPathString[it.first]]);
            }
            memoizationMap.emplace(it.first, result);
        }
    }
};


#endif //HIERARCHY_TESTS_TestingGraphBasic3_H
