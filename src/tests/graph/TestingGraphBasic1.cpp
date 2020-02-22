//
// Created by giacomo on 09/02/20.
//

#include "tests/graph/TestingGraphBasic1.h"

TestingGraphBasic1::TestingGraphBasic1(Graph &ref) : TestingGraph(ref) {}

std::vector<std::vector<double>> TestingGraphBasic1::getVectorRepresentation(const size_t &current) const {
    return memoizationMap.at(current);
}

double
TestingGraphBasic1::similarity(const std::vector<std::vector<double>> &lhs, const std::vector<std::vector<double>> &rhs) const {
    return max_cosine_similarity(lhs, rhs);
}

void TestingGraphBasic1::passGraphDataIfRequired(const Graph &graph) {
    size_t id = 0;
    nodes = 0;
    std::unordered_map<std::string, std::vector<double>> keyValueMap;
    std::unordered_map<size_t, size_t> treeIdToLocalTId;
    for (const auto& cp : this->treeIdToPathString) {
        size_t tmp = t.addChild(string_split_to_sizetvector(cp.second), id);
        treeIdToLocalTId[cp.first] = id;
        nodes += tmp;
    }
    keyValueMap.insert(std::make_pair("", ConceptVector::relevancy_vector(t, std::vector<size_t>{})));
    for (const auto& cp : this->treeIdToPathString) {
        auto x = string_split_to_sizetvector(cp.second);
        keyValueMap.insert(std::make_pair(size_vector_to_string(x), ConceptVector::relevancy_vector(t, x)));
    }
    for (auto it : this->morphismInv) {
        std::vector<std::vector<double>> result;
        for (const size_t& x : it.second) {
            result.emplace_back(keyValueMap[this->treeIdToPathString[it.first]]);
        }
        memoizationMap.emplace(it.first, result);
    }
}
