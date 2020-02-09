//
// Created by giacomo on 09/02/20.
//

#ifndef HIERARCHY_TESTS_TESTINGGRAPHPROPOSAL_H
#define HIERARCHY_TESTS_TESTINGGRAPHPROPOSAL_H


#include <cstdio>
#include <tests/TestingGraph.h>

class TestingGraphProposal : public TestingGraph<std::vector<std::vector<size_t>>> {

    Proposal* prop = nullptr;
    double distance;
    double decay;

public:
    TestingGraphProposal( double distanceFactor, double decayFactor, Graph &ref) : TestingGraph(ref), distance{distanceFactor}, decay{decayFactor} {}



protected:
    std::vector<std::vector<size_t>> getVectorRepresentation(const size_t &current) override {
        std::vector<std::vector<size_t>> result;
        for (const size_t& x : this->morphismInv[current]) {
            result.emplace_back(
                    string_split_to_sizetvector(this->treeIdToPathString[x]));
        }
        return result;
    }

    double similarity(const std::vector<std::vector<size_t>> &lhs, const std::vector<std::vector<size_t>> &rhs) override {
        double similarity = std::numeric_limits<double>::max();
        for (const auto x : lhs) {
            for (const auto y : rhs) {
                double distanceMetric =  prop->rankingMetric(x, y);
                double normalizedDistance = distanceMetric / (1+distanceMetric);
                double similairity = 1 - normalizedDistance;
            }
        }

        return similarity;
    }

    void generateTopKCandidates(PollMap<double, size_t> &map, const size_t &current) override {
        auto allVectors = getVectorRepresentation(current);
        for (auto & x : this->treeToGraphMorphism) {
            double  score = similarity(allVectors, getVectorRepresentation(x.second));
            map.add(score, x.second);
        }
    }

    void passGraphDataIfRequired(const Graph &graph) override {
        //Now, we can read the data from the graph!
        if (prop) {
            delete prop;
            prop = nullptr;
        }
        prop = new Proposal(graph.maxBranch, distance, decay);
    }

};


#endif //HIERARCHY_TESTS_TESTINGGRAPHPROPOSAL_H
