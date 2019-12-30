//
// Created by giacomo on 30/12/19.
//

#ifndef HIERARCHY_TESTS_EXTERNALTESTING_H
#define HIERARCHY_TESTS_EXTERNALTESTING_H

#include <math_utils.h>
#include "Testing.h"

/**
 * Class function used to train external methods. In such methods, we assume that the string associated to the entity
 * is the same that is obtained via size_vector_to_string for representing a path vector as a
 */
class ExternalTesting : public Testing<std::string> {
    std::map<std::string, std::vector<double>> trainedExtenralVectors;

public:
    ExternalTesting(size_t maximumBranchingFactor, size_t maximumHeight, const std::string& filename);

protected:
    void initialize_hierarchy_with_all_paths(const std::vector<std::vector<size_t>> &subgraph_as_paths) override {
        std::cerr << "SKIPPING: the training happened in a external function" << std::endl;
    }

    std::string getVectorRepresentation(const std::vector<size_t> &current) override {
        return size_vector_to_string(current);
    }

    double similarity(const std::string &lhs, const std::string &rhs) override {
        return cosine_similarity(trainedExtenralVectors[lhs], trainedExtenralVectors[rhs]);
    }

    void generateTopKCandidates(PollMap<double, std::string> &map, const std::vector<size_t> &current) override {
        for (auto& it : trainedExtenralVectors) {
            map.add(similarity(it.first, size_vector_to_string(current)), it.first);
        }
    }

    void finalizeDataIngestion() override { }
};


#endif //HIERARCHY_TESTS_EXTERNALTESTING_H
