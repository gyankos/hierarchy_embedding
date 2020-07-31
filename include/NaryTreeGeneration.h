//
// Created by giacomo on 19/05/2020.
//

#ifndef HIERARCHY_TESTS_NARYTREEGENERATION_H
#define HIERARCHY_TESTS_NARYTREEGENERATION_H

#include <unordered_map>
#include <unordered_set>
#include <ostream>

struct NaryTreeGeneration {
    std::unordered_map<size_t, size_t> treeToGraphMorphism;
    std::unordered_map<size_t, std::unordered_set<size_t>> tree;
    std::map<size_t, std::string> treeIdToPathString;
    std::unordered_map<size_t, std::unordered_set<size_t>> morphismInv;
    size_t maximum_branching_factor;
    size_t number_of_vectors_required;

    NaryTreeGeneration() = default;
    NaryTreeGeneration(const NaryTreeGeneration& ) = default;
    NaryTreeGeneration& operator=(const NaryTreeGeneration& ) = default;

    friend std::ostream &operator<<(std::ostream &os, const NaryTreeGeneration &generation) {
        os << "maximum_branching_factor: " << generation.maximum_branching_factor << " number_of_vectors_required: "
           << generation.number_of_vectors_required;
        return os;
    }
};

#endif //HIERARCHY_TESTS_NARYTREEGENERATION_H
