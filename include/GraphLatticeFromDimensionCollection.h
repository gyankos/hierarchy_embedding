//
// Created by giacomo on 19/05/2020.
//

#ifndef HIERARCHY_TESTS_GRAPHLATTICEFROMDIMENSIONCOLLECTION_H
#define HIERARCHY_TESTS_GRAPHLATTICEFROMDIMENSIONCOLLECTION_H

#include <vector>
#include <unordered_map>
#include <fstream>
#include <Graph.h>

class GraphLatticeFromDimensionCollection {
    std::unordered_map<size_t, size_t> nodeMap;
    std::unordered_map<size_t, std::string> labels;

    void nested_loop_join(const std::vector<size_t> &embedding_id, std::vector<size_t> vector_pos, std::vector<Graph> &graph_collections);
    std::vector<size_t> generateNode(size_t v) const;

public:
    std::vector<size_t> embedding_id;
    void generate(const std::string& filename, std::vector<Graph>& graph_collections);
};


#endif //HIERARCHY_TESTS_GRAPHLATTICEFROMDIMENSIONCOLLECTION_H
