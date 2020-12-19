//
// Created by giacomo on 19/12/20.
//

#ifndef HIERARCHY_TESTS_EXPORTGRAPHEMBEDDINGS_H
#define HIERARCHY_TESTS_EXPORTGRAPHEMBEDDINGS_H

#include <Graph.h>

class ExportGraphEmbeddings {
    Graph graph;
    size_t count;
    std::unordered_map<size_t, size_t> map;
    std::map<size_t, std::vector<std::vector<size_t>>> pathEmbeddings;

public:

    std::vector<std::vector<size_t>> getPathRepresentation(const std::string& node_label);
    void prepareForNewHierarchy();
    void addHierarchyEdge(const std::string& child_label, const std::string& parent_label, double score = 1.0);
    void finalizeForEmbeddingGeneration(double distanceFactor = 3, double decayFactor = 2);

};


#endif //HIERARCHY_TESTS_EXPORTGRAPHEMBEDDINGS_H
