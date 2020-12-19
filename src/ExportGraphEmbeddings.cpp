//
// Created by giacomo on 19/12/20.
//

#include <tests/graph/TestingGraphProposal.h>
#include "ExportGraphEmbeddings.h"

void ExportGraphEmbeddings::prepareForNewHierarchy() {
    count = 0;
    graph.clear();
    map.clear();
}

void
ExportGraphEmbeddings::addHierarchyEdge(const std::string &child_label, const std::string &parent_label, double score) {
    graph.addEdgeExternally(child_label, parent_label, count, map, score);
}

void ExportGraphEmbeddings::finalizeForEmbeddingGeneration(double distanceFactor, double decayFactor) {
    // Initializing the embedding generator, directly from the same class as from the DAG paper
    TestingGraphProposal proposal(distanceFactor, decayFactor, graph);
    // Initializing the whole representation of the graph
    proposal.initializeEmbeddingGeneration(nullptr);
    for (const auto& cp : graph.nodeMap) {
        const size_t nodeId = cp.first;
        pathEmbeddings.emplace(nodeId, proposal.getVectorRepresentation(nodeId));
    }
    // clearing the memory with the unecessary elements
}

std::vector<std::vector<size_t>> ExportGraphEmbeddings::getPathRepresentation(const std::string &node_label) {
    size_t node_id = graph.getId(node_label);
    auto it = pathEmbeddings.find(node_id);
    if (it != pathEmbeddings.end()) {
        return it->second;
    } else {
        return {};
    }
}
