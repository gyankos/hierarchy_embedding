//
// Created by giacomo on 19/12/20.
//

#include <ExportGraphEmbeddings.h>

int main() {

    ///test_lattice();

    ExportGraphEmbeddings graph;
    graph.prepareForNewHierarchy();
    graph.addHierarchyEdge("B", "A");
    graph.addHierarchyEdge("C", "A");
    graph.addHierarchyEdge("D", "B");
    graph.addHierarchyEdge("D", "C");
    graph.finalizeForEmbeddingGeneration();

    std::cout << graph.getPathRepresentation("A") << std::endl;
    std::cout << graph.getPathRepresentation("B") << std::endl;
    std::cout << graph.getPathRepresentation("C") << std::endl;
    std::cout << graph.getPathRepresentation("D") << std::endl;

    return 0;
}
