//
// Created by giacomo on 19/05/2020.
//

#include <Graph.h>
#include "GraphLatticeFromDimensionCollection.h"

void GraphLatticeFromDimensionCollection::nested_loop_join(const std::vector<size_t> &embedding_id,
                                                           std::vector<size_t> vector_pos,
                                                           std::vector<Graph> &graph_collections) {
    size_t vps = vector_pos.size();
    if (vps == graph_collections.size()-1) {
        Graph& g = graph_collections.back();
        for (const auto& cp : g.nodeMap) {
            std::vector<size_t> copy{vector_pos};
            copy.emplace_back(cp.first);
            size_t computed_id = 0;
            std::string key = "<";
            for (size_t i = 0, n = embedding_id.size(); i<n; i++) {
                computed_id += embedding_id[i] * copy[i];
                key += graph_collections[i].getName(copy[i]);
                if (i != n-1)
                    key += ',';
            }
            key += ">";
            nodeMap[computed_id] = computed_id;
            ///addNewNode(computed_id);
            labels[computed_id] = key;
            ///fileElementNameToId.put(key, computed_id);
            /*std::cerr <<
                      std::string(vps, '.') << cp.first << "<->" << g.fileElementNameToId.getKey(cp.first) << std::endl;
            std::cerr << std::string(vps+2, ' ') << key << "<->" << computed_id << "=" << copy << '=' << generateNode(embedding_id, computed_id) << std::endl;
        */}
    } else {
        Graph& g = graph_collections[vector_pos.size()];
        for (const auto& cp : g.nodeMap) {
            std::vector<size_t> copy{vector_pos};
            copy.emplace_back(cp.first);
            /*std::cerr <<
                      std::string(vps, '.') << cp.first << "<->" << g.fileElementNameToId.getKey(cp.first) << std::endl;*/
            nested_loop_join(embedding_id, copy, graph_collections);
        }
    }
}

std::vector<size_t> GraphLatticeFromDimensionCollection::generateNode(size_t v) const {
    std::vector<size_t> elements{embedding_id.size(), 0};
    for (size_t i = 0, n = embedding_id.size(); i<n; i++) {
        size_t operand = embedding_id[n-i-1];
        elements[n-i-1] = (size_t)v / operand;
        v = v % operand; /* Likely uses the result of the division. */
    }
    return elements;
}

void GraphLatticeFromDimensionCollection::generate(const std::string &filename, std::vector<Graph>& graph_collections) {
    std::ofstream file{filename};
    // Declaring the embedding for representing the lattice with a single number
    embedding_id.emplace_back(1);
    size_t embedding_id_current = 1;
    for (size_t i = 1, n = graph_collections.size(); i<n; i++) {
        embedding_id_current *= graph_collections[i-1].nodeSize();
        embedding_id.emplace_back(embedding_id_current);
    }

    // Now, generating the cartesian product over the vertices
    nested_loop_join(embedding_id, {}, graph_collections);

    //graph_collections[0].print_graph();
    //graph_collections[1].print_graph();
    //return;

    // Now, generating all the edges among the nodes via cartesian product
    for (const auto& cp1: nodeMap) {
        auto u = cp1.first;
        auto u_vector = generateNode(u); // Generating the original dimensions' nodes from the id

        for (const auto& cp2: nodeMap) {
            auto v = cp2.first;
            if (u == v) continue; // Avoiding creating a hook over a node

            auto v_vector = generateNode(v); // Generating the original dimensions' nodes form the id
            bool doSkip = false;

            // Checking if there is at most one node change
            int candidate_pos = -1;
            for (int i = 0, n = u_vector.size(); i<n; i++) {
                if (u_vector[i] != v_vector[i]) {
                    if (candidate_pos != -1) {
                        doSkip = true;
                        break;
                    } else {
                        candidate_pos = i;
                    }
                }
            }

            // Skipping if there is more than just one component that changes
            if (doSkip) continue;

            // Retrieving the graph over which we need to find the edge
            Graph& g = graph_collections[candidate_pos];

            // Getting if the two elements are connected in the original graph: only in this case, I'm going to create the edge in the lattice
            if (size_t cost = g.hasEdge(u_vector[candidate_pos], v_vector[candidate_pos])) {
                file << labels[v] << " " << labels[u] << " " << cost << std::endl;
            }
        }
    }
}
