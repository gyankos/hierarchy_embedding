//
// Created by giacomo on 31/12/19.
//

#ifndef HIERARCHY_TESTS_GRAPH_H
#define HIERARCHY_TESTS_GRAPH_H

#include <lemon/smart_graph.h>
#include <lemon/connectivity.h>
#include <lemon/dijkstra.h>
#include <lemon/dfs.h>

#include <unordered_set>
#include <unordered_map>
#include <cassert>
#include <string_utils.h>
#include <JLLemma.h>
#include <proposal/Proposal.h>
#include "naryTree.h"
#include "fixed_bimap.h"
#include "cout_utils.h"

/**
 * Pruning the given hierarchy as a DAG using the DFS visit. This is useful when the hierarchies that you receive as an
 * input are not perfect DAGs
 *
 * @tparam Digraph
 * @tparam ArcMap
 *
 * @param digraph
 * @param pruned     edge map for pruning the graph
 * @return
 */
template <typename Digraph, typename ArcMap> bool dag(const Digraph& digraph, ArcMap& pruned) {
    lemon::checkConcept<lemon::concepts::Digraph, Digraph>();
    typedef typename Digraph::Node Node;
    typedef typename Digraph::NodeIt NodeIt;
    typedef typename Digraph::Arc Arc;

    typedef typename Digraph::template NodeMap<bool> ProcessedMap;

    typename lemon::Dfs<Digraph>::template SetProcessedMap<ProcessedMap>::
    Create dfs(digraph);

    ProcessedMap processed(digraph);
    dfs.processedMap(processed);

    dfs.init();
    for (NodeIt it(digraph); it != lemon::INVALID; ++it) {
        if (!dfs.reached(it)) {
            dfs.addSource(it);
            while (!dfs.emptyQueue()) {
                Arc arc = dfs.nextArc();
                Node target = digraph.target(arc);
                if (dfs.reached(target) && !processed[target]) {
                    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << digraph.id(target) << "~~~~~~~~~~~~~~~~" << std::endl;
                    pruned.set(arc, false);
                }
                dfs.processNextArc();
            }
        }
    }
    return true;
}

size_t generateAllPathIds(size_t id, std::unordered_map<size_t, std::unordered_set<size_t>>& adj, const std::vector<size_t>& currentVector, std::map<size_t, std::string> &result, size_t height = 1);

#include <fstream>
#include <iostream>

namespace std {
    template <> struct hash<std::pair<size_t,size_t>> {
            size_t operator()(const pair<size_t, size_t>& p) const
            {
                auto hash1 = hash<size_t>{}(p.first);
                auto hash2 = hash<size_t>{}(p.second);
                return hash1 ^ hash2;
            }
    };
}

class Graph {
    size_t rootId;
    lemon::SmartDigraph g;
    std::string filename;
    std::unordered_map<size_t, size_t> nodeMap;
    std::unordered_map<size_t, size_t> invMap;
    std::unordered_map<std::pair<size_t,size_t>, size_t> transitive_closure_map;
    lemon::SmartDigraph::ArcMap<int> costMap;
    lemon::SmartDigraph::NodeMap<bool> rm;
    size_t maxCurrentNode;
    size_t sourceDfsNode = std::numeric_limits<size_t>::min();
    size_t maxLength;
    //lemon::Dfs<lemon::SmartDigraph> dfs;
    //lemon::Dijkstra<lemon::SmartDigraph> dij;
    std::set<size_t> internalLeafCandidates;

    inline size_t fw_cost(size_t src, size_t dst) {
        return (src == dst) ? 0 : fw_cost(std::make_pair(src,dst));
    }

    inline size_t fw_cost(const std::pair<size_t,size_t>& cp) {
        if (cp.first == cp.second) return 0;
        auto it = transitive_closure_map.find(cp);
        return it != transitive_closure_map.end() ? it->second : std::numeric_limits<size_t>::max();
    }

    void johnsonAlgorithm();

    // Getting which is the maximal branching factor of the graph: this is required for the basic metric and for
    // determining the vectorial size
    size_t isTree = std::numeric_limits<size_t>::min();

    fixed_bimap<std::string, size_t> fileElementNameToId;

    size_t addNewNode(size_t id);
    void addNewEdge(size_t src, size_t dst, int weight);
    Graph();

public:
    Graph(const std::string &filename);


    /**
     *
     * Also, it initializes the internalLeafCandidates from which
     *
     * @param treeToGraphMorphism          Correspondence between the node in adj and the node in the graph, so that we can reconstruct all the vectors to be associated to it
     * @param tree                         Tree representation of the data structure, with new element ids
     * @param treeIdToPathString           path string generated from the adj representation. This can be used to then generate multiple possible elements
     * @param morphismInv                  Inverse morphism representation
     * @return                  Maximum branching factor of the associated tree
     */
    size_t generateNaryTree(std::unordered_map<size_t, size_t>& treeToGraphMorphism, std::unordered_map<size_t, std::unordered_set<size_t>>& tree, std::map<size_t, std::string>& treeIdToPathString,
                            std::unordered_map<size_t, std::unordered_set<size_t>>& morphismInv);
    const std::set<size_t>& getCandidates();
    size_t isTherePath(size_t dst, std::map<size_t, size_t> &dstCandidates);

   // Getting which is the maximal branching factor of the graph: this is required for the basic metric
   // determining the vectorial size
    size_t maxBranch = std::numeric_limits<size_t>::min();
    size_t height;
};


#endif //HIERARCHY_TESTS_GRAPH_H
