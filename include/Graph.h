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

#include <NaryTreeGeneration.h>

class Graph {
    size_t rootId;
    lemon::SmartDigraph g;
    std::string filename;
    std::unordered_map<size_t, size_t> invMap;
    std::unordered_map<std::pair<size_t,size_t>, size_t> transitive_closure_map;
    lemon::SmartDigraph::ArcMap<int> costMap;
    lemon::SmartDigraph::NodeMap<bool> rm;
    size_t maxCurrentNode;
    size_t maxLength;
    //lemon::Dfs<lemon::SmartDigraph> dfs;
    //lemon::Dijkstra<lemon::SmartDigraph> dij;
    std::set<size_t> internalLeafCandidates;
    // Getting which is the maximal branching factor of the graph: this is required for the basic metric and for
    // determining the vectorial size
    size_t isTree = std::numeric_limits<size_t>::min();
    fixed_bimap<std::string, size_t> fileElementNameToId;

    inline size_t fw_cost(size_t src, size_t dst) {
        return (src == dst) ? 0 : fw_cost(std::make_pair(src,dst));
    }

    inline size_t fw_cost(const std::pair<size_t,size_t>& cp) {
        if (cp.first == cp.second) return 0;
        auto it = transitive_closure_map.find(cp);
        return it != transitive_closure_map.end() ? it->second : std::numeric_limits<size_t>::max();
    }

    size_t addNewNode(size_t id);
    void addNewEdge(size_t src, size_t dst, int weight);


    void sub_generation_method(std::unordered_map<size_t, size_t>& treeToGraphMorphism,
                               std::unordered_map<size_t, std::unordered_set<size_t>>& tree,
                               std::unordered_map<size_t, std::unordered_set<size_t>>& morphismInv,
                               size_t& minVectors,
                               size_t& vecAverage,
                               std::map<size_t,size_t>& mostFrequent,
                               lemon::SmartDigraph::Node& node);

    /**
     * Under the assumtion that the current graph has a label <a,b> where a is a unique label from the first graph
     * and b is an unique label from the second graph, returns the id associated to the both graph ids
     *
     * @param g1        First graph creating the lattice
     * @param g2        Second graph creating the lattice
     * @param uid       Id of the current graph
     * @return          Pair of node ids, reconstructing the node from the label
     */
    std::pair<unsigned long, unsigned long> getPair(Graph &g1, Graph &g2, size_t uid) const;

public:

    Graph();
    void addEdgeExternally(const std::string& child, const std::string& parent, size_t& count, std::unordered_map<size_t, size_t>& branchingEvaluator, double score = 1.0);

    size_t hasEdge(size_t src, size_t dst) const;
    std::vector<size_t> embedding_id;
    std::unordered_map<size_t, size_t> nodeMap; // maps a node id to a lemon graph id
    // Getting which is the maximal branching factor of the graph: this is required for the basic metric
    // determining the vectorial size
    size_t maxBranch = std::numeric_limits<size_t>::min();
    size_t height;
    void johnsonAlgorithm(bool storeFile = true);

    /**
     * Testing the assumptions that we can
     *
     * @param g1
     * @param g2
     */
    void test_lattice2(Graph &g1, Graph &g2);

    /**
     * Move constructor: this is more like a copy constructor (as lemon does not support moves)
     *
     * @param x
     */
    Graph(Graph&& x);

    /**
     * Loading the graph from a file
     *
     * childid  parentid   edgeWeight
     *
     * @param filename
     */
    Graph(const std::string &filename);


    std::string getName(size_t id) const;
    size_t getId(const std::string& name);


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
                            std::unordered_map<size_t, std::unordered_set<size_t>>& morphismInv, char* rootNameNode = nullptr, bool doReverse = true);

    NaryTreeGeneration generateNaryTree(char *rootNameNode, bool doReverse = true);

    const std::set<size_t>& getCandidates();

    size_t getCost(size_t src, size_t dst, bool isNotLemonId = false);

    template <typename F> size_t isTherePath(size_t dst, std::map<size_t, size_t> &dstCandidates,  F &callback, int& maxLength) {
        // Converting the outer ids to the internal ones
        size_t dstGraph = nodeMap.at(dst);

            std::unordered_set<size_t> noReachedElements;
            maxLength = -1;
            // bool doInsertInUM = dstCandidates.find(dst) == dstCandidates.end();// perform the insertion only if it hasn't been previously inserted for src.
            for (const auto& graph_id_to_original : invMap) {
                size_t graph_id = graph_id_to_original.first;
                size_t original_id = graph_id_to_original.second;
                if  (dstGraph == graph_id) {
                    callback(original_id, 1);          // --
                    dstCandidates[original_id] = 1;          // --
                    maxLength = std::max(maxLength, 1);
                } else {
                    auto cp = std::make_pair(graph_id, dstGraph);
                    auto it = transitive_closure_map.find(cp);
                    if (it != transitive_closure_map.end()) { // reachability: ignoring the self node.
                        //if (doInsertInUM)
                        dstCandidates[original_id] =  it->second+1;  //--
                        callback(original_id, it->second+1);         //--
                        maxLength = std::max(maxLength, (int)it->second+2);
                    } else {
                        // Observe: I'm not providing the element as a candidate.
                        noReachedElements.emplace(original_id);

                    }
                }

            }
            for (const size_t& r : noReachedElements)
                callback(r, maxLength);
            return fw_cost(rootId, dst);

    }


    size_t getRootId();
    int nodeSize() const;

    /**
     * Prints the graph to the standard output
     */
    void print_graph();

    std::vector<size_t> getChildren(size_t i);

    };


#endif //HIERARCHY_TESTS_GRAPH_H
