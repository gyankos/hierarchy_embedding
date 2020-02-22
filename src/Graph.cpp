//
// Created by giacomo on 31/12/19.
//

#include <sys/stat.h>
#include "Graph.h"
#include "cout_utils.h"

inline bool exists(const std::string& name) {
    struct stat buffer;
    return (stat (name.c_str(), &buffer) == 0);
}

size_t generateAllPathIds(size_t id, std::unordered_map<size_t, std::unordered_set<size_t>> &adj,
                          const std::vector<size_t> &currentVector, std::map<size_t, std::string> &result, size_t height) {
    size_t currentIterationStep = 1;
    result[id] = size_vector_to_string(currentVector);
    size_t maxHeight = height;
    for (const size_t& out : adj[id]) {
        std::vector<size_t> copy{currentVector};
        copy.emplace_back(currentIterationStep++);
        maxHeight = std::max(maxHeight, generateAllPathIds(out, adj, copy, result, height+1));
    }
    return maxHeight;
}

size_t Graph::addNewNode(size_t id) {
    auto it = nodeMap.find(id);
    if (it != nodeMap.end())
        return it->second;
    lemon::SmartDigraph::Node currentNode = g.addNode();
    maxCurrentNode = g.id(currentNode);
    invMap[maxCurrentNode] = id;
    return nodeMap.insert(std::make_pair(id, maxCurrentNode)).first->second;
}

void Graph::addNewEdge(size_t src, size_t dst, int weight) {
    size_t src2 = addNewNode(src), dst2 = addNewNode(dst);
    lemon::SmartDigraph::Node srcNode = g.nodeFromId(src2), dstNode = g.nodeFromId(dst2);
    costMap[g.addArc(srcNode, dstNode)] = weight;
}

size_t Graph::generateNaryTree(std::unordered_map<size_t, size_t> &treeToGraphMorphism,
                               std::unordered_map<size_t, std::unordered_set<size_t>> &tree,
                               std::map<size_t, std::string> &treeIdToPathString,

                               std::unordered_map<size_t, std::unordered_set<size_t>>& morphismInv, char* rootNameNode) {
    //std::unordered_map<size_t, std::unordered_set<size_t>> morphismInv;

    // Ensuring that the loaded dataset is indeed a DAG
    assert(lemon::dag(g));

    // -- Alternative: make the data structure a dag itselfc
    //lemon::SmartDigraph::NodeMap<bool> nm(g, true);
    //lemon::SmartDigraph::ArcMap<bool> em(g, true);
    //dag(g, em); // Pruning the graph
    //auto subGraph = lemon::subDigraph(g, nm, em);
    // --

    // Assuming that the DAG is a hierarchy, and therefore there should be just one root element.
    std::vector<size_t> noParentNodes;
    for (const auto& n : g.nodes()) {
        bool hasParent = false;
        for (const auto& parent : g.inArcs(n)) {
            hasParent = true; break;
        }
        if (!hasParent) noParentNodes.emplace_back(g.id(n));
    }
    // this node shall correspond to "entity"
    //std::cout << invMap[noParentNodes[0]] << std::endl;

    if (!rootNameNode)
        rootId = noParentNodes[0];
    else {
        std::string rootName = rootNameNode;
        rootId = getId(rootNameNode);
       /// std::cout << "Setting " << rootName << " with id " << rootId << "as the root" << std::endl;
    }
    treeToGraphMorphism[rootId] = invMap[rootId];
    morphismInv[invMap[rootId]].insert(rootId);

    // Visiting the DAG in topological order
    std::vector<lemon::SmartDigraph::Node> nodes;
    lemon::LoggerBoolMap<std::back_insert_iterator<std::vector<lemon::SmartDigraph::Node>>> map(std::back_inserter(nodes));
    lemon::topologicalSort(g, map);
    size_t minVectors = std::numeric_limits<size_t>::max();
    size_t vecAverage = 0;
    std::map<size_t,size_t> mostFrequent;

    for (auto it = nodes.rbegin(); it != nodes.rend(); it++) {
        size_t nId = g.id(*it);
        if (nId == 21)
            std::cerr << "DEBUG" << std::endl;

        // Checking whether the current element has multiple parents in the dag: if that's the case, then the current
        // node shall be decomposed in as many vectors as many possible parents, thus making as many vectors as many
        // possible paths from the root.
        size_t isMax = 0;
        for (const auto& parent : g.inArcs(*it)) {
            const auto& it = morphismInv.find(g.id(g.source(parent)));
            assert(it != morphismInv.end());
            isMax+=it->second.size();
            if (isMax > 1) break;
        }

        for (const auto& parent : g.inArcs(*it)) {
            const auto& ita = morphismInv.find(g.id(g.source(parent)));
            assert(ita != morphismInv.end());
            for (const size_t& parentMaps : ita->second) {
                size_t currentVertex = (isMax <= 1) ? nId : ++maxCurrentNode;
                treeToGraphMorphism[currentVertex] = nId;
                morphismInv[nId].insert(currentVertex);
                tree[parentMaps].insert(currentVertex);
                maxBranch = std::max(maxBranch, tree[parentMaps].size());
            }
            isTree = std::max(isTree, morphismInv[nId].size());
            minVectors = std::min(minVectors, morphismInv[nId].size());
        }
        vecAverage += morphismInv[nId].size();
        auto it2 = mostFrequent.insert(std::make_pair(morphismInv[nId].size(), 1));
        if (!it2.second)
            it2.first->second++;

        // Determining which are the leaves: those will be used for the testing
        bool isLeaf = true;
        for (const auto& descendants : g.outArcs(*it)) {
            isLeaf = false;
            break;
        }
        if (isLeaf)
            internalLeafCandidates.insert(nId);
    }

    height = generateAllPathIds(rootId, tree, {}, treeIdToPathString, 1);

    std::set<size_t> mostFrequentValues;
    int mostFrequentVal = -1;
    for (auto& it3 : mostFrequent) {
        if (((int)it3.second) > mostFrequentVal) {
            mostFrequentValues.clear();
            mostFrequentVal = it3.second;
            mostFrequentValues.emplace(it3.first);
        } else if (((int)it3.second) == mostFrequentVal) {
            mostFrequentValues.emplace(it3.first);
        }
    }
/**
    std::cout << " branching factor " << maxBranch << " could be 'decreased' by the JL lemma to at least " << dimension_extimate((size_t)g.nodeNum(), 3.0/std::pow(2.0, height)) <<  std::endl;
    std::cout << " tree height = " << height << std::endl;
    std::cout << " maximum number of vectors per node: " << isTree << std::endl;
    std::cout << " minimum number of vectors per node: " << minVectors << std::endl;
    std::cout << " average number of vectors per node: " << ((vecAverage * 1.0) / (g.nodeNum() * 1.0)) << std::endl;
    std::cout << " most frequent vector values: " << mostFrequentValues << " with frequency " << mostFrequentVal << " which, normalized, is " << (((double)mostFrequentVal)/((double)g.nodeNum())) <<  std::endl;
    std::cout << " #nodes " << g.nodeNum() << std::endl;
    std::cout << " #candidate test leaves " << internalLeafCandidates.size() << ", which are " << internalLeafCandidates << std::endl;

    std::cout << "Johnson (trivial) algorithm for all the possible pairs" << std::endl;*/
    johnsonAlgorithm();
///    std::cout << "... done! " << std::endl;

}



const std::set<size_t> &Graph::getCandidates() {
    return internalLeafCandidates;
}

Graph::Graph() : costMap{g}, /*dfs{g}, dij{g, costMap},*/ rm{g} {
   // dfs.reachedMap(rm);
}

Graph::Graph(const std::string &filename) : Graph{} {
    this->filename = filename;
    std::string child, parent;
    std::ifstream file{filename};
    double score;

    std::unordered_map<size_t, size_t> branchingEvaluator;
    size_t count = 0;

    while (file >> child >> parent >> score) {
        size_t childId, parentId; //ensuring that the root will be the one with id 0
        if ((parentId = fileElementNameToId.putWithKeyCheck(parent, count)) == count) {
            addNewNode(count++);
        }
        if ((childId = fileElementNameToId.putWithKeyCheck(child, count)) == count) {
            addNewNode(count++);
        }
        auto it = branchingEvaluator.find(parentId);
        if (it == branchingEvaluator.end()) {
            maxBranch = std::max(maxBranch, (size_t)1);
            branchingEvaluator.insert(std::make_pair(parentId, 1));
        } else {
            maxBranch = std::max(++it->second, maxBranch);
        }
        addNewEdge(parentId, childId, score);
        /*if (childId == 22708 || parentId == 22708)*/
            std::cout << child << "(" << childId <<  ")  --[" << score << "]--> " << parent << "(" << parentId << ")" <<  std::endl;
    }
}

void Graph::johnsonAlgorithm() {
    size_t count = 0;
    std::string johnson = filename+"_out.csv";

    if (exists(johnson)) {
        ///std::cout << "Loading the precomputed outcome of the johnson Algorithm" << std::endl;
        std::ifstream file{johnson};
        size_t src,dst,w;
        while (file >> src >> dst >> w) {
            transitive_closure_map[std::make_pair(src, dst)] = w;
        }
    } else {
        ///std::cout << "Running Dijkstra for each node in the graph" << std::endl;//In fact, we don't need to run Bellman-Ford's algorithm, as we have no negative edges in here.
        // Still, this algorithm will be more efficient than running Floyd-Washall's algorithm, for the DAG is sparse
        std::ofstream file{johnson};
        for (auto& u : g.nodes()) { // For each vertex u in the graph
                if (!(count % 1000)) {
                    ///std::cout << count+1 << std::endl;
                }
                count++;
                lemon::Dfs<lemon::SmartDigraph> dfs{g};
                lemon::Dijkstra<lemon::SmartDigraph> dij{g, costMap};
                dfs.init();
                dfs.addSource(u);
                dfs.start();
                dij.init();
                dij.addSource(u);
                dij.start();

                auto& d = dij.distMap(); // Get the distance map for the current run of Dijkstra

                for (auto& v : g.nodes()) { // Setting up the distance map
                    if (dij.reached(v) && (v != u)) {
                        transitive_closure_map[std::make_pair(g.id(u), g.id(v))] = d[v];
                        file << g.id(u) << " " << g.id(v) << " " << d[v] << std::endl;
                    }
                }
        }

    }
    ///std::cout << "...done!" << std::endl;
}

size_t Graph::getRootId() {
    return rootId;
}

int Graph::nodeSize() {
    return nodeMap.size();
}

std::vector<size_t> Graph::getChildren(size_t i) {
    std::vector<size_t> toRet;
    auto it = g.outArcs(g.nodeFromId(i));
    auto itx = it.begin(), itf = it.end();
    while (itx != itf) {
        toRet.emplace_back(g.id(itx));
        itx++;
    }
    return toRet;
}

