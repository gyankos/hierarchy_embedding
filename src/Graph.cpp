//
// Created by giacomo on 31/12/19.
//

#include "Graph.h"
#include "cout_utils.h"

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
    std::unordered_map<size_t, size_t>::iterator it = nodeMap.find(id);
    if (it != nodeMap.end())
        return it->second;
    lemon::SmartDigraph::Node currentNode = g.addNode();
    maxCurrentNode = g.id(currentNode);
    invMap[maxCurrentNode] = id;
    return nodeMap.insert(std::make_pair(id, maxCurrentNode)).first->second;
}

void Graph::addNewEdge(size_t src, size_t dst, int weight) {
    lemon::SmartDigraph::Node srcNode = g.nodeFromId(addNewNode(src)), dstNode = g.nodeFromId(addNewNode(dst));
    costMap[g.addArc(srcNode, dstNode)] = weight;
}

size_t Graph::generateNaryTree(std::unordered_map<size_t, size_t> &treeToGraphMorphism,
                               std::unordered_map<size_t, std::unordered_set<size_t>> &tree,
                               std::map<size_t, std::string> &treeIdToPathString,

                               std::unordered_map<size_t, std::unordered_set<size_t>>& morphismInv) {
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

    rootId = noParentNodes[0];
    treeToGraphMorphism[rootId] = invMap[rootId];
    morphismInv[invMap[rootId]].insert(rootId);

    // Visiting the DAG in topological order
    std::vector<lemon::SmartDigraph::Node> nodes;
    lemon::LoggerBoolMap<std::back_insert_iterator<std::vector<lemon::SmartDigraph::Node>>> map(std::back_inserter(nodes));
    lemon::topologicalSort(g, map);


    for (auto it = nodes.rbegin(); it != nodes.rend(); it++) {
        size_t nId = g.id(*it);

        // Checking whether the current element has multiple parents in the dag: if that's the case, then the current
        // node shall be decomposed in as many vectors as many possible parents, thus making as many vectors as many
        // possible paths from the root.
        size_t isMax = 0;
        for (const auto& parent : g.inArcs(*it)) {
            isMax++;
            if (isMax > 1) break;
        }
        //std::cout << nId << " with  isMax = " << isMax <<  std::endl;

        for (const auto& parent : g.inArcs(*it)) {
            size_t currentVertex = (isMax <= 1) ? nId : ++maxCurrentNode;
            treeToGraphMorphism[currentVertex] = nId;
            morphismInv[nId].insert(currentVertex);
            isTree = std::max(isTree, morphismInv[nId].size());
            const auto& it = morphismInv.find(g.id(g.source(parent)));
            assert(it != morphismInv.end());
            for (const size_t& parentMaps : it->second) {
                tree[parentMaps].insert(currentVertex);
                maxBranch = std::max(maxBranch, tree[parentMaps].size());
            }
        }

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

    std::cout << " branching factor " << maxBranch << " could be 'decreased' by the JL lemma to at least " << dimension_extimate((size_t)g.nodeNum(), 3.0/std::pow(2.0, height)) <<  std::endl;
    std::cout << " tree height = " << height << std::endl;
    std::cout << " maximum number of vectors per node: " << isTree << std::endl;
    std::cout << " #nodes " << g.nodeNum() << std::endl;
    std::cout << " #candidate test leaves " << internalLeafCandidates.size() << ", which are " << internalLeafCandidates << std::endl;
}

size_t Graph::isTherePath(size_t dst,
                          std::map<size_t, size_t> &dstCandidates) {
    // Converting the outer ids to the internal ones
    dst = nodeMap[dst];

    // If I haven't run the search from this node, then initialize the search
    if (rootId != sourceDfsNode) {
        sourceDfsNode = rootId;
        dfs.init();
        dfs.addSource(g.nodeFromId(sourceDfsNode));
        dfs.start();
        dij.init();
        dij.addSource(g.nodeFromId(sourceDfsNode));
        dij.start();
        maxLength = std::numeric_limits<size_t>::min();
        std::vector<size_t> noReachedElements;
       // bool doInsertInUM = dstCandidates.find(dst) == dstCandidates.end();// perform the insertion only if it hasn't been previously inserted for src.
        for (const auto& node : g.nodes()) {
            if (dij.reached(node)) { // reachability
                //if (doInsertInUM)
                dstCandidates[invMap[g.id(node)]] =  dfs.dist(node);
                maxLength = std::max(maxLength, (size_t)dfs.dist(node)+1);
            } else {
                noReachedElements.emplace_back(g.id(node));
            }
            for (const auto& notInPath : noReachedElements)
                dstCandidates[notInPath] = maxLength;
        }
    }
    return dij.reached(g.nodeFromId(dst)) ? dij.dist(g.nodeFromId(dst)) : maxLength;
}

const std::set<size_t> &Graph::getCandidates() {
    return internalLeafCandidates;
}

Graph::Graph() : costMap{g}, dfs{g}, dij{g, costMap}, rm{g} {
    dfs.reachedMap(rm);
}

Graph::Graph(std::ifstream &file) : Graph{} {
    std::string child, parent;
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
        /*if (childId == 22708 || parentId == 22708)
            std::cout << child << "(" << childId <<  ")  --[" << score << "]--> " << parent << "(" << parentId << ")" <<  std::endl;*/
    }
}
