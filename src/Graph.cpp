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
                               std::unordered_map<size_t, std::unordered_set<size_t>>& morphismInv,
                               char* rootNameNode, bool doReverse) {
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
    ///std::cerr << getName(invMap[rootId]) << " as " << invMap[rootId] << "-->" << rootId << std::endl;
    treeToGraphMorphism[rootId] = invMap[rootId];
    morphismInv[invMap[rootId]].insert(rootId);

    // Visiting the DAG in topological order
    std::vector<lemon::SmartDigraph::Node> nodes;
    lemon::LoggerBoolMap<std::back_insert_iterator<std::vector<lemon::SmartDigraph::Node>>> map(std::back_inserter(nodes));
    lemon::topologicalSort(g, map);
    size_t minVectors = std::numeric_limits<size_t>::max();
    size_t vecAverage = 0;
    std::map<size_t,size_t> mostFrequent;

    if (doReverse) {
        for (auto it = nodes.rbegin(); it != nodes.rend(); it++) {
            lemon::SmartDigraphBase::Node &node = *it;
            sub_generation_method(treeToGraphMorphism, tree, morphismInv, minVectors, vecAverage, mostFrequent, node);
        }
    } else {
        for (auto it = nodes.begin(); it != nodes.end(); it++) {
            lemon::SmartDigraphBase::Node &node = *it;
            sub_generation_method(treeToGraphMorphism, tree, morphismInv, minVectors, vecAverage, mostFrequent, node);
        }
    }

#if 0
    for (auto it = nodes.rbegin(); it != nodes.rend(); it++) {

        lemon::SmartDigraphBase::Node &node = *it;
        sub_generation_method(treeToGraphMorphism, tree, morphismInv, minVectors,vecAverage, mostFrequent, node);


        size_t nId = g.id(node);
        std::string label = getName(invMap.at(nId));
        std::cout << label << std::endl;

        // Checking whether the current element has multiple parents in the dag: if that's the case, then the current
        // node shall be decomposed in as many vectors as many possible parents, thus making as many vectors as many
        // possible paths from the root.
        size_t isMax = 0;
        for (const auto& parent : g.inArcs(node)) {
            const auto& it = morphismInv.find(g.id(g.source(parent)));
            assert(it != morphismInv.end());
            isMax+=it->second.size();
            if (isMax > 1) break;
        }

        if (isMax > 1) {
            std::cout <<
                      fileElementNameToId
                      .getKey(nId) << std::endl;
        }

        for (const auto& parent : g.inArcs(node)) {
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
        for (const auto& descendants : g.outArcs(node)) {
            isLeaf = false;
            break;
        }
        if (isLeaf)

    } internalLeafCandidates.insert(nId);
#endif

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

#ifdef STATS
    std::cout << " branching factor " << maxBranch << " could be 'decreased' by the JL lemma to at least " << dimension_extimate((size_t)g.nodeNum(), 3.0/std::pow(2.0, height)) <<  std::endl;
    std::cout << " tree height = " << height << std::endl;
    std::cout << " maximum number of vectors per node: " << isTree << std::endl;
    std::cout << " minimum number of vectors per node: " << minVectors << std::endl;
    std::cout << " average number of vectors per node: " << ((vecAverage * 1.0) / (g.nodeNum() * 1.0)) << std::endl;
    std::cout << " most frequent vector values: " << mostFrequentValues << " with frequency " << mostFrequentVal << " which, normalized, is " << (((double)mostFrequentVal)/((double)g.nodeNum())) <<  std::endl;
    std::cout << " #nodes " << g.nodeNum() << std::endl;
    std::cout << " #candidate test leaves " << internalLeafCandidates.size() << ", which are " << internalLeafCandidates << std::endl;
#endif

    std::cout << "Johnson (trivial) algorithm for all the possible pairs" << std::endl;
    johnsonAlgorithm();
    return maxBranch;

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
            //std::cout << child << "(" << childId << "~" << nodeMap[childId] <<  ")  --[" << score << "]--> " << parent << "(" << parentId << "~" << nodeMap[parentId] <<  ")" <<  std::endl;
    }
}

void Graph::johnsonAlgorithm(bool storeFile) {
    size_t count = 0;
    if (!transitive_closure_map.empty()) return;
    std::string johnson = filename+"_out.csv";

    if (storeFile && exists(johnson)) {
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
                        //std::cout << getName(invMap[g.id(u)]) << "-[" << d[v] << "]->" << getName(invMap[g.id(v)]) << std::endl;
                        transitive_closure_map[std::make_pair(g.id(u), g.id(v))] = d[v];
                        if (storeFile) file << g.id(u) << " " << g.id(v) << " " << d[v] << std::endl;
                    }
                }
        }

    }
    ///std::cout << "...done!" << std::endl;
}

size_t Graph::getRootId() {
    return rootId;
}

int Graph::nodeSize() const {
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

Graph::Graph(std::vector<Graph> &graph_collections) : Graph{} {

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
                //addNewEdge(u, v,  1);
                /*std::cout << "for graph #" << candidate_pos << ", ";
                std::cout << getName(u) << "-->" <<
                          getName(v) << " ## " <<  u_vector[candidate_pos] << "~>" << v_vector[candidate_pos]<< std::endl;
*/

                //addNewEdge( u, v,  cost);
                //graph_collections[candidate_pos].print_graph();
                //return ;
                /*std::cout << "for graph #" << candidate_pos << " ";
                std::cout << "{";
                for (size_t i = 0, n = u_vector.size(); i<n; i++) {
                    std::cout << graph_collections[i].getName(u_vector[i]);
                    if (i != (n-1)) std::cout << ", ";
                }
                std::cout << "}  {";
                for (size_t i = 0, n = v_vector.size(); i<n; i++) {
                    std::cout << graph_collections[i].getName(v_vector[i]);
                    if (i != (n-1)) std::cout << ", ";
                }
                std::cout << "}" << std::endl;*/
            }
        }
    }

}

void Graph::print_graph() {
    auto it = g.arcs();
    auto itx = it.begin(), itf = it.end();
    while (itx != itf) {
        std::cout <<
                  fileElementNameToId.getKey(
                          invMap.at(lemon::SmartDigraph::id(g.target(itx)))) << " "
                  <<            fileElementNameToId.getKey(
                          invMap.at(lemon::SmartDigraph::id(g.source(itx)))) << " 1"<< std::endl;
        itx++;
    }
}

size_t Graph::hasEdge(size_t src, size_t dst) const {
    auto hasArc =
            lemon::findArc(g, lemon::SmartDigraph::nodeFromId(src), lemon::SmartDigraph::nodeFromId(dst));

    return hasArc != lemon::INVALID ? costMap[hasArc] : 0;
}

void Graph::nested_loop_join(const std::vector<size_t> &embedding_id, std::vector<size_t> vector_pos,
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
            addNewNode(computed_id);
            fileElementNameToId.put(key, computed_id);
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

std::vector<size_t> Graph::generateNode(size_t v) const {
    std::vector<size_t> elements{embedding_id.size(), 0};
    for (size_t i = 0, n = embedding_id.size(); i<n; i++) {
        size_t operand = embedding_id[n-i-1];
        elements[n-i-1] = (size_t)v / operand;
        v = v % operand; /* Likely uses the result of the division. */
    }
    return elements;
}

Graph::Graph(Graph &&x) : Graph{} {
    //std::cout << "copying...." << std::endl;
    //x.print_graph();
    //std::cout << "-----------" << std::endl;
    rootId = x.rootId;
    filename = x.filename;
    ///nodeMap = x.nodeMap;
    ///invMap = x.invMap;
    ///transitive_closure_map = x.transitive_closure_map;
    ///maxCurrentNode = x.maxCurrentNode;
    ///maxLength = x.maxLength;
    ///internalLeafCandidates = x.internalLeafCandidates;
    maxBranch = x.maxBranch;
    height = x.height;
    isTree = x.isTree;
    ///fileElementNameToId = x.fileElementNameToId;
    {
        std::vector<size_t> nodes;
        for (const auto& cp: x.nodeMap) {
            nodes.emplace_back(cp.first);
        }
        std::sort(nodes.begin(), nodes.end());
        for (size_t g_id : nodes) {
            addNewNode(g_id);
            fileElementNameToId.put(x.fileElementNameToId.getKey(g_id), g_id);
        }
    }

    {
        auto ite = x.g.arcs();
        auto iteb = ite.begin(), itee = ite.end();
        while (iteb != itee) {
            addNewEdge(
                    x.invMap[lemon::SmartDigraph::id(x.g.source(iteb))],
                    x.invMap[lemon::SmartDigraph::id(x.g.target(iteb))],
                    x.costMap[iteb]);
            iteb++;
        }
    }

    //lemon::DigraphCopy(x.g, g).run();
   // print_graph();
   // return ;

}

void Graph::test_lattice2(Graph &g1, Graph &g2) {
    for (const auto& cp_cost : transitive_closure_map) {
        size_t uid = cp_cost.first.first;
        size_t vid = cp_cost.first.second;

        std::pair<unsigned long, unsigned long> u = getPair(g1, g2, uid);
        auto v = getPair(g1, g2, vid); //generateNode(invMap.at(vid));

        //std::cout << getName(uid) << " --> " << getName(vid) << " " << cp_cost.second << " == [[" << u.first << "," << u.second << "]]=" << g1.getCost(u.first, v.first, true) << " + [[" << v.first << "," << v.second <<"]]="<< g2.getCost(u.second, v.second, true) << std::endl;
        assert(cp_cost.second == g1.getCost(u.first, v.first, true) + g2.getCost(u.second, v.second, true));
    }
}

std::pair<unsigned long, unsigned long> Graph::getPair(Graph &g1, Graph &g2, size_t uid) const {
    std::string uid_s = getName(uid);
    uid_s = uid_s.substr(1, uid_s.length()-2);
    unsigned long it = uid_s.find_first_of(',');
    std::cout << uid_s << "  " << uid_s.substr(0, it) << "  " << uid_s.substr(it+1) << std::endl;
    auto u = std::make_pair(g1.getId(uid_s.substr(0, it)), g2.getId(uid_s.substr(it+1)));
    return u;
}

size_t Graph::getCost(size_t src, size_t dst, bool isNotLemonId) {
    if (src == dst) return 0;
    auto it = transitive_closure_map.find(std::make_pair(isNotLemonId ? nodeMap.at(src) : src, isNotLemonId ? nodeMap.at(dst) : dst));
    return it == transitive_closure_map.end() ? 0 : it->second;
}

NaryTreeGeneration Graph::generateNaryTree(char *rootNameNode, bool doReverse) {
    NaryTreeGeneration result;
    result.maximum_branching_factor = generateNaryTree(result.treeToGraphMorphism, result.tree, result.treeIdToPathString, result.morphismInv, rootNameNode, doReverse);
    result.number_of_vectors_required = isTree;
    return result;
}

void Graph::sub_generation_method(std::unordered_map<size_t, size_t> &treeToGraphMorphism,
                                  std::unordered_map<size_t, std::unordered_set<size_t>> &tree,
                                  std::unordered_map<size_t, std::unordered_set<size_t>> &morphismInv,
                                  size_t &minVectors, size_t &vecAverage, std::map<size_t, size_t> &mostFrequent,
                                  lemon::SmartDigraph::Node &node) {
    size_t nId = g.id(node);
    std::string label = getName(invMap.at(nId));
    std::cout << label << std::endl;

    // Checking whether the current element has multiple parents in the dag: if that's the case, then the current
    // node shall be decomposed in as many vectors as many possible parents, thus making as many vectors as many
    // possible paths from the root.
    size_t isMax = 0;
    for (const auto& parent : g.inArcs(node)) {
        ///std::cerr << g.id(parent) << "~" << getName(invMap[g.id(parent)]) << " is the parent of " << g.id(node) << "~" << getName(invMap[g.id(node)]) << std::endl;
        const auto& it = morphismInv.find(g.id(g.source(parent)));
        assert(it != morphismInv.end());
        isMax+=it->second.size();
        if (isMax > 1) break;
    }

    if (isMax > 1) {
        std::cout <<
                  fileElementNameToId
                          .getKey(nId) << std::endl;
    }

    for (const auto& parent : g.inArcs(node)) {
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
    for (const auto& descendants : g.outArcs(node)) {
        isLeaf = false;
        break;
    }
    if (isLeaf)
        internalLeafCandidates.insert(nId);
}

