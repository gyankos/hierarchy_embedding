//
// Created by giacomo on 29/12/19.
//

#ifndef HIERARCHY_TESTS_TESTINGGRAPH_H
#define HIERARCHY_TESTS_TESTINGGRAPH_H

#include <vector>
#include <string>
#include <proposal/proposal_utils.h>
#include <string_utils.h>
#include <cout_utils.h>
#include <map>
#include <set>
#include <iostream>
#include <multithreaded/MultithreadWrap.h>
#include <cassert>
#include "PollMap.h"
#include "stats_utils.h"
#include "multithreaded/thread_pool.h"
#include <Graph.h>
#include "TestingTree.h"

#define         DEBUG           (false)

template<typename T, typename iterableT>
std::vector<std::vector<T>> SplitVector(const iterableT& vec, size_t n)
{
    std::vector<std::vector<T>> outVec;

    size_t length = vec.size() / n;
    size_t remain = vec.size() % n;

    size_t begin = 0;
    size_t end = 0;
    auto it = vec.begin();

    for (size_t i = 0; i < std::min(n, vec.size()); ++i) {
        end += (remain > 0) ? (length + !!(remain--)) : length;
        std::vector<T> element;
        while (begin != end) {
            element.emplace_back(*it++);
            begin++;
        }
        outVec.push_back(element);
        begin = end;
    }

    return outVec;
}

template <typename ForComparison> class TestingGraph {
protected:
    size_t maximumBranchingFactor, maximumHeight;
    Graph& passGraph;

    /// Correspondence between the node in adj and the node in the graph, so that we can reconstruct all the vectors to be associated to it
    std::unordered_map<size_t, size_t> treeToGraphMorphism;
    std::unordered_map<size_t, std::unordered_set<size_t>> morphismInv;

    /// Tree representation of the data structure, with new element ids
    std::unordered_map<size_t, std::unordered_set<size_t>> tree;

    /// path string generated from the adj representation. This can be used to then generate multiple possible elements
    std::map<size_t, std::string> treeIdToPathString;

private:
    double Spearman(const size_t& current, const std::map<size_t, size_t>& currentRankMap) {
        double noCandidates = currentRankMap.size();
        double retrieved = 0.0;
        std::map<size_t, size_t> toretMap;

        PollMap<double, size_t> pollMap{currentRankMap.size()};
        generateTopKCandidates(pollMap, current);
        pollMap.getRankedPoll(toretMap);

        return Spearmam_rank_correlation_coefficient<size_t>(toretMap, currentRankMap);
    }

    double NDCG(const size_t& current, std::map<size_t, size_t>& currentRankMap) {
        double noCandidates = currentRankMap.size();
        double retrieved = 0.0;
        std::map<size_t, size_t> toretMap/*, currentRankMap*/;

        PollMap<double,size_t> pollMap{currentRankMap.size()};
        generateTopKCandidates(pollMap, current);
        const std::map<double, std::set<size_t>> &poll = pollMap.getPoll();

        return normalized_discounted_cumulative_gain<size_t>(poll, currentRankMap);
    }

    double Precision(const size_t& current, const std::map<size_t, size_t>& currentRankMap) {
        double retrieved = 0.0;
        std::map<std::string, size_t> toretMap/*, currentRankMap*/;

        PollMap<double, size_t> pollMap{currentRankMap.size()};
        generateTopKCandidates(pollMap, current);
        std::multimap<double, size_t> mmap;
        pollMap.getPoll(mmap);
        double noCandidates = mmap.size();
        std::set<size_t> k;
        for (auto it = mmap.rbegin(); it != mmap.rend(); it++) {
            if (currentRankMap.find(it->second) != currentRankMap.end()) {
                if(k.insert(it->second).second)
                    retrieved+=1.0;
            }
        }
        assert(retrieved <= noCandidates);
        return retrieved / noCandidates;
    }

    double Precision_narrow(const size_t& current, const std::map<size_t, size_t>& currentRankMap) {
        //std::set<std::string> candidatesSet;
        double retrieved = 0.0;
        std::map<size_t, size_t> toretMap;

        double currentlyGot = 0.0;
        PollMap<double, size_t> pollMap{currentRankMap.size()};
        generateTopKCandidates(pollMap, current);
        pollMap.getNarrowRankedPoll(toretMap);
        std::set<size_t> k;

        for (auto& cp : toretMap) {
            if (currentRankMap.find(cp.first) != currentRankMap.end()) {
                if (k.insert(cp.first).second)
                    retrieved += 1.0;
            }
            currentlyGot += 1.0;
        }
        assert(retrieved <= currentlyGot);
        return retrieved / currentlyGot;
    }


    std::pair<double,double> minmaxCandidates(const size_t& current/*, const std::vector<std::vector<std::vector<size_t>>>& z,*/, const std::map<size_t, size_t>& currentRankMap) {
        double min = std::numeric_limits<double>::max();
        std::vector<size_t> minCandidates;
        double max = std::numeric_limits<double>::min();
        std::vector<size_t> maxCandidates;

        //std::string currentStirng = size_vector_to_string(current);
        //for (auto& candidates : z) {
            for (const std::pair<const unsigned long, unsigned long>& x : currentRankMap) {
                //std::string x_tmp = size_vector_to_string(x);
                if (x.first != current) {
                    double z = similarity(getVectorRepresentation(current), getVectorRepresentation(x.second));/// ???????
                    if (z < min) {
                        minCandidates.clear();
                        minCandidates.emplace_back(x.first);
                        min = z;
                    } else if (z == min) {
                        minCandidates.emplace_back(x.first);
                    }
                    if (z > max) {
                        maxCandidates.clear();
                        maxCandidates.emplace_back(x.first);
                    } else if (z == max) {
                        maxCandidates.emplace_back(x.first);
                    }
                    max = std::max(max, z);
                }
            }
        //}

        double recallToBeZero = 0.0;


        for (const size_t& current : minCandidates) {
            //std::string current = size_vector_to_string(x);
            if (currentRankMap.find(current) != currentRankMap.end())
                recallToBeZero += 1.0;
        }

        return std::make_pair(recallToBeZero/minCandidates.size(), maxCandidates.size()/(1.0+maxCandidates.size()));
    }


public:
    TestingGraph(Graph& ref) : passGraph{ref} {}


    void run(bool isMultithreaded) {
        treeToGraphMorphism.clear();
        morphismInv.clear();
        tree.clear();
        treeIdToPathString.clear();

        // Filling up all the three data structures provided up above.
        passGraph.generateNaryTree(treeToGraphMorphism, tree, treeIdToPathString, morphismInv);
        maximumBranchingFactor = passGraph.maxBranch;
        maximumHeight = passGraph.height;
        passGraphDataIfRequired(passGraph);

        double path_length_size = 0, spearman = 0,  precision = 0,  precision_narrow = 0, ncdg = 0, smallerNotCandidate = 0, recall =0;

        size_t threadNo = (unsigned int)std::thread::hardware_concurrency();
        MultithreadWrap<struct result_map> pool{(unsigned int)std::thread::hardware_concurrency(), isMultithreaded};
        std::vector<std::vector<size_t>> splitNodeCandidatesId = SplitVector<size_t, std::set<size_t>>(passGraph.getCandidates(), threadNo);


        //std::vector<std::future<struct result_map >> futures;
        //thread_pool pool(ls.size());

        std::cout << "Now Computing..." << std::endl;
        // Performing the test for each node in the hierarchy, root excluded
        for (const std::vector<size_t>& y: splitNodeCandidatesId) {
            pool.poolExecute([this, splitNodeCandidatesId](const std::vector<size_t> & candidateSublist) {
                struct result_map maps;
                for (const size_t & candidateId : candidateSublist) {
                    //size_t current_path_length = candidateId.size();
                    //if (current_path_length != maximumHeight) continue;

                    // Increment the path length size by one and, if not present, insert 1.0
                    std::map<size_t, size_t> rankedCandidates;
                    this->passGraph.isTherePath(candidateId, rankedCandidates);
                    maps.path_length_size += 1.0;

                    // Generating all the candidates, represented as indices, as relevant is-a nodes for the current candidate.
                    //const std::vector<std::vector<size_t>> &candidates = generateAllPossibleSubpaths(x);
                    std::pair<double,double> minmax = this->minmaxCandidates(candidateId, /*splitNodeCandidatesId, */rankedCandidates);

                    // Testing the precision at k representation
                    //double currentPrecision = this->Spearman(x, candidates);
                    //if (DEBUG)  std::cout << "Current precision of " << currentPrecision << " for length " << current_path_length << std::endl;
                    maps.spearman += this->Spearman(candidateId, rankedCandidates); // Getting the classification difference between the two
                    maps.precision += this->Precision(candidateId, rankedCandidates);
                    maps.precision_narrow += this->Precision_narrow(candidateId, rankedCandidates);
                    maps.ncdg += this->NDCG(candidateId, rankedCandidates);
                    maps.recall += minmax.first;
                    maps.smallerNotCandidate += minmax.second;
                }

                return maps;
            }, y);
        }

        std::cout << "Summing up the maps..." << std::endl;
        for (auto& x : pool.foreach()) {
            //auto x = y.get();
#define update_local(field, currFuture)           (field) += currFuture .  field
            update_local(path_length_size, x);
            update_local(spearman, x);
            update_local(precision, x);
            update_local(precision_narrow, x);
            update_local(ncdg, x);
            update_local(recall, x);
        }

        std::cout << "Spearman, for top-k                                              " << spearman / path_length_size << std::endl;
        //print_result_maps(path_length_size, spearman);
        std::cout << "Precision, for top-k scores                                      " << precision / path_length_size << std::endl;
        //print_result_maps(path_length_size, precision);
        std::cout << "Precision, for top-k elements                                    " << precision_narrow / path_length_size << std::endl;
        std::cout << "Recall, for non top-k elements (wrongly matching the candidates) " << recall / path_length_size << std::endl;
        //std::cout << "More similar than worst top-1 not candidate, not current element " << smallerNotCandidate / path_length_size << std::endl;
        //print_result_maps(path_length_size, precision_narrow);
    }

protected:
    /// ???? virtual void initialize_hierarchy_with_all_paths(const std::vector<std::vector<size_t>> &subgraph_as_paths) = 0;
    virtual ForComparison getVectorRepresentation(const size_t& current) = 0;
    virtual double similarity(const ForComparison& lhs, const ForComparison& rhs) = 0;
    virtual void  generateTopKCandidates(PollMap<double, size_t>& map, const size_t& current) = 0;

    /**
     * ?????
     * @param graph
     */
    virtual void passGraphDataIfRequired(const Graph& graph) = 0;
};


#endif //HIERARCHY_TESTS_TESTINGTREE_H
