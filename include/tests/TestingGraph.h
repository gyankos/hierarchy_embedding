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

template <typename ForComparison>
struct TestingGraphLambda {
    const size_t& current;
    std::function<ForComparison(size_t)>& getVectorRepresentation;
    std::function<double(const ForComparison&,const ForComparison&)>& similarity;
    std::map<size_t, size_t>& currentRankMap;

    PollMap<double, size_t> pollMap;

    TestingGraphLambda(const size_t& current, std::map<size_t, size_t>& currentRankMap, std::function<ForComparison(size_t)>& getVectorRepresentation,
                       std::function<double(const ForComparison&,const ForComparison&)>& similarity) :
    current{current}, currentRankMap{currentRankMap}, getVectorRepresentation{getVectorRepresentation}, similarity{similarity} {

    }

    size_t map_size = 0;

    double min_minmax = std::numeric_limits<double>::max();
    std::vector<size_t> minCandidates_minmax;
    double max_minmax = std::numeric_limits<double>::min();
    std::vector<size_t> maxCandidates_minmax;

    std::set<size_t>    precision;
    std::set<size_t>    candidates;
    size_t              retrieved_for_precision = 0;
    size_t              noCandidates;


    double              expectedMeanPollRank = 0.0;

    void operator()(size_t first, size_t second) {
        //std::cout << first << "," << second << std::endl;
        double score = similarity(getVectorRepresentation(current), getVectorRepresentation(first));

        // Spearman
        expectedMeanPollRank += second;

        // MinMax computation
        //if (first != current) {
            /// ???????
            if (score < min_minmax) {
                minCandidates_minmax.clear();
                minCandidates_minmax.emplace_back(first);
                min_minmax = score;
            } else if (score == min_minmax) {
                minCandidates_minmax.emplace_back(first);
            }
            if (score > max_minmax) {
                maxCandidates_minmax.clear();
                maxCandidates_minmax.emplace_back(first);
            } else if (score == max_minmax) {
                maxCandidates_minmax.emplace_back(first);
            }
            max_minmax = std::max(max_minmax, score);
        //}
        // MiniMax computation

        // Globals
        pollMap.add(score, first);
        map_size++;
    }

    void finalize(struct result_map& maps) {
        pollMap.resize(currentRankMap.size());

        // MiniMax computation
        double recallToBeZero = 0.0;
        for (const size_t& current : minCandidates_minmax) {
            //std::string current = size_vector_to_string(x);
            if (currentRankMap.find(current) != currentRankMap.end())
                recallToBeZero += 1.0;
        }
        maps.recall += recallToBeZero / minCandidates_minmax.size();
        //maps.smallerNotCandidate +=  maxCandidates_minmax.size() / (1.0 + maxCandidates_minmax.size());
        // MiniMax computation


        // Globals: generating the candidates that are just the ones returned by the algorithm itself
        maps.path_length_size += 1.0;
        //PollMap<double,size_t> pollMap{map_size};
        //generateTopKCandidates(pollMap, current);
        const std::map<double, std::set<size_t>> &castor_et_pollux = pollMap.getPoll();

        // Spearman + NCDG
        double NDCG_pos = 1.0;
        double NDCG_dcg_p = 0.0;
        double NDCG_idcg_p = 0.0;

        expectedMeanPollRank /= ((double)map_size);

        size_t inferredPollCurrentPos = castor_et_pollux.size();
        double inferredMeanPollRank = 0.0;
        size_t overallX = 0;
        std::vector<double> precalculateX;
        std::vector<double>    calculateY;
        double                 calculateYSquaredSummed = 0.0;
        for (const auto& it : castor_et_pollux) {
            size_t N = it.second.size();
            bool doNCDG = false/*(N > 0)*/;
            double div;
            double toSum = 0.0;
            double reli;

            overallX += N;
            for (const auto& x : it.second) {
                if (currentRankMap.find(x) != currentRankMap.end()) {
                    if(candidates.insert(x).second)
                        retrieved_for_precision+=1.0;
                }

                double tmpY = ((double)currentRankMap[x])-expectedMeanPollRank;
                calculateYSquaredSummed += std::pow(tmpY, 2.0);
                calculateY.emplace_back(tmpY);
                precalculateX.emplace_back(inferredPollCurrentPos);
                toSum += (double)currentRankMap[x];
            }
            for (size_t i = 0; i<N; i++)
            inferredMeanPollRank += (inferredPollCurrentPos * N);
            inferredPollCurrentPos--;

        }
        inferredMeanPollRank /= ((double)overallX);

        /// Precision
        maps.precision += ((double)retrieved_for_precision)/((double)map_size);

        double numeratore = 0.0;
        double                 calculateXSquaredSummed = 0.0;
        for (size_t i = 0; i<overallX; i++) {
            double tmpX = (precalculateX[i] - inferredMeanPollRank);
            calculateXSquaredSummed += std::pow(tmpX, 2.0);
            numeratore += tmpX * calculateY[i];
        }

        maps.spearman += (overallX > 1) ? numeratore / std::sqrt(calculateXSquaredSummed * calculateYSquaredSummed) : 1.0-abs(inferredMeanPollRank-expectedMeanPollRank)/(1+abs(inferredMeanPollRank-expectedMeanPollRank));
    }

 };

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

#if 0
    PollMap<double, size_t> getMap(const size_t &current, const std::map<size_t, size_t> &currentRankMap) {
        PollMap<double,size_t> pollMap{currentRankMap.size()};
        generateTopKCandidates(pollMap, current);
        return pollMap;
    }

    double Spearman(const size_t& current, const std::map<size_t, size_t>& currentRankMap) {
        double noCandidates = currentRankMap.size();
        double retrieved = 0.0;
        std::map<size_t, size_t> toretMap;


        PollMap<double, size_t> pollMap = getMap(current, currentRankMap);
        pollMap.getRankedPoll(toretMap);

        return Spearmam_rank_correlation_coefficient<size_t>(toretMap, currentRankMap);
    }

    // GetPoll: requiresIterationOver
    double NDCG(const size_t& current, std::map<size_t, size_t>& currentRankMap) {
        double noCandidates = currentRankMap.size();
        double retrieved = 0.0;

        /// COMMON
        PollMap<double, size_t> pollMap = getMap(current, currentRankMap);
        const std::map<double, std::set<size_t>> &poll = pollMap.getPoll();

        return normalized_discounted_cumulative_gain<size_t>(poll, currentRankMap);
    }


    double Precision(const size_t& current, const std::map<size_t, size_t>& currentRankMap) {
        double retrieved = 0.0;


        PollMap<double, size_t> pollMap = getMap(current, currentRankMap);

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
        PollMap<double, size_t> pollMap = getMap(current, currentRankMap);
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
        double min_minmax = std::numeric_limits<double>::max();
        std::vector<size_t> minCandidates_minmax;
        double max_minmax = std::numeric_limits<double>::min();
        std::vector<size_t> maxCandidates_minmax;

        //std::string currentStirng = size_vector_to_string(current);
        //for (auto& candidates : z) {
            for (const std::pair<const unsigned long, unsigned long>& x : currentRankMap) {
                //std::string x_tmp = size_vector_to_string(x);
                if (x.first != current) {
                    double z = similarity(getVectorRepresentation(current), getVectorRepresentation(x.second));/// ???????
                    if (z < min_minmax) {
                        minCandidates_minmax.clear();
                        minCandidates_minmax.emplace_back(x.first);
                        min_minmax = z;
                    } else if (z == min_minmax) {
                        minCandidates_minmax.emplace_back(x.first);
                    }
                    if (z > max_minmax) {
                        maxCandidates_minmax.clear();
                        maxCandidates_minmax.emplace_back(x.first);
                    } else if (z == max_minmax) {
                        maxCandidates_minmax.emplace_back(x.first);
                    }
                    max_minmax = std::max(max_minmax, z);
                }
            }
        //}

        // PostProcessing
        double recallToBeZero = 0.0;
        for (const size_t& current : minCandidates_minmax) {
            //std::string current = size_vector_to_string(x);
            if (currentRankMap.find(current) != currentRankMap.end())
                recallToBeZero += 1.0;
        }
        return std::make_pair(recallToBeZero / minCandidates_minmax.size(), maxCandidates_minmax.size() / (1.0 + maxCandidates_minmax.size()));
    }
#endif

public:
    TestingGraph(Graph& ref) : passGraph{ref} {}


    void run(bool isMultithreaded, char* rootName = nullptr) {
        treeToGraphMorphism.clear();
        morphismInv.clear();
        tree.clear();
        treeIdToPathString.clear();

        // Filling up all the three data structures provided up above.
        passGraph.generateNaryTree(treeToGraphMorphism, tree, treeIdToPathString, morphismInv, rootName);
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
                std::function<ForComparison(size_t)> lam = this->getVRLambda();
                std::function<double(const ForComparison&,const ForComparison&)> exp = this->exportSimilarity();
                struct result_map maps;
                for (const size_t & candidateId : candidateSublist) {
                    //size_t current_path_length = candidateId.size();
                    //if (current_path_length != maximumHeight) continue;

                    //std::cout << "candidate " << this->passGraph.getName(candidateId) << " #" << candidateId << std::endl;
                    // Increment the path length size by one and, if not present, insert 1.0
                    std::map<size_t, size_t> rankedCandidates;
                    TestingGraphLambda<ForComparison> tgl{candidateId, rankedCandidates, lam, exp};
                    this->passGraph.isTherePath(candidateId, rankedCandidates, tgl);
                    tgl.finalize(maps);
                    //std::cout << "done" << std::endl;

                    // Generating all the candidates, represented as indices, as relevant is-a nodes for the current candidate.
                    //const std::vector<std::vector<size_t>> &candidates = generateAllPossibleSubpaths(x);
                    ///std::pair<double,double> minmax = this->minmaxCandidates(candidateId, /*splitNodeCandidatesId, */rankedCandidates);

                    // Testing the precision at k representation
                    //double currentPrecision = this->Spearman(x, candidates);
                    //if (DEBUG)  std::cout << "Current precision of " << currentPrecision << " for length " << current_path_length << std::endl;
                    //maps.spearman += this->Spearman(candidateId, rankedCandidates); // Getting the classification difference between the two
                    ///double test = this->Precision(candidateId, rankedCandidates);
                    //assert(test == 1.0);
                    ///maps.precision += test;
                    ///maps.precision_narrow += this->Precision_narrow(candidateId, rankedCandidates);
                    //maps.ncdg += this->NDCG(candidateId, rankedCandidates);
                    ///maps.recall += minmax.first;
                    ///maps.smallerNotCandidate += minmax.second;+
                    ///std::cout << maps.precision << std::endl;
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
            //update_local(precision_narrow, x);
            //update_local(ncdg, x);
            update_local(recall, x);
        }

        std::cout << "Spearman, for top-k                                              " << spearman / path_length_size << std::endl;
        //print_result_maps(path_length_size, spearman);
        std::cout << "Precision, for top-k scores                                      " << precision / path_length_size << std::endl;
        //print_result_maps(path_length_size, precision);
        //std::cout << "Precision, for top-k elements                                    " << precision_narrow / path_length_size << std::endl;
        std::cout << "Recall, for non top-k elements (wrongly matching the candidates) " << recall / path_length_size << std::endl;
        //std::cout << "More similar than worst top-1 not candidate, not current element " << smallerNotCandidate / path_length_size << std::endl;
        //print_result_maps(path_length_size, precision_narrow);
    }

protected:
    /// ???? virtual void initialize_hierarchy_with_all_paths(const std::vector<std::vector<size_t>> &subgraph_as_paths) = 0;
    virtual ForComparison getVectorRepresentation(const size_t& current) const = 0;
    std::function<ForComparison(size_t)> getVRLambda() const {
        return [this](size_t x) {
            return getVectorRepresentation(x);
        };
    }

    virtual double similarity(const ForComparison& lhs, const ForComparison& rhs) const = 0;
    std::function<double(const ForComparison&,const ForComparison&)> exportSimilarity() const {
        return [this](const ForComparison& x,const ForComparison& y) {
            return similarity(x,y);
        };
    }

    void  generateTopKCandidates(PollMap<double, size_t>& map, const size_t& current) {
        auto allVectors = getVectorRepresentation(current);
        for (auto & x : this->treeToGraphMorphism) {
            double  score = similarity(allVectors, getVectorRepresentation(x.second));
            map.add(score, x.second);
        }
    }

    /**
     * ?????
     * @param graph
     */
    virtual void passGraphDataIfRequired(const Graph& graph) = 0;
};


#endif //HIERARCHY_TESTS_TESTINGTREE_H
