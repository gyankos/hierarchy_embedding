//
// Created by giacomo on 29/12/19.
//

#ifndef HIERARCHY_TESTS_TESTINGTREE_H
#define HIERARCHY_TESTS_TESTINGTREE_H

#include <chrono>
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
#include <unordered_set>
#include "PollMap.h"
#include "stats_utils.h"
#include "multithreaded/thread_pool.h"
#include "ResultMap.h"

#define         DEBUG           (false)


template <typename ForComparison>
class TestingTree {
protected:
    size_t maximumBranchingFactor, maximumHeight;

private:


public:
    TestingTree(size_t maximumBranchingFactor, size_t maximumHeight) : maximumBranchingFactor(maximumBranchingFactor),
                                                                       maximumHeight(maximumHeight) {}




    struct result_map run(std::vector<std::vector<std::vector<size_t>>>& ls) {
        // Generating all the possible nodes for the hierarchy of choice
        //const std::vector<std::vector<size_t>> &ls = generateCompleteSubgraph(maximumBranchingFactor, maximumHeight);
        //std::cout << "Tree Generation" << std::endl;
        initialize_hierarchy_with_all_paths({{}});
        for (auto& element : ls)
            initialize_hierarchy_with_all_paths(element);
        finalizeDataIngestion();

        //double path_length_size = 0, spearman = 0,  precision_leqK = 0,  precision_narrow = 0, ncdg = 0, smallerNotCandidate = 0, recall_gtK =0;

        MultithreadWrap<struct result_map> pool{(unsigned int)ls.size(), IS_MULTITHREADED};

        //std::vector<std::future<struct result_map >> futures;
        //thread_pool pool(ls.size());

        //std::cout << "Now Computing..." << std::endl;
        std::chrono::time_point<std::chrono::system_clock> now =
                std::chrono::system_clock::now();
        // Performing the test for each node in the hierarchy, root excluded
        for (auto& y: ls) {
            pool.poolExecute
            /*futures.push_back(pool.execute*/([this, ls](const std::vector<std::vector<size_t>>& y) {
                struct result_map maps;

                for (auto& x : y) {
                    size_t current_path_length = x.size();

                    if (current_path_length != maximumHeight) continue; // For simplicity of the testing's sake, only testing the elements starting from the leaves

                    
                    //const auto VRLambda = getVRLambda();
                    //const auto sim = exportSimilarity();
                    std::map<std::string, size_t> currentRankMap;
                    std::map<double, std::set<std::string>> castor_et_pollux;
                    size_t kSplit = 0;
                    ForComparison currentv = this->getVectorRepresentation(x);

                    {
                        // All the possible paths that we consider as valid are the ones in this range
                        const std::vector<std::vector<size_t>> &candidates = generateAllPossibleSubpaths(x);
                        // Generate Expected Ranking

                        for (auto& y : candidates) {
                            std::string tmp = size_vector_to_string(y);
                            size_t expectedRank = candidates.size() - y.size();
                            currentRankMap[tmp] = expectedRank;
                            kSplit = std::max(kSplit, expectedRank);
                        }
                        auto v = this->getVectorRepresentation({});
                        double score = this->similarity(currentv, v);
                        castor_et_pollux[score].emplace(size_vector_to_string({}));
                        for (const auto& all : ls) {
                            for (const auto& toRank : all) {
                                auto v = this->getVectorRepresentation(toRank);
                                double score = this->similarity(currentv, v);
                                castor_et_pollux[score].emplace(size_vector_to_string(toRank));
                            }
                        }
                    }

                    double sumX = 0, sumY = 0;

                    size_t inferredPollCurrentPos = 1;
                    size_t overallX = 0;
                    std::vector<double> calculateX;
                    std::vector<double>    calculateY;
                    auto it = castor_et_pollux.rbegin();

                    std::unordered_set<std::string>    candidates;
                    std::unordered_set<std::string>    candidates_recall;

                    double retrieved_for_precision = 0, recallToBeZero = 0.0;
                    double precisionNumber = 0, recallNumber = 0;

                    double overK = 0;
                    while (it != castor_et_pollux.rend()) { // Iterating over the current implementation's ranked elements from the smaller to the largest
                        size_t N = it->second.size();

                        overallX += N;
                        for (const auto& x : it->second) {
                            auto it2 = currentRankMap.find(x);
                            if (inferredPollCurrentPos > kSplit) { // If I am already past the K elements, then I can only intercept the real positive values with a recall measure
                                recallNumber++;
                                if (it2 != currentRankMap.end()) {
                                    if(candidates_recall.insert(x).second)
                                        recallToBeZero+=1.0;
                                }
                            } else {
                                precisionNumber++;
                                if (it2 != currentRankMap.end()) { // If I am below the K elements, I suppose that I can use this information to reconstruct the path
                                    if(candidates.insert(x).second)
                                        retrieved_for_precision+=1.0;
                                }
                            }

                            // Determining the rank that is provided in the data structure: if it was not there, then it means that it is not a valid candidate, and hence it will be mapped in the last ranked position
                            double expectedPollCurrentPos = (it2 != currentRankMap.end()) ? ((double)it2->second) : (((double)kSplit+1));
                            calculateY.emplace_back(expectedPollCurrentPos);
                            calculateX.emplace_back(inferredPollCurrentPos);
                            sumY += expectedPollCurrentPos;
                            sumX += inferredPollCurrentPos;
                        }
                        inferredPollCurrentPos++;
                        it++;
                    }

                    // Calculating precision and recall
                    maps.precision_leqK += ((double)retrieved_for_precision) / (precisionNumber);
                    maps.recall_gtK    += (recallNumber == 0.0) ? 0.0 : ((double)recallToBeZero) / (recallNumber);

                    // Last, calculating correlation coefficient
                    double numeratore = 0.0;
                    double                 calculateXSquaredSummed = 0.0;
                    double                 calculateYSquaredSummed = 0.0;
                    sumY /= ((double)overallX);
                    sumX /= ((double)overallX);
                    for (size_t i = 0; i<overallX; i++) {
                        double tmpX = (((double) calculateX[i]) - sumX);
                        double tmpY = (((double) calculateY[i]) - sumY);
                        calculateXSquaredSummed += std::pow(tmpX, 2.0);
                        calculateYSquaredSummed += std::pow(tmpX, 2.0);
                        numeratore += (tmpX) * (tmpY);
                    }
                    maps.spearman += (numeratore > 0) ? numeratore / std::sqrt(calculateXSquaredSummed * calculateYSquaredSummed) : 0;
                    maps.path_length_size += 1;

                    /*maps.path_length_size += 1.0;
                    std::pair<double,double> minmax = this->minmaxCandidates(x, ls, candidates);
                    maps.spearman += this->Spearman(x, candidates); // Getting the classification difference between the two
                    maps.precision_leqK += this->Precision(x, candidates);*/
                }

                return maps;
            }, y);
        }

        std::chrono::time_point<std::chrono::system_clock> ende =
                std::chrono::system_clock::now();

        //std::cout << "Summing up the maps..." << std::endl;
        struct result_map mappa;
        for (auto& x : pool.foreach()) {
            mappa += x;
        }


        // floating-point duration: no duration_cast needed
        std::chrono::duration<double, std::milli> fp_ms = ende - now;
        mappa.milliseconds = fp_ms.count();
        mappa.spearman /= mappa.path_length_size;
        mappa.precision_leqK /= mappa.path_length_size;
        mappa.recall_gtK /= mappa.path_length_size;
        return mappa;
    }

    /*void
    print_result_maps(const std::map<size_t, double> &keys_with_sizes,
                      const std::map<size_t, double> &keys_with_values) const {
        for(auto it_m1 = keys_with_sizes.cbegin(), end_m1 = keys_with_sizes.cend(),
                    it_m2 = keys_with_values.cbegin(), end_m2 = keys_with_values.cend();
            it_m1 != end_m1 &&  it_m2 != end_m2; it_m1++, it_m2++)
        {
            std::cout << it_m1->first << " = " << it_m2->second / it_m1->second << std::endl;
        }
    }*/

protected:
    virtual void initialize_hierarchy_with_all_paths(const std::vector<std::vector<size_t>> &subgraph_as_paths) = 0;

    virtual ForComparison getVectorRepresentation(const std::vector<size_t>& current) = 0;
    virtual double similarity(const ForComparison& lhs, const ForComparison& rhs) = 0;
    virtual void  generateTopKCandidates(PollMap<double,std::string>& map, const std::vector<size_t>& current) = 0;
    virtual void finalizeDataIngestion() = 0;
};


#endif //HIERARCHY_TESTS_TESTINGTREE_H
