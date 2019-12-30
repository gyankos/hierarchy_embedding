//
// Created by giacomo on 30/12/19.
//

#include "tests/TestingProposal.h"

TestingProposal::TestingProposal(size_t maximumBranchingFactor, double distanceFactor, double decayFactor,
                                 size_t maxHeight) : Testing{maximumBranchingFactor, maxHeight}, prop{(double)maximumBranchingFactor, distanceFactor, decayFactor} {}

void TestingProposal::initialize_hierarchy_with_all_paths(const std::vector<std::vector<size_t>> &subgraph_as_paths) {
    allPossiblePaths.insert(allPossiblePaths.end(), subgraph_as_paths.begin(), subgraph_as_paths.end()); // Copying the possible paths
    //allPossiblePaths = subgraph_as_paths;
    return;
}

std::vector<size_t> TestingProposal::getVectorRepresentation(const std::vector<size_t> &current) {
    return current;
}

double TestingProposal::similarity(const std::vector<size_t> &lhs, const std::vector<size_t> &rhs) {
    double distanceMetric =  prop.rankingMetric(lhs, rhs);
    double normalizedDistance = distanceMetric / (1+distanceMetric);
    double similairity = 1 - normalizedDistance;
    //std::cerr << "sim(" << size_vector_to_string(lhs) << ", " << size_vector_to_string(rhs)<< ") = " << similairity << " via distance " << distanceMetric << " normalized as " << normalizedDistance << std::endl;
    return similairity;
}

void TestingProposal::generateTopKCandidates(PollMap<double, std::string> &map, const std::vector<size_t> &current) {
    std::string currentString{size_vector_to_string(current)};
    for (auto & it : allPossiblePaths) {
        //if (it.first != currentString) {
        double  score = similarity(current, it);
        // std::cout << "sim(" << currentString << "," << it.first << ")=" << score << std::endl ;
        map.add(score, size_vector_to_string(it));
        //}
    }
    //map.printNarrow();
}

double test_my_implementation(size_t maximumBranchingFactor, size_t distanceFactor, size_t decayFactor,
                              const std::vector<std::vector<size_t>> &ls) {
    Proposal hd(maximumBranchingFactor, distanceFactor, decayFactor);
    std::cout << "Starting the testing..." << std::endl;
    double tests = 0.0;
    double correct = 0.0;
    for (int i = 0, n = ls.size(); i<n; i++) {
        std::cout << i << std::endl;
        auto veci = ls[i];
        for (int j = 0; j<i; j++) {
            auto vecj = ls[j];
            bool consistenceInference = hd.isConsistent(veci, vecj);
            tests++;
            if ((subArrayOf(veci, vecj)) || (subArrayOf(vecj, veci))) {
                if (!consistenceInference)
                    std::cerr << ("ERROR! the algorithm says that is not consistent, while it should be") << std::endl;
                else
                    correct++;
            } else {
                if (consistenceInference) {
                    std::cerr << veci << " vs. " << vecj << std::endl;
                    bool consistenceInference = hd.isConsistent(veci, vecj);
                    std::cerr << ("ERROR! the algorithm says that is consistent, while it is not") << std::endl;
                } else
                    correct++;
            }
        }
    }
    return (correct/tests);
}

void test_my_proposal() {
    //testing_basic_implementation(); return 0;
    /*
 * Spearman, for top-k                                              0.99999
Precision, for top-k scores                                      0.00640205
Precision, for top-k elements                                    1
Recall, for non top-k elements (wrongly matching the candidates) 3.72461e-312
More similar than worst top-1 not candidate, not current element 0.5
 */

    size_t maximumBranchingFactor = 5;
    double distanceFactor = 3;
    int decayFactor = 2;

    size_t maximumHeight = 4;
    std::cout << "Generating the complete tree" << std::endl;
    std::vector<std::vector<std::vector<size_t>> > ls = generateCompleteSubgraph(maximumBranchingFactor, maximumHeight);

    //size_t maximumBranchingFactor, double distanceFactor, double decayFactor, size_t  maxHeight
    TestingProposal tepee{maximumBranchingFactor, distanceFactor, (double)decayFactor, maximumHeight};
    tepee.run(ls);

}
