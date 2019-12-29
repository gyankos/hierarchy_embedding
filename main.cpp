#include <iostream>
#include <vector>
#include <proposal/Proposal.h>
#include "cout_utils.h"
#include <concept_vector/TestingBasic1.h>
#include <concept_vector/TestingBasic2.h>
#include <concept_vector/TestingBasic3.h>


double test_my_implementation(size_t maximumBranchingFactor, size_t distanceFactor, size_t decayFactor, const std::vector<std::vector<size_t>> &ls) {
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







int main() {
    size_t maximumBranchingFactor = 5;
    double distanceFactor = 3;
    int decayFactor = 2;

    size_t maximumHeight = 4;
    std::cout << "Generating the complete tree" << std::endl;
    std::vector<std::future<std::vector<std::vector<size_t>>>> tmp = generateCompleteSubgraph(maximumBranchingFactor, maximumHeight);
    std::vector<std::vector<std::vector<size_t>>> ls;
    for (auto& element : tmp)
        ls.emplace_back(element.get());

    TestingBasic1 testingBasic1{maximumBranchingFactor, maximumHeight};
    testingBasic1.run(ls);

    TestingBasic2 testingBasic2{maximumBranchingFactor, maximumHeight, 0.75};
    testingBasic2.run(ls);

    TestingBasic3 testingBasic3{maximumBranchingFactor, maximumHeight, 0.75, 0.5};
    testingBasic3.run(ls);

    TestingBasic3 testingBasic4{maximumBranchingFactor, maximumHeight, 1, 0.5};
    testingBasic4.run(ls);
    //std::cout << generateAllPossibleSubpaths({1,2,3}) << std::endl;
    //std::cout << test_basic_4(ls) << std::endl;
    //std::cout << test_my_implementation(maximumBranchingFactor, distanceFactor, decayFactor, ls) << std::endl;
    return 0;
}