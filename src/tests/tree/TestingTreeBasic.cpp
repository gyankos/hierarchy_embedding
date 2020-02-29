//
// Created by giacomo on 30/12/19.
//

#include "tests/tree/TestingTreeBasic.h"

void testing_basic_implementation() {

    /*
     * Generating the complete tree
Tree Generation
Now Computing...
Summing up the maps...
Spearman, for top-k                                              0.879905
Precision, for top-k scores                                      0.318625
Precision, for top-k elements                                    0.407143
Recall, for non top-k elements (wrongly matching the candidates) 3.73957e-312
More similar than worst top-1 not candidate, not current element 0.5
Tree Generation
Now Computing...
Summing up the maps...
Spearman, for top-k                                              0.889868
Precision, for top-k scores                                      0.556664
Precision, for top-k elements                                    0.428571
Recall, for non top-k elements (wrongly matching the candidates) 3.73957e-312
More similar than worst top-1 not candidate, not current element 0.5
Tree Generation
Now Computing...
Summing up the maps...
Spearman, for top-k                                              0.925894
Precision, for top-k scores                                      0.726778
Precision, for top-k elements                                    0.428571
Recall, for non top-k elements (wrongly matching the candidates) 3.73957e-312
More similar than worst top-1 not candidate, not current element 0.5
Tree Generation
Now Computing...
Summing up the maps...
Spearman, for top-k                                              0.910532
Precision, for top-k scores                                      0.926283
Precision, for top-k elements                                    0.428571
Recall, for non top-k elements (wrongly matching the candidates) 3.73957e-312
More similar than worst top-1 not candidate, not current element 0.5

     */

    size_t maximumBranchingFactor = 5;
    double distanceFactor = 3;
    int decayFactor = 2;

    size_t maximumHeight = 4;

    // TODO: lightweight
    /*size_t maximumBranchingFactor = 2;
     double distanceFactor = 3;
     int decayFactor = 2;

     size_t maximumHeight = 2;*/
    std::cout << "Generating the complete tree" << std::endl;
    /*std::vector<std::future<std::vector<std::vector<size_t>>>> tmp =*/ generateCompleteSubgraph(maximumBranchingFactor, maximumHeight);
    std::vector<std::vector<std::vector<size_t>>> ls = generateCompleteSubgraph(maximumBranchingFactor, maximumHeight);
    /*for (auto& element : tmp)
        ls.emplace_back(element.get());*/

    TestingTreeBasic1 testingBasic1{maximumBranchingFactor, maximumHeight};
    testingBasic1.run(ls);

    TestingTreeBasic2 testingBasic2{maximumBranchingFactor, maximumHeight, 0.75};
    testingBasic2.run(ls);

    TestingTreeBasic3 testingBasic3{maximumBranchingFactor, maximumHeight, 0.75, 0.5};
    testingBasic3.run(ls);

    TestingTreeBasic3 testingBasic4{maximumBranchingFactor, maximumHeight, 1, 0.5};
    testingBasic4.run(ls);
}

void basic_testing(size_t maximumBranchingFactor, size_t maximumHeight, std::ofstream& file) {
    std::vector<std::vector<std::vector<size_t>>> ls = generateCompleteSubgraph(maximumBranchingFactor, maximumHeight);
    std::cout << maximumBranchingFactor << "," << maximumHeight << std::endl;

    file << maximumBranchingFactor << "," << maximumHeight << ",Relevance Vector,";
    TestingTreeBasic1 testingBasic1{maximumBranchingFactor, maximumHeight};
    file << testingBasic1.run(ls) << std::endl;

    file << maximumBranchingFactor << "," << maximumHeight << ",Local Density (beta=0.75),";
    TestingTreeBasic2 testingBasic2{maximumBranchingFactor, maximumHeight, 0.75};
    file <<  testingBasic2.run(ls)<< std::endl;

    file << maximumBranchingFactor << "," << maximumHeight << ",Multiple Descent (beta=0.75 alpha=0.5),";
    TestingTreeBasic3 testingBasic3{maximumBranchingFactor, maximumHeight, 0.75, 0.5};
    file << testingBasic3.run(ls)<< std::endl;

    file << maximumBranchingFactor << "," << maximumHeight << ",Multiple Descent (beta=1 alpha=0.5),";
    TestingTreeBasic3 testingBasic4{maximumBranchingFactor, maximumHeight, 1, 0.5};
    file << testingBasic4.run(ls)<< std::endl;
}
