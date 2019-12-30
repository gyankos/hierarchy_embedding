//
// Created by giacomo on 30/12/19.
//

#include "tests/TestingBasic.h"

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

    TestingBasic1 testingBasic1{maximumBranchingFactor, maximumHeight};
    testingBasic1.run(ls);

    TestingBasic2 testingBasic2{maximumBranchingFactor, maximumHeight, 0.75};
    testingBasic2.run(ls);

    TestingBasic3 testingBasic3{maximumBranchingFactor, maximumHeight, 0.75, 0.5};
    testingBasic3.run(ls);

    TestingBasic3 testingBasic4{maximumBranchingFactor, maximumHeight, 1, 0.5};
    testingBasic4.run(ls);
}
