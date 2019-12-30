//
// Created by giacomo on 30/12/19.
//

#include "tests/TestingLearning.h"

TestingLearning::TestingLearning(size_t maximumBranchingFactor, size_t maximumHeight, size_t vectorDimension)
        : Testing(maximumBranchingFactor,
                  maximumHeight), dimEmbedding(vectorDimension <= 0 ? std::max((size_t)2, maximumBranchingFactor) : vectorDimension) {
    ee_engine_ = nullptr;
    trainer = nullptr;
}

void TestingLearning::initialize_hierarchy_with_all_paths(const std::vector<std::vector<size_t>> &element) {
    for (auto& x : element/*.get()*/) {
        std::string x_val = size_vector_to_string(x);
        allNodes.emplace(x_val);
        for (auto& y : generateAllPossibleSubpaths(x)) {
            std::string y_val = size_vector_to_string(y);
            positiveExamples[x_val].emplace(y_val);
            positiveExamples[y_val].emplace(x_val);
        }
        // Getting the id of the new element that is now added: id
        num_entities_and_classes += tree.addChild(x, id);
        // Adding the id of the associated element, and associating that to the string. That is going ot be used to determine the best common ancestor.
        bimap.put(x_val, id);
    }
}

void TestingLearning::finalizeDataIngestion() {
    // Init category hierarchy
    std::cout << "Initializing the EEEL Engine with some tweaks" << std::endl;
    ee_engine_ = new EEEngine(FULL, dimEmbedding, tree, bimap, num_entities_and_classes);
    //EEEngine engine{FULL, dimEmbedding, tree, bimap, num_entities_and_classes};

    std::cout << "Generating the negative examples" << std::endl;
    std::map<std::string, std::set<std::string>> negativeExamples;
    for (const std::string& x : allNodes) {
        const std::set<std::string>& set = positiveExamples[x];
        std::set<std::string>& difference = negativeExamples[x];
        std::set_difference(set.begin(), set.end(), allNodes.begin(), allNodes.end(), std::inserter(difference, difference.begin()));
    }

    std::vector<Datum> batch = ee_engine_->generateData(positiveExamples, negativeExamples);

    // Finished to generate the
    std::cout << "Initializing the EEEL Trainer" << std::endl;
    trainer = new HierarchyLearning (FULL, dimEmbedding, num_entities_and_classes, num_entities_and_classes);

    //batch_learning(batch, trainer, 10);
    //for (int i = 0; i<10; i++)
    trainer->Solve_single(batch);
    //std::cout << trainer.ComputeObjective_single(batch) << std::endl;
}

size_t TestingLearning::getVectorRepresentation(const std::vector<size_t> &current) {
    return bimap.getValue(size_vector_to_string(current));
}

double TestingLearning::similarity(const size_t &lhs, const size_t &rhs) {
    Path path = ee_engine_->entity_category_hierarchy().FindPathBetweenEntities(lhs, rhs);
    path.RefreshAggrDistMetric(trainer->copyCategories());
    return trainer->ComputeDist(lhs, rhs, path);
}

void TestingLearning::generateTopKCandidates(PollMap<double, std::string> &map, const std::vector<size_t> &current) {
    size_t  entity = getVectorRepresentation(current);
    for (size_t e_id = 0; e_id < num_entities_and_classes; ++e_id) {
        double dist = similarity(entity, e_id);
        map.add(dist, bimap.getKey(e_id));
    }
}

void testing_learning_method() {

#if 0
    size_t maximumBranchingFactor = 5;
    double distanceFactor = 3;
    int decayFactor = 2;
    size_t maximumHeight = 3;

    std::cout << "Generating the complete tree" << std::endl;
    naryTree tree{0};
    size_t id = 0, num_entities_and_classes = 1; //num_entities_and_classes is initialized with 1, because the root will be added automatically

    std::cout << "Generating the positive examples" << std::endl;
    std::map<std::string, std::set<std::string>> positiveExamples;
    fixed_bimap<std::string, size_t>             bimap;             // Bijection between path as a string and the id
    std::set<std::string> allNodes;                                 // Containing all the nodes within the hierarchy

    // Adding the root node
    bimap.put("", 0);
    for (auto& element : generateCompleteSubgraph(maximumBranchingFactor, maximumHeight)) {
        for (auto& x : element/*.get()*/) {
            std::string x_val = size_vector_to_string(x);
            allNodes.emplace(x_val);
            for (auto& y : generateAllPossibleSubpaths(x)) {
                std::string y_val = size_vector_to_string(y);
                positiveExamples[x_val].emplace(y_val);
                positiveExamples[y_val].emplace(x_val);
            }
            // Getting the id of the new element that is now added: id
            num_entities_and_classes += tree.addChild(x, id);
            // Adding the id of the associated element, and associating that to the string. That is going ot be used to determine the best common ancestor.
            bimap.put(x_val, id);
        }
    }

    // Init category hierarchy
    std::cout << "Initializing the EEEL Engine with some tweaks" << std::endl;
    EEEngine engine{FULL, maximumBranchingFactor, tree, bimap, num_entities_and_classes};

    std::cout << "Generating the negative examples" << std::endl;
    std::map<std::string, std::set<std::string>> negativeExamples;
    for (const std::string& x : allNodes) {
        const std::set<std::string>& set = positiveExamples[x];
        std::set<std::string>& difference = negativeExamples[x];
        std::set_difference(set.begin(), set.end(), allNodes.begin(), allNodes.end(), std::inserter(difference, difference.begin()));
    }

    std::vector<Datum> batch = engine.generateData(positiveExamples, negativeExamples);

    // Finished to generate the
    std::cout << "Initializing the EEEL Trainer" << std::endl;
    HierarchyLearning trainer{FULL, maximumBranchingFactor, num_entities_and_classes, num_entities_and_classes};

    //batch_learning(batch, trainer, 10);
    for (int i = 0; i<10; i++)
        trainer.Solve_single(batch);
    std::cout << trainer.ComputeObjective_single(batch) << std::endl;
#endif
    size_t maximumBranchingFactor = 5;
    double distanceFactor = 3;
    int decayFactor = 2;

    size_t maximumHeight = 4;

    // TODO: lightweight
    /*size_t maximumBranchingFactor = 2;
     double distanceFactor = 3;
     int decayFactor = 2;

     size_t maximumHeight = 2;*/
    std::vector<std::vector<std::vector<size_t>>> ls = generateCompleteSubgraph(maximumBranchingFactor, maximumHeight);
    TestingLearning testing_dimension_as_branching{maximumBranchingFactor, maximumHeight, maximumBranchingFactor};
    testing_dimension_as_branching.run(ls);

    TestingLearning testing_dimension_mid_branching_100{maximumBranchingFactor, maximumHeight,
                                                        (maximumBranchingFactor + 100) / 2};
    testing_dimension_mid_branching_100.run(ls);

    TestingLearning testing_dimension_100{maximumBranchingFactor, maximumHeight, 100};
    testing_dimension_100.run(ls);
}

void batch_learning(const std::vector<Datum> &batch, HierarchyLearning &trainer, int iterations) {
    for (int i = iterations; i > 0; i--) {
        // split vector into sub-vectors each of size n
        size_t n = std::thread::hardware_concurrency() * i;

        // determine number of sub-vectors of size n
        int size = (batch.size() - 1) / n + 1;

        // create array of vectors to store the sub-vectors
        std::vector<Datum> vec[size];

        // each iteration of this loop process next set of n elements
        // and store it in a vector at k'th index in vec
        for (int k = 0; k < size; ++k)
        {
            // get range for next set of n elements
            auto start_itr = std::next(batch.cbegin(), k*n);
            auto end_itr = std::next(batch.cbegin(), k*n + n);

            // allocate memory for the sub-vector
            vec[k].resize(n);

            // code to handle the last sub-vector as it might
            // contain less elements
            if (k*n + n > batch.size()) {
                end_itr = batch.cend();
                vec[k].resize(batch.size() - k*n);
            }

            // copy elements from the input range to the sub-vector
            std::copy(start_itr, end_itr, vec[k].begin());
        }

        for (size_t i = 0; i<size; i++) {
            trainer.Solve_single(vec[i]);
        }

        double  objFinal = 0.0;
        for (size_t i = 0; i<size; i++) {
            objFinal += trainer.ComputeObjective_single(vec[i]);
        }
        std::cout << objFinal << std::endl;

    }
}
