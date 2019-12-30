//
// Created by giacomo on 28/12/19.
//

#ifndef HIERARCHY_TESTS_HIERARCHYLEARNING_H
#define HIERARCHY_TESTS_HIERARCHYLEARNING_H

#include <set>
#include <string>
#include <vector>
#include <cmath>
#include <random>
#include "Blob.h"
#include "Datum.h"
#include "math_utils.h"

struct Entity {
    std::string name;
    double frequency;
    double frequency_updated;
};


/**
 * Working implementation rewritten from the source code available at: https://github.com/ZhitingHu/EEEL/issues/6
 *
 * Zhiting Hu, "Entity Hierarchy Embedding" Proceedings of the 53rd Annual Meeting of the Association for Computational Linguisticsand the 7th International Joint Conference on Natural Language Processing, pages 1292–1300,Beijing, China, July 26-31, 2015.c©2015 Association for Computational Linguistics
 *
 * URL: https://www.aclweb.org/anthology/P15-1125.pdf
 */
class HierarchyLearning {
    std::vector<Entity> entities;

    //double freq_sum_;
    double learning_rate = 0.100000001;
    //double kRandInitRange = 0.05;
    double kEpsilon = 1e-20;
    //size_t num_neg_sample = 50;
    size_t num_epoch_on_batch = 10, num_iter_on_entity = 5, num_iter_on_category = 5, dim_embedding;
    size_t num_entity, num_category;
    std::vector<Blob> entity_grads, category_grads, entitiesB, categoriesB;


    std::set<size_t> updated_entities_;
    std::set<size_t> updated_categories_;

public:
    HierarchyLearning(DistMetricMode m, size_t dimEmbedding, size_t numEntity, size_t numCategory);
    HierarchyLearning&operator=(const HierarchyLearning&) = default;
    HierarchyLearning(const HierarchyLearning&) = default;

    void  Solve_single(std::vector<Datum>& minibatch);
    double ComputeObjective_single(std::vector<Datum>& val_batch);

    std::vector<Blob> copyCategories() {
        return categoriesB;
    }

    double ComputeDist(const size_t entity_from, const size_t entity_to,const Path& path);

private:

    void ComputeEntityGradient(Datum& datum);
    void AccumulateEntityGradient(const double coeff,
                                  const Blob& dist_metric, const size_t entity_from, const size_t entity_to,
                                  Blob& grad);
    void ComputeCategoryGradient(Datum& datum);
    void AccumulateCategoryGradient(const float coeff,
                                    const int entity_from, const int entity_to, Blob& grad);
    void SGDUpdateEntity(const double lr);
    void SGDUpdateCategory(const double lr);

    void BuildNoiseDistribution();
};


#endif //HIERARCHY_TESTS_HIERARCHYLEARNING_H
