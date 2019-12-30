//
// Created by giacomo on 28/12/19.
//

#include "learning/HierarchyLearning.h"


double fastsigmoid(const double x) { return (double)0.5 * (Tanh::getInstance().fasttanh(0.5 * x) + 1.); }

void HierarchyLearning::BuildNoiseDistribution() {
    /*freq_sum_ = 1.0;
    for (int e_id = 1; e_id < entities.size(); ++e_id) {
        entities[e_id].frequency_updated = 1.0;
        entities[e_id].frequency_updated += entities[e_id - 1].frequency_updated;
    }*/
    //freq_sum_ = entities.size();
}

void HierarchyLearning::SGDUpdateCategory(const double lr) {
    auto set_it_ = updated_categories_.begin();
    for (; set_it_ != updated_categories_.end(); ++set_it_) {
        const int cate_id = *set_it_;
        categoriesB[cate_id].Accumulate(category_grads[cate_id], lr);
        categoriesB[cate_id].Rectify();
        // Clear gradients
        category_grads[cate_id].ClearData();
    }
}

void HierarchyLearning::SGDUpdateEntity(const double lr) {
    auto set_it_ = updated_entities_.begin();
    for (; set_it_ != updated_entities_.end(); ++set_it_) {
        const int entity_id = *set_it_;
        entitiesB[entity_id].Accumulate(entity_grads[entity_id], lr);
        // Projection
        entitiesB[entity_id].Normalize();
        // Clear gradients
        entity_grads[entity_id].ClearData();
    }
}

void HierarchyLearning::AccumulateEntityGradient(const double coeff, const Blob &dist_metric, const size_t entity_from,
                                                 const size_t entity_to, Blob &grad) {

    //float* grad_vec = grad.mutable_data();
    //const float* entity_from_vec = entitiesB[entity_from]->data();
    //const float* entity_to_vec = entitiesB[entity_to]->data();
    /*if (entity::Context::dist_metric_mode() == entity::Context::DIAG)*/ {
        //const float* dist_metric_mat = dist_metric->data();
        for (size_t i = 0; i < dim_embedding; ++i) {
            grad.incr_data_at(i, coeff * (2 * dist_metric.get(i))
                                 * (entitiesB[entity_from].get(i) - entitiesB[entity_to].get(i)));
        }
    }/* else if (entity::Context::dist_metric_mode() == entity::Context::EDIAG) {
            LOG(FATAL) << "do not support DistMetricMode::EDIAG right now.";
        } else if (entity::Context::dist_metric_mode() == entity::Context::FULL) {
            LOG(FATAL) << "do not support DistMetricMode::FULL right now.";
        } else {
            LOG(FATAL) << "Unkown Distance_Metric_Mode";
        }*/
}

double HierarchyLearning::ComputeDist(const size_t entity_from, const size_t entity_to, const Path &path) {
    // (x-y)^{T} M (x-y) = sum_ij (x_i - y_i) * (x_j - y_j) * M_ij
    double dist = 0;
    /*if (entity::Context::dist_metric_mode() == entity::Context::DIAG)*/ {
        //const float* dist_metric_mat = path->aggr_dist_metric()->data();
        for (int i = 0; i < dim_embedding; ++i) {
            dist += (entitiesB[entity_from].get(i) - entitiesB[entity_to].get(i))
                    * (entitiesB[entity_from].get(i) - entitiesB[entity_to].get(i))
                    * path.aggr_dist_metric().get(i);
        }
    }
    return dist;
}

HierarchyLearning::HierarchyLearning(DistMetricMode m, size_t dimEmbedding, size_t numEntity, size_t numCategory)
        : dim_embedding(dimEmbedding), num_entity(numEntity), num_category(numCategory) {

    // BuildNoiseDistribution
    for (size_t i = 0; i<num_entity; i++) {
        entity_grads.emplace_back(m, dim_embedding);
        entitiesB.emplace_back(m, dim_embedding);                                           // RandInit
    }
    for (size_t i = 0; i<num_category; i++) {
        category_grads.emplace_back(m, dim_embedding, dim_embedding);
        categoriesB.emplace_back(m, dim_embedding, dim_embedding);         // RandInit
    }

    // RandInit
    std::random_device rd;
    std::mt19937 rng_engine(rd());
    std::uniform_real_distribution<double> dis(0, 1.0);
    for (int e_idx = 0; e_idx < num_entity; ++e_idx) {
        for (int i = 0; i < dim_embedding; ++i) {
            entitiesB[e_idx].init_data_at(dis(rng_engine), i);
        }
        entitiesB[e_idx].Normalize();
    }
    for (int c_idx = 0; c_idx < num_category; ++c_idx) {
        for (int i = 0; i < dim_embedding; ++i) {
            for (int j = 0; j < dim_embedding; ++j) {
                categoriesB[c_idx].init_data_at(dis(rng_engine), i, j);
            }
        }
    }

    BuildNoiseDistribution();
}

void
HierarchyLearning::AccumulateCategoryGradient(const float coeff, const int entity_from, const int entity_to, Blob &grad) {
    for (size_t i = 0; i < dim_embedding; ++i) {
        grad.incr_data_at(i, coeff * (entitiesB[entity_from].get(i) - entitiesB[entity_to].get(i))
                             * (entitiesB[entity_from].get(i)  - entitiesB[entity_to].get(i)));
    }
}

void HierarchyLearning::Solve_single(std::vector<Datum> &minibatch) {
    float update_coeff = learning_rate / minibatch.size();
    for (size_t epoch = 0; epoch < num_epoch_on_batch; ++epoch) {
        // 1) Metric aggregation (lines 4 -- 8)
        for (int d = 0; d < minibatch.size(); ++d) {
            // Postive data metric aggregation
            minibatch[d].category_path().RefreshAggrDistMetric(categoriesB);

            // Negative samples metric aggregation
            std::vector<Path>& neg_category_paths = minibatch[d].neg_category_paths();
            for (int p_idx = 0; p_idx < neg_category_paths.size(); ++p_idx) {
                neg_category_paths[p_idx].RefreshAggrDistMetric(categoriesB);
            }
        }

        // Optimize entity embedding
        for (int iter = 0; iter < num_iter_on_entity; ++iter) {
            for (int d = 0; d < minibatch.size(); ++d) {
                Datum& datum = minibatch[d];

                if (iter > 0) {
                    datum.ClearEntityGrads();
                }

                ComputeEntityGradient(datum);

                // Accumulate entity gradients
                entity_grads[datum.entity_i()].Accumulate(
                        datum.entity_i_grad(), 1.0);
                entity_grads[datum.entity_o()].Accumulate(
                        datum.entity_o_grad(), 1.0);
                if (epoch == 0 && iter == 0) {
                    updated_entities_.insert(datum.entity_i());
                    updated_entities_.insert(datum.entity_o());
                }
                for (int neg_idx = 0; neg_idx < datum.neg_category_paths().size(); ++neg_idx) {
                    entity_grads[datum.neg_entity(neg_idx)].Accumulate(datum.neg_entity_grad(neg_idx), 1.0);
                    if (epoch == 0 && iter == 0) {
                        updated_entities_.insert(datum.neg_entity(neg_idx));
                    }
                }
            } // end of minibatch

            /// Update entity vectors
            SGDUpdateEntity(update_coeff);
        } // end of optimizing entity embedding

        //LOG(ERROR) << "optimize entity vector done.";

        // Optimize category embedding
        for (int iter = 0; iter < num_iter_on_category; ++iter) {
            for (int d = 0; d < minibatch.size(); ++d) {
                Datum& datum = minibatch[d];
                if (iter > 0) {
                    // Refresh path aggregated distance metric
                    datum.category_path().RefreshAggrDistMetric(categoriesB);
                    std::vector<Path>& neg_category_paths = datum.neg_category_paths();
                    for (int p_idx = 0; p_idx < neg_category_paths.size(); ++p_idx) {
                        neg_category_paths[p_idx].RefreshAggrDistMetric(categoriesB);
                    }

                    // Clear category grads
                    datum.ClearCategoryGrads();
                }

                ComputeCategoryGradient(datum);

                // Accumulate category gradients
                const std::vector<Blob>& datum_category_grads = datum.category_grads();
                auto map_it_ = datum.category_index().begin();
                for (; map_it_ != datum.category_index().end(); ++map_it_) {
                    category_grads[map_it_->first].Accumulate(datum_category_grads[map_it_->second], 1.0);
                    if (epoch == 0 && iter == 0) {
                        updated_categories_.insert(map_it_->first);
                    }
                }
            } // end of minibatch

            /// Update category metrics
            SGDUpdateCategory(update_coeff);
        }
    } // end of epoches

    updated_entities_.clear();
    updated_categories_.clear();
}

void HierarchyLearning::ComputeEntityGradient(Datum &datum) {
    // on e_i
    const int entity_i = datum.entity_i();
    Blob& entity_i_grad = datum.entity_i_grad();
    double coeff = fastsigmoid((-1.0) * ComputeDist(
            entity_i, datum.entity_o(), datum.category_path())) - 1.0;
    AccumulateEntityGradient(coeff, datum.category_path().aggr_dist_metric(),
                             entity_i, datum.entity_o(), entity_i_grad);
    // on e_o
    // = (-1) * gradient_on_e_i, so simply do the copy
    datum.entity_o_grad().CopyFrom(entity_i_grad, -1.0);

    // neg_samples
    for (int neg_idx = 0; neg_idx < datum.neg_category_paths().size(); ++neg_idx) {
        Blob& neg_entity_grad = datum.neg_entity_grad(neg_idx);
        coeff = 1.0 - fastsigmoid(ComputeDist(
                entity_i, datum.neg_entity(neg_idx),
                datum.neg_category_path(neg_idx)));

        AccumulateEntityGradient(
                (-1.0) * coeff, datum.neg_category_path(neg_idx).aggr_dist_metric(),
                entity_i, datum.neg_entity(neg_idx), neg_entity_grad);
    }
    // Accumulate (-1) * gradient_on_neg_samples to gradient_on_e_i
    for (int neg_idx = 0; neg_idx < datum.neg_category_paths().size(); ++neg_idx) {
        entity_i_grad.Accumulate(datum.neg_entity_grad(neg_idx), -1.0);
    }
}

void HierarchyLearning::ComputeCategoryGradient(Datum &datum) {
    // process (e_i,  e_o)
    Path& category_path = datum.category_path();
    const std::vector<size_t >& category_nodes = category_path.category_nodes();
    const size_t entity_i = datum.entity_i();
    const size_t entity_o = datum.entity_o();
    double coeff = fastsigmoid((-1.0) * ComputeDist(
            entity_i, entity_o, datum.category_path())) - 1.0;
    for (int c_idx = 0; c_idx < category_nodes.size(); ++c_idx) {
        // use weighted coeff
        const float weighted_coeff = coeff
                                     * category_path.category_node_weight(category_nodes[c_idx]);
        AccumulateCategoryGradient(weighted_coeff, entity_i, entity_o,
                                   datum.category_grad(category_nodes[c_idx]));
    }

    // process (e_i, negative samples)
    for (int path_idx = 0; path_idx < datum.neg_category_paths().size(); ++path_idx) {
        Path &neg_path = datum.neg_category_path(path_idx);
        const std::vector<size_t>& neg_category_nodes = neg_path.category_nodes();
        const int neg_entity = datum.neg_entity(path_idx);
        float coeff = 1.0 - fastsigmoid(ComputeDist(
                entity_i, neg_entity, neg_path));
        for (int c_idx = 0; c_idx < neg_category_nodes.size(); ++c_idx) {
            // use weighted coeff
            const float weighted_coeff = coeff
                                         * neg_path.category_node_weight(neg_category_nodes[c_idx]);
            AccumulateCategoryGradient(weighted_coeff, entity_i, neg_entity,
                                       datum.category_grad(neg_category_nodes[c_idx]));
        }
    }
}

double HierarchyLearning::ComputeObjective_single(std::vector<Datum> &val_batch) {
    double obj = 0;
    for (int d = 0; d < val_batch.size(); ++d) {
        Datum& datum = val_batch[d];
        datum.category_path().RefreshAggrDistMetric(categoriesB);
        double datum_obj = 0;
        int entity_i = datum.entity_i();
        datum_obj += log(std::max(kEpsilon, fastsigmoid((-1.0) * ComputeDist(
                entity_i, datum.entity_o(), datum.category_path()))));
        for (int neg_idx = 0; neg_idx < datum.neg_category_paths().size(); ++neg_idx) {
            datum.neg_category_path(neg_idx).RefreshAggrDistMetric(categoriesB);
            datum_obj += log(std::max(kEpsilon, fastsigmoid(ComputeDist(
                    entity_i, datum.neg_entity(neg_idx),
                    datum.neg_category_path(neg_idx)))));
        }
        obj += datum_obj;
    }
    return obj;
}
