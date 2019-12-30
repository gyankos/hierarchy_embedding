//
// Created by giacomo on 28/12/19.
//

#include "learning/Datum.h"

Datum::Datum(DistMetricMode metric, const size_t dim_embedding, const size_t entity_i, const size_t entity_o,
             const size_t num_neg_sample, const size_t count)
        : m{metric}, entity_i_(entity_i), dim_embedding{dim_embedding}, category_path_{metric, dim_embedding},
                                         entity_o_(entity_o), count_(count), entity_i_grad_{metric, dim_embedding}, entity_o_grad_{metric, dim_embedding} {
    neg_entity_id_.resize(num_neg_sample);
    //neg_entity_grads_.resize(num_neg_sample);
    for (size_t i = 0; i < num_neg_sample; ++i) {
        neg_entity_grads_.emplace_back(metric, dim_embedding);
    }
}

void Datum::AddPath(Path &path) {
    category_path_ = path;
    const std::vector<size_t>& category_path_nodes = category_path_.category_nodes();
    for (size_t c_idx = 0; c_idx < category_path_nodes.size(); ++c_idx) {
        category_index_[category_path_nodes[c_idx]] = c_idx;
        category_grads_.emplace_back(m, dim_embedding, dim_embedding);
    }
}

Path &Datum::category_path() {
    return category_path_;
}

void Datum::AddNegSample(size_t neg_idx, size_t neg_entity_id, Path &path) {
    neg_entity_id_[neg_idx] = neg_entity_id;
    neg_category_paths_.push_back(path);

    const std::vector<size_t>& neg_category_path_nodes = path.category_nodes();
    for (size_t c_idx = 0; c_idx < neg_category_path_nodes.size(); ++c_idx) {
        const size_t category_id = neg_category_path_nodes[c_idx];
        if (category_index_.find(category_id) == category_index_.end()) {
            category_index_[category_id] = category_grads_.size();
            category_grads_.emplace_back(m, dim_embedding, dim_embedding);
        }
    }
}

void Datum::ClearEntityGrads() {
    entity_i_grad_.ClearData();
    entity_o_grad_.ClearData();
    for (size_t neg_idx = 0; neg_idx < neg_entity_grads_.size(); ++neg_idx) {
        neg_entity_grads_[neg_idx].ClearData();
    }
}

void Datum::ClearCategoryGrads() {
    for (size_t c_idx = 0; c_idx < category_grads_.size(); ++c_idx) {
        category_grads_[c_idx].ClearData();
    }
}
