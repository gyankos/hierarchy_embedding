//
// Created by giacomo on 29/12/19.
//

#include "learning/Path.h"

void Path::RefreshAggrDistMetric(const std::vector<Blob> &categories) {
    aggr_dist_metric_.ClearData();
    for (int c_idx = 0; c_idx < category_nodes_.size(); ++c_idx) {
        const int category_id = category_nodes_[c_idx];
        aggr_dist_metric_.Accumulate(categories[category_id],
                                     category_node_weights_[category_id]);
    }
}

void Path::AddCategoryNode(const int category_id, const float weight) {
    if (category_node_weights_.find(category_id) ==
        category_node_weights_.end()) {
        category_nodes_.push_back(category_id);
    }
    category_node_weights_[category_id] += weight;
}

void Path::ScaleCategoryWeights(const double scale, const double scale_2) {
    for (auto it = category_node_weights_.begin(); it != category_node_weights_.end(); ++it) {
        it->second *= scale * scale_2; //weighted
    }
    scale_2_ = scale_2;
}

double Path::category_node_weight(size_t category_id) const {
    return category_node_weights_.find(category_id)->second;
}

Path::Path(const Path &x) : category_nodes_{x.category_nodes_}, category_node_weights_{x.category_node_weights_}, aggr_dist_metric_{x.aggr_dist_metric_} {
    scale_2_ = x.scale_2_;
}

Path &Path::operator=(const Path &x) {
    category_node_weights_ = x.category_node_weights_;
    category_nodes_ = x.category_nodes_;
    aggr_dist_metric_ = x.aggr_dist_metric_;
    scale_2_ = x.scale_2_;
}
