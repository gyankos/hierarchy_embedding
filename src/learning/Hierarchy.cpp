//
// Created by giacomo on 29/12/19.
//

#include <queue>
#include <string_utils.h>
#include "learning/Hierarchy.h"

/**
 * TODO: the implementation of the path-finding method is incomplete.
 *   See the paper for the correct method
 */

Path Hierarchy::FindPathBetweenEntities(int entity_from, int entity_to) {
    std::set<int> common_ancestors;
    Path path{m, dim_embedding};

    FindCommonAncestors(entity_from, entity_to, common_ancestors);

    ExpandPathFromCommonAncestors(entity_from, common_ancestors,path);
    ExpandPathFromCommonAncestors(entity_to, common_ancestors,path);

    return path;
}


void Hierarchy::ExpandPathFromCommonAncestors(const int entity_idx, const std::set<int> &common_ancestors, Path &path) {

    const std::map<size_t, std::vector<size_t>>& entity_ancestor_hierarchy
            = entity_ancestor_hierarchies_[entity_idx];

    std::map<size_t , std::vector<size_t> >::const_iterator hierarchy_it;
    std::set<size_t> processed_nodes;
    int cur_idx;

    // expand from common ancestors
    std::queue<size_t> unprocessed_nodes;
    for (size_t common_ancestor : common_ancestors) {
        unprocessed_nodes.push(common_ancestor);
    }

    while (!unprocessed_nodes.empty()) {
        cur_idx = unprocessed_nodes.front();
        unprocessed_nodes.pop();
        // check if being processed, 'cause it is a DAG
        if (processed_nodes.find(cur_idx) != processed_nodes.end()) {
            continue;
        }
        processed_nodes.insert(cur_idx);
        path.AddCategoryNode(nodes_[cur_idx].id(), 1.0);

        hierarchy_it = entity_ancestor_hierarchy.find(cur_idx);
        if (hierarchy_it == entity_ancestor_hierarchy.end()) {
            continue;
        }
        const std::vector<size_t>& cur_children  = hierarchy_it->second;
        for (int c_idx_i : cur_children) {
            //const int c_idx = c_idx_i;
            unprocessed_nodes.push(c_idx_i);
        }
    } // end of expansion
}

/**
 * @param [IN] entity_from
 * @param [IN] entity_to
 * @param [OUT] common_ancestors:
 *     common ancestors that have disjoint paths to entity_from
 *     and entity_to
 *
 */
void Hierarchy::FindCommonAncestors(int entity_from, int entity_to,
                                    std::set<int>& common_ancestors) {

    const std::vector<size_t> &from = string_split_to_sizetvector(bijection.getKey(entity_from), "_");
    const std::vector<size_t> &to = string_split_to_sizetvector(bijection.getKey(entity_to), "_");
    std::vector<size_t> ancestor;
    for (size_t i = 0, n = std::min(from.size(), to.size()); i<n; i++) {
        if (from[i] == to[i]) {
            ancestor.emplace_back(from[i]);
        } else
            break;
    }
    common_ancestors.emplace(bijection.getValue(size_vector_to_string(ancestor)));
}

Hierarchy::Hierarchy(DistMetricMode m, size_t dimEmbedding, fixed_bimap<std::string, size_t>& bimap)
        : m(m), dim_embedding(dimEmbedding), bijection{bimap} {}

void Hierarchy::InitHierarchy(int num_entity_and_category) {
    //Node *node;
    // entity node
    for (int nidx = 0; nidx < num_entity_and_category; ++nidx){
        //node = new Node(nidx, -1);  // init levels of entities to -1
        nodes_.emplace_back(nidx, -1);
    }
    // category node
    for (int nidx = 0; nidx < num_entity_and_category; ++nidx){
        //node = new Node(nidx, -1);
        nodes_.emplace_back(nidx, -1);
    }
    //entity_ancestor_weights_.resize(num_entity);
    entity_ancestor_hierarchies_.resize(num_entity_and_category);
}

void Hierarchy::AddAncestors(const int entity_id, std::set<size_t> &ancestors) {
    std::map<size_t, std::vector<size_t>> ancestor_hierarchy;
    for (auto it = ancestors.begin(); it != ancestors.end(); ++it) {
        const std::vector<int>& parents = nodes_[*it].parent_idx();
        for (int p_idx : parents) {
            if (ancestors.find(p_idx) != ancestors.end()) {
                (ancestor_hierarchy)[p_idx].push_back(*it);
            }
        }
    }
    entity_ancestor_hierarchies_[entity_id] = ancestor_hierarchy;
}

