// Date: 2014.10.26

#include <time.h>
#include "learning/EEE_Engine.h"


//TODO

// Constructor
EEEngine::EEEngine(DistMetricMode metric, const size_t dim_embedding, naryTree &tree,
                   fixed_bimap<std::string, size_t> &bimap)
        : entity_category_hierarchy_{metric, dim_embedding, bimap}, treeRef{tree}, bijection{bimap} {
    num_neg_sample_ = 50;
    num_entity_category = 0;
}

EEEngine::~EEEngine() {}

void EEEngine::InitEntityCategories(size_t num_entities_and_categories) {
    for (size_t entity_idx = 0; entity_idx < num_entities_and_categories; entity_idx++) {
        //const int entity_idx = tokens[0];
        Node *entity_node = entity_category_hierarchy_.node(entity_idx);
        //for (int p_idx = 1; p_idx < tokens.size(); ++p_idx) {
        //const int category_id = tokens[p_idx];
        const int category_idx = entity_idx + num_entities_and_categories;
        // add parent categories to entity
        entity_node->AddParent(category_idx);
        // add child entity to category @hzt
        entity_category_hierarchy_.node(category_idx)->AddChild(entity_idx);
    }
    num_entity_category = num_entities_and_categories;
}


void EEEngine::ReadHierarchyIdWithLevel(naryTree &tree, size_t levelId) {
    const size_t entity_idx = tree.getId();
    const size_t category_idx = entity_idx + num_entity_category;

    // Getting the pre-allocated category information
    Node *category_node = entity_category_hierarchy_.node(category_idx);

    // Each entity has itself as an ancestor, represented as a category node.
    std::set<size_t> self_as_ancestor;
    self_as_ancestor.emplace(category_idx);
    entity_category_hierarchy_.AddAncestors(entity_idx, self_as_ancestor);

    // Iterating over all the possible children within the tree subroot
    for (size_t i = 1, n = tree.getChildNo(); i <= n; i++) {
        naryTree &next = tree.getIthChild(i - 1);
        const size_t child_category_idx = next.getId() + num_entity_category;
        category_node->AddChild(child_category_idx);
        entity_category_hierarchy_.node(child_category_idx)->AddParent(category_idx);
        category_node->set_level(levelId);

        // Continuing recursively with the other sub-tree nodes.
        ReadHierarchyIdWithLevel(next, levelId + 1); // recursive call
    }
}

