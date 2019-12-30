//
// Created by giacomo on 29/12/19.
//

#ifndef HIERARCHY_TESTS_HIERARCHY_H
#define HIERARCHY_TESTS_HIERARCHY_H


#ifndef ENTITY_HIERARCHY_HPP_
#define ENTITY_HIERARCHY_HPP_


#include <vector>
#include <map>
#include <set>
#include <stdint.h>
#include <utility>
#include <learning/Path.h>
#include <fixed_bimap.h>
#include "learning/Node.h"


    class Hierarchy {
    public:
        Hierarchy(DistMetricMode m, size_t dimEmbedding, fixed_bimap <std::string, size_t>& bimap);
        ~Hierarchy() {};

        Path FindPathBetweenEntities(int entity_from, int entity_to);

        void InitHierarchy(int num_entity_and_category);;


        Node* node(int idx) { return &nodes_[idx]; }

        void AddAncestors(const int entity_id, std::set<size_t> &ancestors);

        //const std::map<int, float>& entity_ancestor_weights(const int entity_idx);

        DistMetricMode m;
        size_t dim_embedding;
    private:
    // Using the bijection to determine which is the common ancestor for the two elements, eventually the root
        void FindCommonAncestors(int entity_from, int entity_to,
                                 std::set<int>& common_ancestors);

        void ExpandPathFromCommonAncestors(const int entity_idx, const std::set<int> &common_ancestors, Path &path);

    private:
        // [0, num_entity): entity nodes
        // [num_entity, num_entity + num_category): category nodes
        // Note: entity_id = index in nodes_;
        //       category_id = index in nodes_ - num_entity
        std::vector<Node> nodes_;

        // dim = num_entity
        // each entry map<int, float> is: ancestor_category_idx => weight
        //std::vector<std::map<int, float>> entity_ancestor_weights_;
        // dim = num_entity
        // each entry map<int, vector<int> > is:
        // ancestor_category_idx => child category idx of entity's ancestor
        std::vector<std::map<size_t , std::vector<size_t>>> entity_ancestor_hierarchies_;
        fixed_bimap<std::string, size_t>& bijection;

        //int num_entity_;
        //int num_category_;
        //int max_node_level_;

        //set<int>::const_iterator set_it_;
        //map<int, float>::const_iterator map_cit_;
    };


#endif



#endif //HIERARCHY_TESTS_HIERARCHY_H
