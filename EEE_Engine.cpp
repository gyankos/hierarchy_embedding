// Date: 2014.10.26

#include <string>
#include <fstream>
#include <iterator>
#include <vector>
#include <stdint.h>
#include <thread>
#include <time.h>
#include <omp.h>
#include <math.h>
#include <algorithm>
#include <cassert>
#include <sstream>
#include "EEE_Engine.h"


//TODO

    // Constructor
    EEEngine::EEEngine(DistMetricMode metric, const size_t dim_embedding) : /*train_data_{metric, dim_embedding},*/ entity_category_hierarchy_{metric, dim_embedding} {
        num_neg_sample_ = 50;
        num_train_data_ = 0;
        num_entity_category = 0;
        //num_category_ = 0;
    }

    EEEngine::~EEEngine() { }

void EEEngine::ReadEntityCategoryFile(size_t num_entities_and_categories) {
        //LOG(INFO) << "Reading " << filename;

        /*std::ifstream entity_category_file(num_entities_and_categories.c_str());
        if (!entity_category_file.is_open()) {
            //LOG(FATAL) << "fail to open:" << filename;
        }*/

        for (size_t entity_idx = 0; entity_idx < num_entities_and_categories; entity_idx++) {
            //const int entity_idx = tokens[0];
            Node* entity_node = entity_category_hierarchy_.node(entity_idx);
            //for (int p_idx = 1; p_idx < tokens.size(); ++p_idx) {
                //const int category_id = tokens[p_idx];
                const int category_idx = entity_idx + num_entities_and_categories;
                // add parent categories to entity
                entity_node->AddParent(category_idx);
                // add child entity to category @hzt
                entity_category_hierarchy_.node(category_idx)->AddChild(entity_idx);
            //}
        }

        //int counter = 0;
        //std::string line;
        //while (getline(entity_category_file, line)) {
            /*std::istringstream iss(line);
            std::vector<int> tokens(
                    (std::istream_iterator<int>(iss)), std::istream_iterator<int>());
*/
            // entity_idx in hierarchy = entity_id

            /*++counter;
            if (counter % 200000 == 0) {
                //cout << "." << std::flush;
            }*/
        //}
        //cout << endl;

        //entity_category_file.close();
        //assert(counter == num_entity_);
        num_entity_category = num_entities_and_categories;
    }


#if 0
    void EEEngine::ReadEntityAncestorFile_txt(const std::string &filename) {
        std::ifstream ancestor_file(filename.c_str());
        assert(ancestor_file.is_open());


        //LOG(INFO) << "Reading " << filename;
        //int counter = 0;
        int num_field, entity_id, ancestor_id;
        float rev_ancestor_weight;
        while (ancestor_file >> num_field) {
            ancestor_file >> entity_id;

            //std::map<int, float> ancestor_weight_map;
            for (int idx = 1; idx < num_field; ++idx) {
                ancestor_file >> ancestor_id >> rev_ancestor_weight;

                /*(ancestor_weight_map)[ancestor_id + num_entity_category]
                        = (1.0 / rev_ancestor_weight);*/
            }
            entity_category_hierarchy_.AddAncestors(
                    entity_id, ancestor_weight_map);

            /*counter++;
            if (counter % 200000 == 0) {
                //cout << "." << flush;
            }*/
        }
        //cout << endl;
        ancestor_file.close();
    }
#endif

   /* void EEEngine::BuildNoiseDistribution() {
        CHECK(entity_freq_.size() == num_entity_) << " ";

        entity_freq_[0] = pow(entity_freq_[0], 1.33333333);
        for (int e_id = 1; e_id < num_entity_; ++e_id) {
            entity_freq_[e_id] = pow(entity_freq_[e_id], 1.33333333);
            entity_freq_[e_id] += entity_freq_[e_id - 1];
        }
        freq_sum_ = entity_freq_[num_entity_ - 1];

#ifdef DEBUG
        LOG(INFO) << "BuildNoiseDistribution: sum = " << freq_sum_;
#endif
        LOG(INFO) << "Build Noise Distribution Done.";
    }*/

#if 0
    void EEEngine::ThreadCreateMinibatch(const std::vector<int> &next_minibatch_data_idx,
                                         std::vector<Datum> &next_minibatch) {
        const int batch_size = next_minibatch_data_idx.size();

        for (int d_idx = 0; d_idx < batch_size; ++d_idx) {
            const int data_idx = next_minibatch_data_idx.at(d_idx);
            Datum datum = train_data_.datum(data_idx);

            // compute path between entity_i and entity_o
            Path entity_pair_path
                    = entity_category_hierarchy_.FindPathBetweenEntities(
                            datum.entity_i(), datum.entity_o());
            datum.AddPath(entity_pair_path);

            // sample negative entities
            SampleNegEntities(datum);

            next_minibatch[d_idx] = datum;
        }
    }


    void EEEngine::SampleNegEntities(Datum &datum) {
        const int entity_i = datum.entity_i();
        const std::set<int>& pos_entities
                = train_data_.positive_entities(entity_i);

        for (int neg_sample_idx = 0; neg_sample_idx < num_neg_sample_;
             ++neg_sample_idx) {
            int neg_entity = RandSampleNegEntity();
            while (neg_entity == entity_i ||
                   pos_entities.find(neg_entity) != pos_entities.end() ||
                   train_data_.positive_entities(neg_entity).find(entity_i)
                   != train_data_.positive_entities(neg_entity).end()) {
                neg_entity = RandSampleNegEntity();
            }
            // Generate path between entity_i and neg_sample
            Path neg_path = entity_category_hierarchy_.FindPathBetweenEntities(
                    entity_i, neg_entity);

            datum.AddNegSample(neg_sample_idx, neg_entity, neg_path);
        }
    }

    inline int EEEngine::RandSampleNegEntity() {
        double rn = ((double)rand() / RAND_MAX) * freq_sum_;
        return std::upper_bound(entity_freq_.begin(), entity_freq_.end(), rn)
               - entity_freq_.begin();
        //return rand() % num_entity_;
    }
#endif
#if 0
    inline void EEEngine::CopyMinibatch(const vector<Datum*>& source,
                                        vector<Datum*>& target) {
        for (int idx = 0; idx < source.size(); ++idx) {
            target[idx] = source[idx];
        }
    }


    inline void EEEngine::ClearMinibatch(vector<Datum*>& minibatch) {
        // clear neg samples and related grads
        //for (int d_idx = 0; d_idx < minibatch.size(); ++d_idx){
        //  minibatch[d_idx]->ClearNegSamples();
        //}

        // destroy datum's
        for (int d_idx = 0; d_idx < minibatch.size(); ++d_idx){
            delete minibatch[d_idx];
        }
    }
#endif

    void EEEngine::Start() {
#if 0
        int thread_id = 0;

        entity::Context& context = entity::Context::get_instance();
        const int client_id = context.get_int32("client_id");
        const int num_client = context.get_int32("num_client");
        const int num_thread = context.get_int32("num_thread");
        const int num_iter = context.get_int32("num_iter");
        const int batch_size = context.get_int32("batch_size");
        const int eval_interval = context.get_int32("eval_interval");
        const int num_iter_per_eval = context.get_int32("num_iter_per_eval");
        const int snapshot = context.get_int32("snapshot");
        const string& output_file_prefix = context.get_string("output_file_prefix");
        const string& resume_path = context.get_string("resume_path");
        const int resume_iter = context.get_int32("resume_iter");

        int eval_counter = 0;

        BuildNoiseDistribution();

        // EEEL solver initialization
        Solver eeel_solver(num_entity_, num_category_);
        if (resume_path.length() > 0) {
            LOG(INFO) << "Resume from " << resume_path;
            CHECK_GE(resume_iter, 0);
            eeel_solver.Restore(resume_path, resume_iter);
        } else {
            eeel_solver.RandInit();
        }

        // workload manager configuration
        WorkloadManagerConfig workload_mgr_config;
        workload_mgr_config.thread_id = thread_id;
        workload_mgr_config.client_id = client_id;
        workload_mgr_config.num_clients = num_client;
        workload_mgr_config.num_threads = num_thread;
        workload_mgr_config.batch_size = batch_size;
        workload_mgr_config.num_data = num_train_data_;
        WorkloadManager workload_mgr(workload_mgr_config);

        // pre-computed minibatch
        vector<int> next_minibatch_data_idx(workload_mgr.GetBatchSize());
        vector<Datum*> next_minibatch(workload_mgr.GetBatchSize());
        // current-used minibatch
        vector<Datum*> minibatch(workload_mgr.GetBatchSize());
        // pre-computed test minibatch
        vector<int> next_test_minibatch_data_idx(workload_mgr.GetBatchSize());
        vector<Datum*> next_test_minibatch(workload_mgr.GetBatchSize());
        // current-used test minibatch
        vector<Datum*> test_minibatch(workload_mgr.GetBatchSize());
        bool test = false;
        thread minibatch_creator;


        // skip initial minibatches
        for (int skip_idx = 0; skip_idx < resume_iter; ++skip_idx) {
            workload_mgr.IncreaseDataIdxByBatchSize();
        }
        // create the first minibatch
        workload_mgr.GetBatchDataIdx(workload_mgr.GetBatchSize(),
                                     next_minibatch_data_idx);
        LOG(INFO) << "Segfault here.";
        ThreadCreateMinibatch(&next_minibatch_data_idx, &next_minibatch);
        workload_mgr.IncreaseDataIdxByBatchSize();

        clock_t t_start = clock();
        double wait_time = 0;
        double test_time = 0;

        // Train
        const int start_iter = (resume_iter > 1 ? resume_iter : 1);
        for (int iter = start_iter; iter <= num_iter; ++iter) {
            // get current minibatch from pre-computed one
            CopyMinibatch(next_minibatch, minibatch);
            // pre-compute next minibatch
            workload_mgr.GetBatchDataIdx(workload_mgr.GetBatchSize(),
                                         next_minibatch_data_idx);
            //LOG(INFO) << "start thread";
            minibatch_creator = thread(&EEEngine::ThreadCreateMinibatch, this,
                                       &next_minibatch_data_idx, &next_minibatch);
            //ThreadCreateMinibatch(&next_minibatch_data_idx, &next_minibatch);
            //LOG(INFO) << "end thread";

            // TODO
            //LOG(INFO) << workload_mgr.GetDataIdx() << " "
            //    << ((double)(clock() - t_start) / CLOCKS_PER_SEC);

            // optimize based on current minibatch
            eeel_solver.Solve(minibatch);

            //LOG(INFO) << "SOLVE done.";

            ClearMinibatch(minibatch);

            if (iter % eval_interval == 0 || iter == start_iter) {
                clock_t t_start_test = clock();

                float obj = 0;
                int cur_batch_start_idx = workload_mgr.GetDataIdx();
                // get the first test minibatch
                minibatch_creator.join();
                CopyMinibatch(next_minibatch, test_minibatch);
                for (int test_iter = 0; test_iter < num_iter_per_eval; ++test_iter) {

                    // TODO
                    //LOG(INFO) << "test " << cur_batch_start_idx << " "
                    //    << ((double)(clock() - t_start) / CLOCKS_PER_SEC);

                    if (test_iter < num_iter_per_eval - 1) {
                        // pre-compute next test minibatch
                        cur_batch_start_idx
                                = workload_mgr.GetNextBatchStartIdx(cur_batch_start_idx);
                        workload_mgr.GetBatchDataIdx(workload_mgr.GetBatchSize(),
                                                     next_test_minibatch_data_idx, cur_batch_start_idx);
                        minibatch_creator = thread(&EEEngine::ThreadCreateMinibatch, this,
                                                   &next_test_minibatch_data_idx, &next_test_minibatch);
                        //ThreadCreateMinibatch(&next_test_minibatch_data_idx, &next_test_minibatch);
                    }

                    // use current minibatch for test
                    obj += eeel_solver.ComputeObjective(test_minibatch);

                    // preserve the first test minibatch for training usage
                    if (test_iter > 0) {
                        ClearMinibatch(test_minibatch);
                    }

                    if (test_iter < num_iter_per_eval - 1) {
                        minibatch_creator.join();
                        CopyMinibatch(next_test_minibatch, test_minibatch);
                    }
                }

                obj /= (workload_mgr.GetBatchSize() * num_iter_per_eval);
                LOG(ERROR) << iter << "," << obj << ","
                           << ((double)(clock() - t_start) / CLOCKS_PER_SEC) << ","
                           << wait_time << "," << test_time; //<< "," << find_time << "," << expand_time;
                ++eval_counter;
                test = true;

                test_time += ((double)(clock() - t_start_test) / CLOCKS_PER_SEC);
                //LOG(INFO) << "test time," << test_time;
            }

            if (iter % snapshot == 0 && iter > start_iter) {
                eeel_solver.Snapshot(output_file_prefix, iter);
            }

            workload_mgr.IncreaseDataIdxByBatchSize();
            if (!test) {
                //LOG(INFO) << "wait " << wait_time;
                clock_t t_start_wait = clock();
                minibatch_creator.join();
                wait_time += ((double)(clock() - t_start_wait) / CLOCKS_PER_SEC);
            } else {
                test = false;
            }
        } // end of training

        if (num_iter % snapshot != 0) {
            eeel_solver.Snapshot(output_file_prefix, num_iter);
        }
#endif
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
    for (size_t i = 1, n = tree.getChildNo(); i<=n; i++) {
        naryTree& next = tree.getIthChild(i-1);
        const size_t child_category_idx = next.getId() + num_entity_category;
        category_node->AddChild(child_category_idx);
        entity_category_hierarchy_.node(child_category_idx)->AddParent(category_idx);
        category_node->set_level(levelId);

        // Continuing recursively with the other sub-tree nodes.
        ReadHierarchyIdWithLevel(next, levelId+1); // recursive call
    }
}

