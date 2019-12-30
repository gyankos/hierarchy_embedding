cmake_minimum_required(VERSION 3.15)
project(hierarchy_tests)
find_package(Threads)
set(CMAKE_CXX_STANDARD 17)
include_directories(include)

add_executable(hierarchy_tests main.cpp src/math_utils.cpp src/learning/HierarchyLearning.cpp src/learning/Blob.cpp src/learning/Datum.cpp src/learning/Path.cpp src/proposal/Proposal.cpp include/proposal/Proposal.h src/proposal/proposal_utils.cpp include/proposal/proposal_utils.h src/naryTree.cpp include/naryTree.h src/concept_vector/ConceptVector.cpp include/concept_vector/ConceptVector.h src/cout_utils.cpp include/cout_utils.h src/string_utils.cpp include/string_utils.h include/Testing.h include/PollMap.h include/stats_utils.h src/thread_pool.cpp include/thread_pool.h src/concept_vector/TestingBasic1.cpp include/concept_vector/TestingBasic1.h src/concept_vector/TestingBasic2.cpp include/concept_vector/TestingBasic2.h src/concept_vector/TestingBasic3.cpp include/concept_vector/TestingBasic3.h src/learning/Dataset.cpp include/learning/Dataset.h src/learning/EEE_Engine.cpp include/learning/EEE_Engine.h src/learning/Hierarchy.cpp include/learning/Hierarchy.h src/learning/Node.cpp include/learning/Node.h include/fixed_bimap.h)
target_link_libraries(hierarchy_tests pthread)