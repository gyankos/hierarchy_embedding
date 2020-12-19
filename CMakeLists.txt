cmake_minimum_required(VERSION 3.15)
project(hierarchy_tests)

find_package(Threads)                                       #For multithreading
find_package (Eigen3 3.3 REQUIRED NO_MODULE)                #For par2hier paper using Truncated SVD

set(CMAKE_CXX_STANDARD 17)
include_directories(include)

add_subdirectory(submodules/lemon)
include_directories(submodules/lemon)
include_directories(cmake-build-debug/submodules/lemon)

add_library(EuclideanEmbedding src/Graph.cpp src/naryTree.cpp src/cout_utils.cpp)

add_executable(hierarchy_tests main.cpp src/JLLemma.cpp src/math_utils.cpp src/learning/HierarchyLearning.cpp src/learning/Blob.cpp src/learning/Datum.cpp src/learning/Path.cpp src/proposal/Proposal.cpp include/proposal/Proposal.h src/proposal/proposal_utils.cpp include/proposal/proposal_utils.h src/naryTree.cpp include/naryTree.h src/concept_vector/ConceptVector.cpp include/concept_vector/ConceptVector.h src/cout_utils.cpp include/cout_utils.h src/string_utils.cpp include/string_utils.h include/tests/TestingTree.h include/PollMap.h include/stats_utils.h src/multithreaded/thread_pool.cpp include/multithreaded/thread_pool.h src/tests/tree/TestingTreeBasic1.cpp include/tests/tree/TestingTreeBasic1.h src/tests/tree/TestingTreeBasic2.cpp include/tests/tree/TestingTreeBasic2.h src/tests/tree/TestingTreeBasic3.cpp include/tests/tree/TestingTreeBasic3.h src/learning/Dataset.cpp include/learning/Dataset.h src/learning/EEE_Engine.cpp include/learning/EEE_Engine.h src/learning/Hierarchy.cpp include/learning/Hierarchy.h src/learning/Node.cpp include/learning/Node.h include/fixed_bimap.h include/multithreaded/MultithreadWrap.h src/tests/tree/TestingTreeProposal.cpp include/tests/tree/TestingTreeProposal.h src/tests/tree/TestingTreeLearning.cpp include/tests/tree/TestingTreeLearning.h src/tests/tree/TestingTreeBasic.cpp src/tests/tree/TestingTreeExternal.cpp include/tests/tree/TestingTreeExternal.h src/par2hier/Par2Hier.cpp include/par2hier/Par2Hier.h src/JLLemma.cpp include/JLLemma.h src/Graph.cpp include/Graph.h src/tests/tree/TestingTreePar.cpp include/tests/tree/TestingTreePar.h src/tests/graph/TestingGraphProposal.cpp src/tests/graph/TestingGraphBasic1.cpp src/tests/graph/TestingGraphBasic2.cpp src/tests/graph/TestingGraphBasic3.cpp include/tests/graph/TestingGraphProposal.h src/tests/graph/TestingGraphExternal.cpp include/tests/graph/TestingGraphExternal.h src/GraphLatticeFromDimensionCollection.cpp src/ExportGraphEmbeddings.cpp include/ExportGraphEmbeddings.h)
target_link_libraries(hierarchy_tests pthread Eigen3::Eigen lemon)