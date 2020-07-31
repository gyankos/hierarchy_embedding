/*
 *  branching factor 35 could be 'decreased' by the JL lemma to at least 6.59772e+06
 tree height = 10
 maximum number of vectors per node: 2
 minimum number of vectors per node: 1
 average number of vectors per node: 1
 most frequent vector values: { 1 } with frequency 1168 which, normalized, is 0.998291
 #nodes 1170
 #candidate test leaves 879, which are { 1, 5, 7, 9, 11, 13, 15, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 45, 47, 49, 51, 53, 55, 57, 59, 62, 64, 68, 70, 72, 74, 77, 80, 83, 86, 88, 90, 94, 96, 101, 103, 105, 107, 109, 111, 113, 114, 115, 120, 122, 126, 127, 129, 131, 135, 137, 139, 143, 144, 146, 149, 151, 153, 155, 157, 158, 159, 160, 162, 164, 165, 166, 168, 171, 173, 176, 178, 180, 182, 184, 186, 187, 189, 191, 193, 195, 196, 197, 199, 201, 204, 205, 209, 211, 212, 213, 215, 216, 218, 220, 221, 222, 226, 227, 229, 231, 233, 234, 235, 236, 237, 239, 241, 242, 243, 247, 249, 250, 251, 252, 254, 255, 257, 258, 260, 264, 265, 267, 268, 271, 272, 273, 275, 277, 278, 279, 281, 283, 287, 290, 292, 294, 296, 297, 298, 299, 301, 303, 304, 306, 307, 309, 310, 312, 314, 315, 319, 321, 322, 323, 325, 327, 328, 329, 330, 332, 337, 338, 340, 341, 342, 344, 347, 348, 350, 351, 354, 356, 357, 358, 360, 362, 363, 364, 365, 366, 367, 369, 370, 372, 373, 374, 376, 378, 380, 381, 383, 384, 387, 388, 389, 391, 393, 397, 399, 401, 402, 404, 406, 407, 408, 409, 410, 412, 414, 415, 416, 418, 420, 421, 422, 425, 426, 427, 428, 429, 432, 433, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 447, 448, 449, 451, 453, 454, 455, 456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 467, 468, 469, 470, 471, 472, 473, 474, 476, 477, 478, 479, 480, 481, 483, 484, 486, 488, 489, 490, 491, 493, 494, 495, 496, 497, 499, 501, 503, 504, 505, 506, 507, 510, 514, 515, 517, 519, 520, 521, 522, 523, 524, 527, 529, 531, 533, 534, 535, 536, 537, 538, 540, 542, 543, 544, 545, 546, 547, 548, 551, 552, 553, 554, 555, 556, 557, 559, 560, 562, 563, 564, 565, 566, 567, 568, 569, 570, 571, 572, 574, 576, 577, 578, 579, 581, 583, 584, 585, 586, 587, 588, 589, 593, 595, 596, 598, 600, 601, 603, 605, 606, 607, 608, 609, 612, 613, 614, 616, 618, 619, 620, 621, 622, 623, 624, 625, 626, 627, 628, 629, 630, 631, 632, 633, 635, 636, 637, 638, 639, 640, 641, 643, 645, 646, 647, 648, 649, 650, 651, 653, 654, 655, 656, 657, 658, 659, 660, 661, 663, 664, 665, 666, 667, 668, 669, 670, 671, 672, 673, 674, 675, 676, 677, 678, 679, 680, 681, 682, 683, 685, 686, 687, 688, 690, 691, 692, 695, 696, 697, 698, 700, 701, 704, 705, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 718, 719, 720, 721, 722, 724, 725, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 737, 738, 740, 741, 742, 743, 744, 745, 746, 747, 748, 749, 751, 752, 753, 755, 757, 758, 760, 762, 763, 764, 765, 766, 767, 768, 769, 772, 774, 775, 776, 777, 778, 779, 780, 781, 782, 783, 784, 785, 786, 787, 788, 789, 790, 791, 792, 793, 794, 795, 796, 797, 798, 799, 800, 801, 802, 803, 804, 805, 806, 807, 809, 810, 811, 812, 813, 814, 816, 817, 818, 819, 820, 821, 822, 823, 824, 825, 826, 827, 828, 829, 830, 832, 833, 834, 835, 836, 837, 838, 839, 840, 841, 842, 843, 844, 845, 846, 847, 848, 849, 850, 851, 852, 853, 854, 855, 856, 857, 858, 860, 861, 862, 864, 865, 866, 867, 869, 870, 871, 872, 874, 875, 876, 877, 878, 879, 880, 881, 882, 883, 884, 886, 887, 888, 889, 890, 891, 892, 893, 894, 895, 896, 897, 899, 900, 901, 902, 903, 904, 905, 906, 907, 908, 909, 910, 911, 912, 913, 914, 915, 916, 917, 918, 919, 920, 921, 923, 924, 926, 927, 928, 929, 930, 931, 932, 933, 934, 935, 937, 938, 939, 940, 941, 942, 943, 944, 945, 946, 947, 948, 949, 950, 951, 952, 953, 954, 955, 956, 957, 958, 959, 961, 962, 963, 964, 965, 966, 967, 968, 969, 970, 972, 973, 974, 975, 976, 977, 978, 979, 980, 981, 982, 983, 984, 985, 986, 987, 988, 990, 991, 992, 993, 994, 996, 997, 998, 999, 1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019, 1020, 1021, 1022, 1023, 1024, 1025, 1026, 1027, 1028, 1029, 1030, 1031, 1032, 1033, 1034, 1035, 1036, 1037, 1038, 1039, 1040, 1041, 1042, 1043, 1044, 1045, 1046, 1047, 1048, 1049, 1050, 1051, 1052, 1053, 1054, 1055, 1056, 1057, 1058, 1059, 1060, 1061, 1062, 1063, 1064, 1065, 1066, 1067, 1068, 1069, 1070, 1071, 1072, 1073, 1074, 1075, 1076, 1077, 1078, 1079, 1080, 1081, 1082, 1083, 1084, 1085, 1086, 1087, 1088, 1089, 1090, 1091, 1092, 1093, 1094, 1095, 1096, 1097, 1098, 1099, 1100, 1101, 1102, 1103, 1104, 1105, 1106, 1107, 1108, 1109, 1110, 1111, 1112, 1113, 1114, 1115, 1116, 1117, 1118, 1119, 1120, 1121, 1122, 1123, 1124, 1125, 1126, 1127, 1128, 1129, 1130, 1131, 1132, 1133, 1134, 1135, 1136, 1137, 1138, 1139, 1140, 1141, 1142, 1143, 1144, 1145, 1146, 1147, 1148, 1149, 1150, 1151, 1152, 1153, 1154, 1155, 1156, 1157, 1158, 1159, 1160, 1161, 1162, 1163, 1164, 1165, 1166, 1167, 1168, 1169 }
Johnson (trivial) algorithm for all the possible pairs
Loading the precomputed outcome of the johnson Algorithm
...done!
... done!
Taking some time to memoize the computations...
...done!
Now Computing...
 */

#include <iostream>

#include <tests/tree/TestingTreeLearning.h>
#include <tests/tree/TestingTreeProposal.h>
#include "tests/tree/TestingTreeExternal.h"


#include <tests/TestingGraph.h>
#include <tests/tree/TestingTreePar.h>

#define POINCARRE_TREE_TESTING_BASIC_PATH       ("/media/giacomo/Data/hierarchy_paper/projects/poincare-embeddings/results/usecase_")
#define POINCARRE_GRAP_TESTING_EMBED_PATH       ("/media/giacomo/Data/hierarchy_paper/projects/poincare-embeddings/output/mammals.pth.184.txt")

/**
 * This generic function computes all the operations in function for all the possibile parameters
 * @tparam T
 * @param maxBranch
 * @param maxHeight
 * @param function
 */
template <class T> void perform_test(size_t maxBranch, size_t maxHeight, T function) {
    std::ofstream file{};
    file.open("test2.csv", std::ios::out | std::ios::app);
    for (size_t i =/* 2*/2; i<=maxBranch; i+=2) {
        for (size_t j = /*2*/2; j<=maxHeight; j++) {
            function(i, j, file);
        }
    }
}

void write_poincarre_files(size_t maximumBranchingFactor, size_t maximumHeight, std::ofstream& file) {
    std::string filename{"usecase_" + std::to_string(maximumBranchingFactor) + "_" + std::to_string(maximumHeight) + ".csv"};
    std::cerr << "write_poincarre_files: writing the current run to " << filename << std::endl;
    std::ofstream myfile (filename);
    if (myfile.is_open())
    {
        std::vector<std::vector<std::vector<size_t>> > ls = generateCompleteSubgraph(maximumBranchingFactor, maximumHeight);
        std::map<std::string, std::set<std::string>> calculateClosure;

        for (const auto& x : ls) {
            for (std::vector<size_t> y : x) {
                std::string current = size_vector_to_string(y);
                y.pop_back();
                while (!y.empty()) {
                    std::string upper = size_vector_to_string(y);
                    calculateClosure[current].insert(upper);
                    y.pop_back();
                }
                calculateClosure[current].insert("e");
            }
        }

        myfile << "id1,id2,weight" << std::endl;
        for (const auto& x : calculateClosure) {
            for (const auto& y : x.second) {
                myfile << x.first << ',' << y << ",1" << std::endl;
            }
        }
        myfile.close();
    }
}

void test_proposal(size_t maximumBranchingFactor, size_t maximumHeight, std::ofstream& file) {
    std::cout <<maximumBranchingFactor << ',' << maximumHeight << std::endl;
    file << maximumBranchingFactor << ',' << maximumHeight << ",EuclideanEmbedding,";
    std::vector<std::vector<std::vector<size_t>> > ls = generateCompleteSubgraph(maximumBranchingFactor, maximumHeight);
    //size_t maximumBranchingFactor, double distanceFactor, double decayFactor, size_t  maxHeight
    TestingTreeProposal tepee{maximumBranchingFactor, 3, (double)2.0, maximumHeight};
    file << tepee.run(ls) << std::endl;
}

void external_testing(size_t maximumBranchingFactor, size_t maximumHeight, std::ofstream& file) {
    std::cout <<maximumBranchingFactor << ',' << maximumHeight << std::endl;
    std::vector<std::vector<std::vector<size_t>> > ls = generateCompleteSubgraph(maximumBranchingFactor, maximumHeight);
    file << maximumBranchingFactor << "," << maximumHeight << ",Poincarré k=7,";
    {
        TestingTreeExternal p{maximumBranchingFactor, maximumHeight, POINCARRE_TREE_TESTING_BASIC_PATH + std::to_string(maximumBranchingFactor) + "_" + std::to_string(maximumHeight) + ".csv.7pth.txt"};
        file << p.run(ls) << std::endl;
    }
    file << maximumBranchingFactor << "," << maximumHeight << ",Poincarré k=50,";
    {
        TestingTreeExternal p{maximumBranchingFactor, maximumHeight, POINCARRE_TREE_TESTING_BASIC_PATH + std::to_string(maximumBranchingFactor) + "_" + std::to_string(maximumHeight) + ".csv.50pth.txt"};
        file << p.run(ls) << std::endl;
    }
}

void testing_learning(size_t maximumBranchingFactor, size_t maximumHeight, std::ofstream& file) {
    std::cout <<maximumBranchingFactor << ',' << maximumHeight << std::endl;
    std::vector<std::vector<std::vector<size_t>>> ls = generateCompleteSubgraph(maximumBranchingFactor, maximumHeight);

    file << maximumBranchingFactor << "," << maximumHeight << ",EEEL k=7," ;
    TestingTreeLearning testing_dimension_as_branching{maximumBranchingFactor, maximumHeight, 7};
    file << testing_dimension_as_branching.run(ls) << std::endl;

    file << maximumBranchingFactor << "," << maximumHeight << ",EEEL k=50," ;
    TestingTreeLearning testing_dimension_mid_branching_100{maximumBranchingFactor, maximumHeight,
                                                            50};
    file << testing_dimension_mid_branching_100.run(ls)<< std::endl;
}

void parhier_testing(size_t maximumBranchingFactor, size_t maximumHeight, std::ofstream& file) {
    double distanceFactor = 3;
    double decayFactor = 2;
    std::cout << maximumBranchingFactor << "," << maximumHeight << std::endl;
              auto ls = generateCompleteSubgraph(maximumBranchingFactor, maximumHeight);
    {
        file << maximumBranchingFactor << "," << maximumHeight << ",ParHier Cluster k=3,";
        TestingTreePar testing{maximumBranchingFactor, maximumHeight, 3, CLUSTER, 3.0, 2.0};
        file <<  testing.run(ls) << std::endl;
    }
    {
        file << maximumBranchingFactor << "," << maximumHeight << ",ParHier Sum k=3,";
        TestingTreePar testing{maximumBranchingFactor, maximumHeight, 3, SUM, 3.0, 2.0};
        file <<  testing.run(ls) << std::endl;
    }
}

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <JLLemma.h>
#include <tests/tree/TestingTreeBasic.h>
#include <tests/graph/TestingGraphBasic1.h>
#include <tests/graph/TestingGraphBasic2.h>
#include <tests/graph/TestingGraphBasic3.h>
#include <tests/graph/TestingGraphExternal.h>
#include <GraphLatticeFromDimensionCollection.h>
//#include <tests/graph/TestingGraphBasic1.h>
#include "Graph.h"
#include "tests/graph/TestingGraphProposal.h"

void mammals_graph_tests() {
    /// XXX: important: cannot test with Basic metric, because it takes more than 64GB of primary memory to represent that: too computational inefficient to represent

    std::cerr << "At this stage, I am assuming that you have already run the Poincarré embeddings over the Mammals dataset " << std::endl
              << "provided by the same repository. As an alternative, you can use the mammal.txt that is provided by the current" << std::endl
              << "repo. The current run assumes that the trained embeddings for the Poincarré are given in the following" << std::endl
              << "path: " << std::endl << POINCARRE_GRAP_TESTING_EMBED_PATH << std::endl << "If not, please change the macro 'POINCARRE_GRAP_TESTING_EMBED_PATH' with your value of choice" << std::endl;


    Graph g{"mammals.txt"};
    double distanceFactor = 3;
    double decayFactor = 2;
    {
        std::cout << "Relevance" << std::endl;
        TestingGraphBasic1 proposal(g);
        proposal.run(true);
    }
    {
        std::cout << "Local Density (beta=0.75) " << std::endl;
        TestingGraphBasic2 proposal(g, 0.75);
        proposal.run(true);
    }
    {
        std::cout << "Multiple Descent (beta=0.75, alpha=0.5) " << std::endl;
        TestingGraphBasic3 proposal(g, 0.75, 0.5);
        proposal.run(true);
    }
    {
        std::cout << "Multiple Descent (beta=0.75, alpha=1) " << std::endl;
        TestingGraphBasic3 proposal(g, 0.75, 1.0);
        proposal.run(true);
    }
    {
        std::cout << "External 50" << std::endl;
        TestingGraphExternal proposal(g,POINCARRE_GRAP_TESTING_EMBED_PATH);
        proposal.run(true);
    }
    {
        std::cout << "Proposed" << std::endl;
        TestingGraphProposal proposal(distanceFactor, decayFactor, g);
        proposal.run(true);
    }
}

void complete_tree_benchmarking() {

    /*std::cout << "Proposal" << std::endl;
    perform_test(7, 5, test_proposal);

    std::cout << "Poincarré" << std::endl;
    std::cerr << "Warning: with respect to the Poincarré testing, you need to first generate the trees in a format that " << std::endl <<
              "is suitable to the Poincarré python project." << std::endl;
    std::cerr << "This data is going to be generated now in your local directory." << std::endl;
    perform_test(7, 5, write_poincarre_files);
    std::cerr << "Now, I assume that you have already run the tests in the Poincarré project, and that you extracted the " << std::endl <<
              "embeddings on the .pth files generated by such project using the read_cases.py script. At this stage," << std::endl <<
              "I am assuming that all the resulting .txt files associating each node to an embedding are given in the" << std::endl <<
              "following base path:" << std::endl << POINCARRE_TREE_TESTING_BASIC_PATH << std::endl;
    std::cerr << "If that is not the case, then please change the macro named 'POINCARRE_TREE_TESTING_BASIC_PATH', recompile and re-run the code" << std::endl;
    perform_test(7, 5, external_testing);

    std::cout << "Basic" << std::endl;
    perform_test(7, 5, basic_testing);
    std::cout << "ParHier" << std::endl;
    perform_test(7, 5, parhier_testing);*/
    std::cout << "Learning" << std::endl;
    perform_test(7, 5, testing_learning);
}

/**
 * Next step, on testing lattices
 */
void test_lattice() {
    std::vector<Graph> graph_vector;
    //graph_vector.reserve(2); // TODO: quando avviene il move del grafo, succede un errore di copia che erra la creazione degli archi
    std::cout << " --- Graph 1 --- " << std::endl;
    graph_vector.emplace_back("hierarchy1.txt");
    std::cout << " --- Graph 2 --- " << std::endl;
    graph_vector.emplace_back("hierarchy2.txt");

    //graph_vector[0].print_graph();
    //return 0;
    std::cout << " --- Graph Lattice --- " << std::endl;

    // Generating the graph from the lattice
    std::string lattice_file = "hierarchy_lattice.txt";
    GraphLatticeFromDimensionCollection lattice_generator;
    lattice_generator.generate(lattice_file, graph_vector);
    Graph lattice{lattice_file};
    lattice.embedding_id = lattice_generator.embedding_id;
    //Graph lattice{graph_vector};
    //lattice.print_graph();

    graph_vector[0].johnsonAlgorithm(false);
    NaryTreeGeneration gv1_results{graph_vector[0].generateNaryTree(nullptr) };
    std::cout << gv1_results << std::endl;
    //
    double cmp0 = gv1_results.maximum_branching_factor * std::pow(gv1_results.number_of_vectors_required, 2);

    graph_vector[1].johnsonAlgorithm(false);
    NaryTreeGeneration gv2_results{graph_vector[1].generateNaryTree(nullptr)};
    std::cout << gv2_results << std::endl;
    double cmp1 = gv2_results.maximum_branching_factor * std::pow(gv2_results.number_of_vectors_required, 2);

    lattice.johnsonAlgorithm(false);
    NaryTreeGeneration worstCase{ lattice.generateNaryTree(nullptr)};
    std::cout << worstCase << std::endl;

    double cmpWC = worstCase.maximum_branching_factor * std::pow(worstCase.number_of_vectors_required, 2);

    lattice.test_lattice2(graph_vector[0], graph_vector[1]);

    std::cout << "Given that by the former test the two approaches are equivalent, the decomposed approach will make us do a speedup of at most " << cmpWC << "/" << "(" << cmp0 << "+" << cmp1 << ") = " << (cmpWC)/(cmp0+cmp1) << " when all the vectors have the same number of elements (e.g., on the leaves)" << std::endl;
}

int main() {

    complete_tree_benchmarking();

    return 0;
}
