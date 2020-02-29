//
// Created by giacomo on 29/12/19.
//

/*
 *  branching factor 35 could be 'decreased' by the JL lemma to at least 6.59772e+06
 tree height = 10
 maximum number of vectors per node: 2
 minimum number of vectors per node: 1
 average number of vectors per node: 1.00769
 most frequent vector values: { 1 } with frequency 1159 which, normalized, is 0.990598
 #nodes 1170
 #candidate test leaves 879, which are { 1, 5, 7, 9, 11, 13, 15, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 45, 47, 49, 51, 53, 55, 57, 59, 62, 64, 68, 70, 72, 74, 77, 80, 83, 86, 88, 90, 94, 96, 101, 103, 105, 107, 109, 111, 113, 114, 115, 120, 122, 126, 127, 129, 131, 135, 137, 139, 143, 144, 146, 149, 151, 153, 155, 157, 158, 159, 160, 162, 164, 165, 166, 168, 171, 173, 176, 178, 180, 182, 184, 186, 187, 189, 191, 193, 195, 196, 197, 199, 201, 204, 205, 209, 211, 212, 213, 215, 216, 218, 220, 221, 222, 226, 227, 229, 231, 233, 234, 235, 236, 237, 239, 241, 242, 243, 247, 249, 250, 251, 252, 254, 255, 257, 258, 260, 264, 265, 267, 268, 271, 272, 273, 275, 277, 278, 279, 281, 283, 287, 290, 292, 294, 296, 297, 298, 299, 301, 303, 304, 306, 307, 309, 310, 312, 314, 315, 319, 321, 322, 323, 325, 327, 328, 329, 330, 332, 337, 338, 340, 341, 342, 344, 347, 348, 350, 351, 354, 356, 357, 358, 360, 362, 363, 364, 365, 366, 367, 369, 370, 372, 373, 374, 376, 378, 380, 381, 383, 384, 387, 388, 389, 391, 393, 397, 399, 401, 402, 404, 406, 407, 408, 409, 410, 412, 414, 415, 416, 418, 420, 421, 422, 425, 426, 427, 428, 429, 432, 433, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 447, 448, 449, 451, 453, 454, 455, 456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 467, 468, 469, 470, 471, 472, 473, 474, 476, 477, 478, 479, 480, 481, 483, 484, 486, 488, 489, 490, 491, 493, 494, 495, 496, 497, 499, 501, 503, 504, 505, 506, 507, 510, 514, 515, 517, 519, 520, 521, 522, 523, 524, 527, 529, 531, 533, 534, 535, 536, 537, 538, 540, 542, 543, 544, 545, 546, 547, 548, 551, 552, 553, 554, 555, 556, 557, 559, 560, 562, 563, 564, 565, 566, 567, 568, 569, 570, 571, 572, 574, 576, 577, 578, 579, 581, 583, 584, 585, 586, 587, 588, 589, 593, 595, 596, 598, 600, 601, 603, 605, 606, 607, 608, 609, 612, 613, 614, 616, 618, 619, 620, 621, 622, 623, 624, 625, 626, 627, 628, 629, 630, 631, 632, 633, 635, 636, 637, 638, 639, 640, 641, 643, 645, 646, 647, 648, 649, 650, 651, 653, 654, 655, 656, 657, 658, 659, 660, 661, 663, 664, 665, 666, 667, 668, 669, 670, 671, 672, 673, 674, 675, 676, 677, 678, 679, 680, 681, 682, 683, 685, 686, 687, 688, 690, 691, 692, 695, 696, 697, 698, 700, 701, 704, 705, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 718, 719, 720, 721, 722, 724, 725, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 737, 738, 740, 741, 742, 743, 744, 745, 746, 747, 748, 749, 751, 752, 753, 755, 757, 758, 760, 762, 763, 764, 765, 766, 767, 768, 769, 772, 774, 775, 776, 777, 778, 779, 780, 781, 782, 783, 784, 785, 786, 787, 788, 789, 790, 791, 792, 793, 794, 795, 796, 797, 798, 799, 800, 801, 802, 803, 804, 805, 806, 807, 809, 810, 811, 812, 813, 814, 816, 817, 818, 819, 820, 821, 822, 823, 824, 825, 826, 827, 828, 829, 830, 832, 833, 834, 835, 836, 837, 838, 839, 840, 841, 842, 843, 844, 845, 846, 847, 848, 849, 850, 851, 852, 853, 854, 855, 856, 857, 858, 860, 861, 862, 864, 865, 866, 867, 869, 870, 871, 872, 874, 875, 876, 877, 878, 879, 880, 881, 882, 883, 884, 886, 887, 888, 889, 890, 891, 892, 893, 894, 895, 896, 897, 899, 900, 901, 902, 903, 904, 905, 906, 907, 908, 909, 910, 911, 912, 913, 914, 915, 916, 917, 918, 919, 920, 921, 923, 924, 926, 927, 928, 929, 930, 931, 932, 933, 934, 935, 937, 938, 939, 940, 941, 942, 943, 944, 945, 946, 947, 948, 949, 950, 951, 952, 953, 954, 955, 956, 957, 958, 959, 961, 962, 963, 964, 965, 966, 967, 968, 969, 970, 972, 973, 974, 975, 976, 977, 978, 979, 980, 981, 982, 983, 984, 985, 986, 987, 988, 990, 991, 992, 993, 994, 996, 997, 998, 999, 1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019, 1020, 1021, 1022, 1023, 1024, 1025, 1026, 1027, 1028, 1029, 1030, 1031, 1032, 1033, 1034, 1035, 1036, 1037, 1038, 1039, 1040, 1041, 1042, 1043, 1044, 1045, 1046, 1047, 1048, 1049, 1050, 1051, 1052, 1053, 1054, 1055, 1056, 1057, 1058, 1059, 1060, 1061, 1062, 1063, 1064, 1065, 1066, 1067, 1068, 1069, 1070, 1071, 1072, 1073, 1074, 1075, 1076, 1077, 1078, 1079, 1080, 1081, 1082, 1083, 1084, 1085, 1086, 1087, 1088, 1089, 1090, 1091, 1092, 1093, 1094, 1095, 1096, 1097, 1098, 1099, 1100, 1101, 1102, 1103, 1104, 1105, 1106, 1107, 1108, 1109, 1110, 1111, 1112, 1113, 1114, 1115, 1116, 1117, 1118, 1119, 1120, 1121, 1122, 1123, 1124, 1125, 1126, 1127, 1128, 1129, 1130, 1131, 1132, 1133, 1134, 1135, 1136, 1137, 1138, 1139, 1140, 1141, 1142, 1143, 1144, 1145, 1146, 1147, 1148, 1149, 1150, 1151, 1152, 1153, 1154, 1155, 1156, 1157, 1158, 1159, 1160, 1161, 1162, 1163, 1164, 1165, 1166, 1167, 1168, 1169 }

 */

#ifndef HIERARCHY_TESTS_TESTINGGRAPH_H
#define HIERARCHY_TESTS_TESTINGGRAPH_H

#include <vector>
#include <string>
#include <proposal/proposal_utils.h>
#include <string_utils.h>
#include <cout_utils.h>
#include <map>
#include <set>
#include <iostream>
#include <multithreaded/MultithreadWrap.h>
#include <cassert>
#include "PollMap.h"
#include "stats_utils.h"
#include "multithreaded/thread_pool.h"
#include <Graph.h>
#include "TestingTree.h"
#include "TestingGraphLambda.h"
#include <chrono>

#define         DEBUG           (false)

template<typename T, typename iterableT>
std::vector<std::vector<T>> SplitVector(const iterableT& vec, size_t n)
{
    std::vector<std::vector<T>> outVec;

    size_t length = vec.size() / n;
    size_t remain = vec.size() % n;

    size_t begin = 0;
    size_t end = 0;
    auto it = vec.begin();

    for (size_t i = 0; i < std::min(n, vec.size()); ++i) {
        end += (remain > 0) ? (length + !!(remain--)) : length;
        std::vector<T> element;
        while (begin != end) {
            element.emplace_back(*it++);
            begin++;
        }
        outVec.push_back(element);
        begin = end;
    }

    return outVec;
}
template <typename ForComparison> class TestingGraph {
protected:
    size_t maximumBranchingFactor, maximumHeight;
    Graph& passGraph;

    /// Correspondence between the node in adj and the node in the graph, so that we can reconstruct all the vectors to be associated to it
    std::unordered_map<size_t, size_t> treeToGraphMorphism;
    std::unordered_map<size_t, std::unordered_set<size_t>> morphismInv;

    /// Tree representation of the data structure, with new element ids
    std::unordered_map<size_t, std::unordered_set<size_t>> tree;

    /// path string generated from the adj representation. This can be used to then generate multiple possible elements
    std::map<size_t, std::string> treeIdToPathString;

public:
    TestingGraph(Graph& ref) : passGraph{ref} {}


    void run(bool isMultithreaded, char* rootName = nullptr) {
        treeToGraphMorphism.clear();
        morphismInv.clear();
        tree.clear();
        treeIdToPathString.clear();

        // Filling up all the three data structures provided up above.
        passGraph.generateNaryTree(treeToGraphMorphism, tree, treeIdToPathString, morphismInv, rootName);
        maximumBranchingFactor = passGraph.maxBranch;
        maximumHeight = passGraph.height;
        passGraphDataIfRequired(passGraph);

        double path_length_size = 0, spearman = 0,  precision_leqK = 0,  precision_narrow = 0, ncdg = 0, smallerNotCandidate = 0, recall_gtK =0;

        size_t threadNo = (unsigned int)std::thread::hardware_concurrency();
        MultithreadWrap<struct result_map> pool{(unsigned int)std::thread::hardware_concurrency(), isMultithreaded};
        std::vector<std::vector<size_t>> splitNodeCandidatesId = SplitVector<size_t, std::set<size_t>>(passGraph.getCandidates(), threadNo);

        // Performing the test for each node in the hierarchy, root excluded

        std::chrono::time_point<std::chrono::system_clock> now =
                std::chrono::system_clock::now();
        for (const std::vector<size_t>& y: splitNodeCandidatesId) {
            pool.poolExecute([this, splitNodeCandidatesId](const std::vector<size_t> & candidateSublist) {
                std::function<ForComparison(size_t)> lam = this->getVRLambda();
                std::function<double(const ForComparison&,const ForComparison&)> exp = this->exportSimilarity();
                struct result_map maps;
                for (const size_t & candidateId : candidateSublist) {
                    std::map<size_t, size_t> rankedCandidates;
                    TestingGraphLambda<ForComparison> tgl{candidateId, rankedCandidates, lam, exp};
                    int maxLength;
                    this->passGraph.isTherePath(candidateId, rankedCandidates, tgl,maxLength);
                    tgl.finalize(maps, maxLength);
                }

                return maps;
            }, y);
        }
        std::chrono::time_point<std::chrono::system_clock> ende =
                std::chrono::system_clock::now();
        std::chrono::duration<double, std::milli> fp_ms = ende - now;

        std::cout << "Summing up the maps..." << std::endl;
        for (auto& x : pool.foreach()) {
            //auto x = y.get();
#define update_local(field, currFuture)           (field) += currFuture .  field
            update_local(path_length_size, x);
            update_local(spearman, x);
            update_local(precision_leqK, x);
            update_local(recall_gtK, x);
        }

        std::cout << "Spearman, for top-k                                              " << spearman / path_length_size << std::endl;
        std::cout << "Precision, for top-k scores                                      " << precision_leqK / path_length_size << std::endl;
        std::cout << "Recall, for non top-k elements (wrongly matching the candidates) " << recall_gtK / path_length_size << std::endl;
        std::cout << "Time (milliseconds) " << fp_ms.count() / path_length_size << std::endl;
    }

protected:
    virtual ForComparison getVectorRepresentation(const size_t& current) const = 0;
    std::function<ForComparison(size_t)> getVRLambda() const {
        return [this](size_t x) {
            return this->getVectorRepresentation(x);
        };
    }

    virtual double similarity(const ForComparison& lhs, const ForComparison& rhs) const = 0;
    std::function<double(const ForComparison&,const ForComparison&)> exportSimilarity() const {
        return [this](const ForComparison& x,const ForComparison& y) {
            return this->similarity(x,y);
        };
    }

    /**
     * Method used to precompute all the vectorial representations
     *
     * @param graph
     */
    virtual void passGraphDataIfRequired(const Graph& graph) = 0;
};


#endif //HIERARCHY_TESTS_TESTINGTREE_H
