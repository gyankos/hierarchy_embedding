//
// Created by giacomo on 30/12/19.
//

#ifndef HIERARCHY_TESTS_TESTINGTREEBASIC_H
#define HIERARCHY_TESTS_TESTINGTREEBASIC_H

#include "TestingTreeBasic1.h"
#include "TestingTreeBasic2.h"
#include "TestingTreeBasic3.h"
#include <fstream>
#include <ostream>

void testing_basic_implementation();

void basic_testing(size_t maximumBranchingFactor, size_t maximumHeight, std::ofstream &file);

#endif //HIERARCHY_TESTS_TESTINGTREEBASIC_H
