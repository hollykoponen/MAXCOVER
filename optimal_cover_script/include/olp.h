//
// Created by vi.
//

#ifndef OPTIMAL_COVERS_OLP_H
#define OPTIMAL_COVERS_OLP_H

// INCLUDES ====================================================================

#include <memory>
#include <vector>
#include <set>
#include <string>
#include "RMQ/RMQ_succinct.hpp"

// UTILS =======================================================================

int get_max(std::vector<std::vector<int>> &vec);

// COMPUTE =====================================================================

std::vector<int> compute_rank(std::vector<int>&);
int lc(
    RMQ_succinct &,
    const std::vector<int> &,
    const std::vector<int>&,
    const int&,
    const int&
);
std::vector<std::set<std::pair<int, int>>> compute_runs(
    std::vector<int>&,
    std::vector<int>&,
    std::vector<int>&,
    std::vector<int>&
);
std::vector<int> compute_r1(std::vector<int>&, std::vector<int>&);
std::vector<int> compute_rm(
    std::vector<int>&,
    std::vector<int>&
);
std::vector<int> exrun(
    const int&,
    const int&,
    const int&,
    std::vector<std::set<std::pair<int, int>>>&
);
int compute_frequency(std::vector<int>&, int&, int&);
int compute_olpi(
    int&,
    int&,
    int&,
    std::vector<int>&,
    std::vector<int>&,
    std::vector<std::set<std::pair<int, int>>>&,
    std::string&
);

#endif //OPTIMAL_COVERS_OLP_H
