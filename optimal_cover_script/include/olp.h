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

// COMPUTE OLP_nlogn O(nlogn) Implementation ===========================================================

void compute_OLP_nlogn_stack(
    int i,
    int Sorted_LI, 
    int Sorted_UI,
    std::pair<int, int> &top,
    std::vector<int> &SA,
    std::vector<int> &LCP,
    std::vector<int> &OLP_nlogn,
    std::stack<std::pair<int,int>> &stack,
    std::map<int, std::vector<int>> &runsHT,
    std::vector<std::set<std::pair<int, int>>> runs  
    );
void compute_sorted_range(
    int &Sorted_LI, 
    int &Sorted_UI, 
    int r1, 
    int rm
    );
void compute_Ru(
    int &index, 
    int &i, // (i,j) is the range pair (r1 .. rm) for u
    int j, 
    int &sorted_i, 
    int &sorted_j,
    std::vector<int> &SA,
    std::vector<int> &LCP,
    std::map<int,std::vector<int>> &runsHT,
    std::vector<std::set<std::pair<int, int>>> runs 
    );
int compute_OLP_nlogn_at_index(
    int index,
    std::vector<int> &LCP,
    std::map<int, std::vector<int>> &runsHT
    );
void sort(
    int index, 
    int i, 
    int j, 
    int prev_i, 
    int prev_j, 
    int len_u,
    std::vector<int> &SA,
    std::vector<int> &LCP,
    std::map<int, std::vector<int>> &runsHT,
    std::vector<std::set<std::pair<int, int>>> runs 
    ); 
void compute_eruns(
    int index, 
    int i, 
    int j, 
    std::vector<int> &LCP,
    std::vector<int> &SA_temp,
    std::map<int, std::vector<int>> &runsHT,
    std::vector<std::set<std::pair<int, int>>> runs 
    );
void merge_compute_eruns(
    int i1, 
    int j1, 
    int i2, 
    int j2,
    std::vector<int> &SA_temp,
    std::vector<int> &LCP,
    std::map<int, std::vector<int>> &runsHT,
    int index,
    int len_u,
    std::vector<std::set<std::pair<int, int>>> runs 
    );
std::vector<int> exrun_nlogn(
    const int &i,
    const int &j,
    const int &lcp_i,
    std::vector<std::set<std::pair<int, int>>> &runs
);

// COMPUTE OLP Quadratic Implementation ===========================================================

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
int compute_frequency(
    std::vector<int> r, 
    int i_prime, 
    int j_prime
);
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
