#ifndef OPTIMAL_COVERS_OLP_H
#define OPTIMAL_COVERS_OLP_H

// INCLUDES ====================================================================

#include "globals.h"
#include <memory>
#include "RMQ/RMQ_succinct.hpp"

using namespace std;

// UTILS =======================================================================

int get_max(vector<vector<int>> &vec);
void reset_SA_temp();
void set_OLP_nlogn_variables(
    int &lcp_i_source, 
    int top_index_source
    );

// COMPUTE OLP_nlogn O(nlogn) Implementation ===========================================================

vector<int> compute_OLP_nlogn(vector<set<pair<int, int>>> runs_src);
void remove_non_eligible_runs();
void eligible_run(  
    int i,
    int j,
    int &sorted_i,
    int &sorted_j
    );
void copy_sort_compute(int x, int y);
void compute_eruns(
    int i, 
    int j
    );
void merge_compute_eruns(
    int i1, 
    int j1, 
    int i2, 
    int j2
    );
int compute_OLP_nlogn_at_index();
// void compute_sorted_range(
//     int r1, 
//     int rm
//     );

// COMPUTE OLP Quadratic Implementation ===========================================================

void runs_for_exrun();
vector<int> compute_olp(vector<set<pair<int, int>>> runs_src);
vector<int> compute_rank(vector<int> &sa);
int lc(
    RMQ_succinct &,
    const std::vector<int> &,
    const std::vector<int>&,
    const int&,
    const int&
    );
vector<set<pair<int, int>>> compute_runs(
    vector<int> &lcp,
    vector<int> &lcs,
    vector<int> &rank,
    vector<int> &rev_rank
    );
void compute_R1();
void compute_RM();
vector<int> exrun(
    const int &i,
    const int &j,
    const int &lcp_i
    );
int compute_frequency(
    vector<int> r, 
    int i_prime, 
    int j_prime
    );
int compute_olpi(
    int &index,
    int &r1,
    int &rm
    ); 

#endif //OPTIMAL_COVERS_OLP_H
