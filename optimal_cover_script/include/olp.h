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
void copy_SA_to_SA_temp(int x, int y);
void PrintSA_temp(int x, int y);

// COMPUTE OLP_nlogn O(nlogn) Implementation ===========================================================

// vector<int> compute_OLP_nlogn(vector<set<pair<int, int>>> runs_src);
vector<int> compute_OLP_nlogn();
void remove_non_eligible_runs();
void compute_Ru(  
    int i,
    int j,
    int &sorted_i,
    int &sorted_j
    );
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
int compute_OLPi_nlogn();
// void compute_sorted_range(
//     int r1, 
//     int rm
//     );

// COMPUTE OLP Quadratic IMPROVED Implementation ========================================================

vector<int> compute_olp_improved();
int maxborder(string u, int n);
vector<int> borderarray(string text, int n);

// COMPUTE OLP Quadratic Implementation ===========================================================

// vector<int> compute_olp(vector<set<pair<int, int>>> runs_src);
vector<int> compute_olp();
vector<int> compute_rank(vector<int> &sa);
int lc(
    RMQ_succinct &,
    const std::vector<int> &,
    const std::vector<int>&,
    const int&,
    const int&
    );
// vector<set<pair<int, int>>> compute_runs(
void compute_runs(
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
