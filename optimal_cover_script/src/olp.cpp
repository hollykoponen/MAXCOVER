// INCLUDES =======================================================================================
#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <set>
#include <stack>
#include <tuple>
#include <string>
#include <map>
#include "RMQ/RMQ_succinct.hpp"
#include "olp.h"

using namespace std;

// UTILS ==========================================================================================

int get_max(std::vector<std::vector<int>> &vec) {
    return (*std::max_element(
        vec.begin(), vec.end(),
        [](const std::vector<int> &a, const std::vector<int> &b) -> int {
            return a[3] < b[3];
        }
    ))[0];
}

// COMPUTE OLP_nlogn O(nlogn) Implementation ===========================================================

// TODO: Update header file to include OLP_nlogn functions
// TODO: Update makefile to ensure it compiles with new functions
// TODO: Double check for errors in code / converting from counting from 0 instead of 1


void compute_OLP_nlogn_stack(
    int i,
    int Sorted_LI, 
    int Sorted_UI,
    std::pair<int, int> &top,
    std::vector<int> &SA,
    std::vector<int> &LCP,
    std::vector<int> &OLP_nlogn,
    std::stack<std::pair<int,int>> &stack,
    std::map<int, std::array<int, 3>> &runsHT 
    ) {
    
    while (!stack.empty()) {
        if (LCP[top.first = 0]) { OLP_nlogn[top.first] = 0; }
        else {
            compute_Ru(top.first, top.second, i-1, Sorted_LI, Sorted_UI, SA, LCP, runsHT); 
            OLP_nlogn[top.first] = compute_OLP_nlogn_at_index(top.first, LCP, runsHT);
            compute_sorted_range(Sorted_LI, Sorted_UI, top.second, i-1); 
        }
        stack.pop();
    }
}

void compute_sorted_range(
    int &Sorted_LI, 
    int &Sorted_UI, 
    int r1, 
    int rm
    ) {
    
    if (Sorted_LI == 0 && Sorted_UI == 0) {
        Sorted_LI = r1; Sorted_UI = rm;
    }
    else {
        if (Sorted_UI < r1) { Sorted_UI = rm; } // Sorted range is before an unsorted range
        else if (rm < Sorted_LI) { Sorted_LI = r1; } // Unsorted range is before the sorted range
        else { Sorted_LI = r1; Sorted_UI = rm; }
    }
}

void compute_Ru(
    int index, 
    int i, // (i,j) is the range pair (r1 .. rm) for u
    int j, 
    int sorted_i, 
    int sorted_j,
    std::vector<int> &SA,
    std::vector<int> &LCP,
    std::map<int,std::array<int, 3>> runsHT) { 
    
    int len_u = LCP[index];
    if (runsHT.size() != 0) { // Remove non-eligible runs of u from runsHT
        for (auto it = runsHT.begin(); it != runsHT.end();) { // r = (i, j, p)
            if (it->first > len_u) { 
                it = runsHT.erase(it); // |v| is actually len of u; Delete all slots with period > |v| in runsHT
            }
            else {
                it++;
            }
        }
    }
    sort(index, i, j, sorted_i, sorted_j, len_u, SA, LCP, runsHT);
}


int compute_OLP_nlogn_at_index(
    int index,
    std::vector<int> &LCP,
    std::map<int, std::array<int, 3>> runsHT
    ) {
    
    int len_u = LCP[index]; 
    int OLPi = 0;
    
    for (auto const& r : runsHT ) { // for each r in hash table r = (i, j, p); 
        if (r.second[2] < len_u) { // if p < l;
            int fru = compute_frequency(r, r.sp, r.sp + len_u - 1) // TODO: Check Typing; r.sp is the starting position of the NRE in r
            OLPi = OLPi + (len_u - r[2]) * (fru - 1)
        }
    }
    return OLPi;
}

void sort(
    int index, 
    int i, 
    int j, 
    int prev_i, 
    int prev_j, 
    int len_u,
    std::vector<int> &SA,
    std::vector<int> &LCP,
    std::map<int, std::array<int, 3>> runsHT
    ) {
    
    std::vector<int> SA_temp = std::vector<int>(j-i+1, 0);
    std::copy(SA.begin() + i, SA.begin() + j + 1, SA_temp); // copy elements SA[i..j] to SA_temp[i..j]
    
    if (prev_i != 0 
      && prev_j != 0 
      && (i < prev_i || prev_j < j)
        ) {
        
        if (i < prev_i) {
            std::sort(SA_temp.begin() + i, SA_temp.begin() + prev_i); // Sort elements of SA_temp[i..prev_i -1] in ascending order, +1 end
            compute_eruns(index, i, prev_i-1, LCP, SA_temp, runsHT);     
            }
        if (prev_j < j) {
            std::sort(SA_temp.begin() + prev_j+1, SA_temp.begin() + j + 1); // Sort elements of SA_temp[prev_j+1..j] in ascending order, +1 end
            compute_eruns(index, prev_j+1, j, LCP, SA_temp, runsHT);
            }
        if (i < prev_i) { merge_compute_eruns(i, prev_i-1, prev_i, prev_j, SA_temp, runsHT); } 
        if (prev_j < j) { merge_compute_eruns(i, prev_j, prev_j + 1, j, SA_temp, runsHT); } 
        }
}

void compute_eruns(
    int index, 
    int i, 
    int j, 
    std::vector<int> &LCP,
    std::vector<int> &SA_temp,
    std::map<int, std::array<int, 3>> runsHT
    ) {
    
    int len_u = LCP[index];
    int k = i;

    while (k < j) {
        if (SA_temp[k+1] - SA_temp[k] < len_u) { 
            std::array<int, 3> r = exrun(SA_temp[k], SA_temp[k+1] + len_u - 1); // TODO: Need to implement optimized Exrun; r = (i', j', p'); TODO: Ensure runs & variable names are called correctly; Call Compute_runs beforehand to determine the runs needed to be passed into exrun
            int fru = compute_frequency(r, SA_temp[k], SA_temp[k] + len_u - 1); 
            runsHT.insert({r[2], r}); // Add r to hashtable runsHT
            k = k + fru - 1; 
        }
        else { k++; }
    }
}


void merge_compute_eruns(
    int i1, 
    int j1, 
    int i2, 
    int j2,
    std::vector<int> &SA_temp,
    std::map<int, std::array<int, 3>> runsHT
    ) { 
    
    int len_1 = j1 - i1 + 1; // length of sublist SA_temp[i1..j1]
    int len_2 = j2 - i2 + 1; // length of sublist SA_temp[i2..j2]
    std::vector<int> L = std::vector<int>(len_1, 0) ;
    std::vector<int> R = std::vector<int> (len_1, 0);

    for (int i = 1; i <= len_1; i++){ // TODO: Check j = 0 to i < len_1 instead of 1 to i <= len_1 ? 
        L[i] = SA_temp[i1 + i - 1];
    }
    for (int j = 1; j <= len_2; j++){ // TODO: Check j = 0 to j < len_2 instead jof 1 to j <= len_2 ?
        R[j] = SA_temp[i2 + j - 1];
    }

    int i = 1; 
    int j = 1; 
    
    for (int k = i1; k <= j2; k++) { 
        if( i <= len_1 && j <= len_2) {
            if (L[i] <= R[j]) {
                SA_temp[k] = L[i];
                i++;
            }
            else {
                SA_temp[k] = R[j]; 
                j++;
            }
        }
        else {
            if (i = len_1 + 1) {
                SA_temp[k] = R[j];
                j++;
            }
            if (j = len_2 + 1) { 
                SA_temp[k] = L[i];
                i++; 
            }
        }
        if( (k > i1) 
          && (SA_temp[k] - SA_temp[k-1] < lcp_i)) { //TODO: Check: Is this lcp_i that I labelled before? or some other l? or another len?
            std::array<int, 3> r = exrun(SA_temp[k-1], SA_temp[k]+lcp_i-1); // TODO: Check exrun exists; TODO: Ensure runs & variable names are called correctly; Call Compute_runs beforehand to determine the runs needed to be passed into exrun
            runsHT.insert({r[2], r}); // Hash function to add r to hashtable runsHT
        }
    }
}

// COMPUTE OLP Quadratic Implementation ===========================================================


std::vector<int> compute_rank(std::vector<int> &sa) {
//    std::cout << "compute_rank" << std::endl;
    std::vector<int> rank(sa.size());

    for (int i = 0; i < sa.size(); ++i){
        rank[sa[i]] = i;
    }

    return rank;
}

int lc(
    RMQ_succinct &rmq,
    const std::vector<int> &lc_arr,
    const std::vector<int> &r_arr,
    const int &a,
    const int &b
) {
    return rmq.query(
        std::min(r_arr[a], r_arr[b]) + 1,
        std::max(r_arr[a], r_arr[b])
    );
}

std::vector<std::set<std::pair<int, int>>> compute_runs(
    std::vector<int> &lcp,
    std::vector<int> &lcs,
    std::vector<int> &rank,
    std::vector<int> &rev_rank
) {
    // create runs list
    std::vector<std::vector<int>> runs(lcp.size(), std::vector<int>(3));
    int size = lcp.size();
    if (lcp.size() < 120) {
        std::vector<int> tmp(120, 0);
        lcp.insert(lcp.end(), tmp.begin(), tmp.end());
    }
    if (lcs.size() < 120) {
        std::vector<int> tmp(120, 0);
        lcs.insert(lcs.end(), tmp.begin(), tmp.end());
    }

    RMQ_succinct rmq_lcp = RMQ_succinct(lcp.data(), lcp.size());
    RMQ_succinct rmq_lcs = RMQ_succinct(lcs.data(), lcs.size());

    std::vector<std::set<std::pair<int, int>>> unique_runs(size);

    int top = 0;
    for (int per = 1; per <= floor(size / 2); per++) {
//        std::cout << "(period)" << per << " <= " << floor(size / 2) << "(floor(size / 2))" << std::endl;
        int pos = per - 1;
        while (pos + per < size) {
//            std::cout << "(pos + per)" << pos + per << " < " << size << "(size)" << std::endl;
            int right = lcp[lc(rmq_lcp, lcp, rank, pos, pos + per)];
            int left = lcs[lc(rmq_lcs, lcs, rev_rank, size - pos - 1, size - pos - per - 1)];
            //std::cout << "right=" << right << "; left=" << left << "; period=" << per << std::endl;

            if (left + right > per) {
                std::pair<int, int> run_pair = std::make_pair(
                    pos - left + 1,
                    pos + per + right - 1
                );
//                runs[top][0] = pos - left + 1;
//                runs[top][1] = pos + per + right - 1;
//                runs[top][2] = per;
//                top++;
                unique_runs[per].insert(run_pair);
            }

            pos = pos + per;

        }
    }

//    std::cout << "runs calculated" << std::endl;

    if (size < 120) {
        lcp.resize(size);
        lcs.resize(size);
    }

//    std::cout << "Make UniqueMatrix" << std::endl;
//    std::vector<std::vector<std::pair<int, int>>> unique_runs(lcp.size());
//    for (auto item : runs) {
//        std::pair<int, int> run_pair = std::make_pair(item[0], item[1]);
//        if (item[2] != 0 && std::find(
//            unique_runs[item[2]].begin(), unique_runs[item[2]].end(), run_pair
//        ) == unique_runs[item[2]].end()) {
//            unique_runs[item[2]].push_back(run_pair);
//        }
//
//    }
//    std::cout << "UniqueMatrix Created" << std::endl;

    return unique_runs;
}

std::vector<int> compute_r1(
    std::vector<int> &lcp,
    std::vector<int> &rsf
) {
   std::vector<int> r1 = std::vector<int>(lcp.size());
   std::stack<std::tuple<int, int>> st;
   int lastr1 = 0;

   r1[0] = -1;
   st.push(std::make_tuple(0, 0));

    for (int i = 0; i < lcp.size(); i++) {
        /*if (rsf[i] == 0) {
            r1[i] = -1;
        } else if (rsf[i] != 0) {*/
            if ((lcp[i] != 0) && (lcp[i] == lcp[std::get<0>(st.top())]) ) {
                r1[i] = std::get<1>(st.top());
            }
            else if (lcp[i] > lcp[std::get<0>(st.top())]) {
                st.push(std::make_tuple(i, i-1));
                r1[i] = std::get<1>(st.top());
            }
            else {
                while (lcp[i] < lcp[std::get<0>(st.top())]) {
                    lastr1 = std::get<1>(st.top());
                    st.pop();
                }
                if ((lcp[i] != 0) && (lcp[i] == lcp[std::get<0>(st.top())]) ) {
                    r1[i] = std::get<1>(st.top());
                }
                else if (lcp[i] == 0) {
                    st.push(std::make_tuple(0, -1));
                    r1[i] = -1;
                }
                else if (lcp[i] > lcp[std::get<0>(st.top())]) {
                    st.push(std::make_tuple(i, lastr1));
                    r1[i] = lastr1;
                }
            }
        //}
    }
//    std::cout << "r1 = [";
//    for (int i = 0; i < r1.size(); i++) {
//        std::cout << r1[i] << ", ";
//    }
//    std::cout << "]" << std::endl;

    return r1;
}

std::vector<int> compute_rm(
    std::vector<int> &r1,
    std::vector<int> &rsf
) {
    std::vector<int> rm = std::vector<int>(rsf.size());

    for (int i = 0; i < rsf.size(); ++i) {
        if (rsf[i] != 0) {
            rm[i] = r1[i] + rsf[i] - 1;
        }
    }
    
//    std::cout << "rm = [";
//    for (int i = 0; i < rm.size(); i++) {
//        std::cout << rm[i] << ", ";
//    }
//    std::cout << "]" << std::endl;

    return rm;
}


// TODO: Check that exrun runs correctly for both quadratic and O(nlogn) implementation;
// Answer: This is not optimized version. TODO: Update later for optimization to O(nlogn) version
std::vector<int> exrun(
    const int &i,
    const int &j,
    const int &lcp_i,
    std::vector<std::set<std::pair<int, int>>> &runs
){
    std::vector<int> r = std::vector<int>(3);
    for (int p = 0; p <= ((j - i + 1) / 2); ++p) {
        if (runs[p].size() > 0) {
            for (auto run : runs[p]) {
                r[0] = std::get<0>(run);
                r[1] = std::get<1>(run);
                r[2] = p;
                // std::cout << "i=" << i << "; j=" << j << "; p=" << p << "; r[0]=" << r[0] << "; r[1]=" << r[1] << "; r[2]=" << r[2] << "; lcp[index]=" << lcp_i << std::endl;
                if (r[0] <= i && j <= r[1]) {
                    return r;
                }
            }
        }
    }
}

int compute_frequency(std::vector<int> &r, int i_prime, int j_prime) {
    int i = r[0];
    int j = r[1];
    int p = r[2];
    int abs_r = j - i + 1;
    int abs_u = j_prime - i_prime + 1;
    int abs_b_prime = abs_r % p;
    int e = floor(abs_r / p);
    int abs_x;

    if (((i_prime - i) % p) == 0) {
        abs_x = 0;
    } else {
        abs_x = (p - (i_prime - i)) % p;
    }

    int abs_y = (abs_u - abs_x) % p;
    int e_prime = (abs_u - (abs_x + abs_y)) / p;
    int delta = e - e_prime;

    if (abs_x == 0) {
        if (abs_y <= abs_b_prime) {
            return (delta + 1);
        } else {
            return delta;
        }
    } else {
        if (abs_y <= abs_b_prime) {
            return delta;
        } else {
            return delta - 1;
        }
    }
}

int compute_olpi(
    int &index,
    int &r1,
    int &rm,
    std::vector<int> &lcp,
    std::vector<int> &sa,
    // std::vector<std::vector<std::pair<int, int>>> &runs,
    std::vector<std::set<std::pair<int, int>>> &runs,
    std::string &input
) {
    int olpi = 0;
    std::vector<int> sa_temp(rm - r1 + 1, 0);

    std::copy(sa.begin() + r1, sa.begin() + rm + 1, sa_temp)    // copy elements of SA[r1..rm] to SA_temp[r1..rm]
    //for (int q = r1; q <= rm; ++q) {
    //    sa_temp[q] = sa[q];
    //}

    // sort SA_temp in ascending order using mergesort O(nlogn)
    //    mergeSort(sa_temp, r1, rm); // NOTE: why r1 and rm here
    if (r1 < rm) {
        std::sort(sa_temp.begin() + r1, sa_temp.begin() + rm + 1);
    }

//    std::cout << "i = " << index << " r1: " << r1 << " rm: " << rm << "; sa_temp = [";
//    for (auto s : sa_temp) {
//        std::cout << s << ",";
//    }
//    std::cout << " ]" << std::endl;

    //std::cout << "SA_temp=";
    /*    for (int i = 0; i < sa_temp.size(); i++) {
            std::cout << sa_temp[i];
        }*/

    int k = r1;
    while (k < rm) {
        // determines if there exists an overlap
        if ((sa_temp[k + 1] - sa_temp[k]) < lcp[index]) {
            // if the difference between r[j] and r[j+1] is less than the length of the substring
            // r = (i, j, p) // a few runs: (1, 13, 6); (5, 9, 2); (7,17,4)
            std::vector<int> r = exrun(sa_temp[k], sa_temp[k + 1] + lcp[index] - 1, lcp[index], runs);
            // the frequency of the substring occuring in the run
            int fru = compute_frequency(r, sa_temp[k], sa_temp[k] + lcp[index] - 1);

            // compute the overlap which is the freq of the substring minus one, multipled by the amount of overlap (which is length of the substring minus the period)
            olpi = olpi + (lcp[index] - r[2]) * (fru-1);
            k = k + fru-1;
        }
        else {
            k++;
        }
    }

    return olpi;
}
