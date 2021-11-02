//
// Created by vi.
//

// INCLUDES =======================================================================================
#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <set>
#include <stack>
#include <tuple>
#include <string>
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

// COMPUTE OLP* O(nlogn) Implementation ===========================================================

// TODO: Update header file to include OLP* functions
// TODO: Update makefile to ensure it compiles with new functions
// TODO: Double check for errors in code / converting from counting from 0 instead of 1

// TODO: Ensure procedure modifies the specified variables, or does it need to be passed in by reference?
void compute_OLP*_stack(std::stack<std::pair<int,int>> stack) {
    while (!stack.empty()) do {
        if(LCP[top[0] = 0) { OLP*[top[0]] = 0 };
        else {
            compute_Ru(top[0], top[1], i-1, Sorted_LI, Sorted_UI); // TODO: Check function exists
            OLP*[top[0]] = compuite_OLP*_at_index(top[0]); // TODO: Check function exists
            compute_sorted_range(Sorted_LI, Sorted_UI, top[1], i-1); 
        }
        stack.pop();
    }
}

// TODO: Ensure procedure modifies the specified variables, or does it need to be passed in by reference?
void compute_sorted_range(int Sorted_LI, int Sorted_UI, int r1, int rm){
    if (Sorted_LI = 0 && SortedUI = 0) {
        Sorted_LI = r1; Sorted_UI = rm;
    }
    else {
        if (Sorted_UI < r1) { Sorted_UI = rm; } // Sorted range is before an unsorted range
        else if (rm < Sorted_LI) { Sorted_LI = r1; } // Unsorted range is before the sorted range
        else { Sorted_LI = r1; Sorted_UI = rm; }
    }
}

// TODO: Ensure procedure modifies the specified variables, or does it need to be passed in by reference?
void compute_Ru(int index, int i, int j, sorted_i, sorted_j){ // TODO: Add variable typing; (i,j) is the range pair (r1 .. rm) for u
    int lcp_i = LCP[index];
    if (runsHT.slots != 0) { // Remove non-eligible runs of u from runsHT
        // TODO: Delete all slots with period > |v| in runsHT
    }
    sort(index, i, j, sorted_i, sorted_j, lcp_i); // TODO: Check function exists
}

// TODO: Ensure procedure modifies the specified variables, or does it need to be passed in by reference?
int commpute_OLP*_at_index(int index){
    lcp_i = LCP[index]; OLPi = 0;
    for (int r : runsHT) { // for each r in runs HT and p < l; r = (i, j, p); 
        if (p < lcp_i) {
            fru = compute_frequency(i, j, p, r.sp, r.sp + lcp_i - 1) // TODO: Check Typing; r.sp is the starting position of the NRE in r
            OLPi = OLPi + (lcp_i - p) * (fru - 1)
        }
    }
    return OLPi;
}

// TODO: Ensure procedure modifies the specified variables, or does it need to be passed in by reference?
void sort(int index, int i, int j, prev_i, prev_j, int lcp_i){ // TODO: Add variable typing; 
    // TODO: Copy elements SA[i..j] to SA*[i..j]
    if (prev_i != 0 
      && prev_j != 0 
      && (i < prev_i || prev_j < j) {
        if (i < prev_i) {
            // TODO: Sort elements of SA*[i..prev_i -1] in ascending order
            compute_eruns(index, i, prev_i-1); // TODO: Check function exists      
        }
        if (prev_j < j) {
            // TODO: Sort elements of SA*[prev_j+1..j] in ascending order
            compute_eruns(inddex, prev_j+1, j); // TODO: Check function exists
        }
        if (i < prev_i) { merge_compute_eruns(i, prev_i-1, prev_i, prev_j); } // TODO: Check function exists
        if (prev_j < j) { merge_compute_eruns(i, prev_j, prev_j + 1, j); } // TODO: Check function exists
    }
}

// TODO: Ensure procedure modifies the specified variables, or does it need to be passed in by reference?
void compute_eruns(int index, int i, int j){
    int lcp_i = LCP[index];
    int k = i;
    while (k < j) do {
        if (SA*[k+1] - SA*[k] < lcp_i) { // TODO: Ensure has access to SA*
            r = exrun(SA*[k], SA*[k+1] + lcp_i - 1); // TODO: Need to implement Exrun; r = (i', j', p'); TODO: Ensure runs & variable names are called correctly; Call Compute_runs beforehand to determine the runs needed to be passed into exrun
            fru = compute_frequency(r[0], r[1], r[2], SA*[k], SA*[k] + lcp_i - 1); // TODO: Check function exists
            runsHT.insert(r); // Add r to hashtable runsHT // TODO: ensure this is correctly called
            k = k + fru - 1; 
        }
        else { k++; }
    }
}

// TODO: Ensure procedure modifies the specified variables, or does it need to be passed in by reference?
void merge_compute_eruns(i1, j1, i2, j2){ // TODO: Check function has access to SA*
    len_1 = j1 - i1 + 1; // length of sublist SA*[i1..j1]
    len_2 = j2 - i2 + 1; // length of sublist SA*[i2..j2]
    std::vector<int> L = std::vector<int>(len_1, 0);
    std::vector<int> R = std::vector<int>(len_1, 0);
    for (i = 1; i <= len_1; i++){ // TODO: Check j = 0 to i < len_1 instead of 1 to i <= len_1 ? 
        L[i] = SA*[i1 + i - 1];
    }
    for (j = 1; j <= len_2; j++){ // TODO: Check j = 0 to j < len_2 instead jof 1 to j <= len_2 ?
        R[j] = SA*[i2 + j - 1];
    }
    int i = 1; int j = 1; 
    for (k = i1; k <= j2; k++) { 
        if( i <= len_1 && j <= len_2) {
            if (L[i] <= R[j]) {
                SA*[k] = L[i];
                i++;
            }
            else {
                SA*[k] = R[j]; 
                j++;
            }
        }
        else {
            if (i = len_1 + 1) {
                SA*[k] = R[j];
                j++;
            }
            if (j = len_2 + 1) { 
                SA*[k] = L[i];
                i++; 
            }
        }
        if( (k > i1) 
          && (SA*[k] - SA*[k-1] < lcp_i)) { //TODO: Check: Is this lcp_i that I labelled before? or some other l? or another len?
            r = exrun(SA*[k-1], SA*[k]+lcp_i-1); // TODO: Check exrun exists; TODO: Ensure runs & variable names are called correctly; Call Compute_runs beforehand to determine the runs needed to be passed into exrun
            runsHT.insert(r); // Hash function to add r to hashtable runsHT
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
//    std::vector<std::vector<std::pair<int, int>>> &runs,
    std::vector<std::set<std::pair<int, int>>> &runs,
    std::string &input
) {
    int olpi = 0;
    std::vector<int> sa_temp(input.size());

    // copy elements of SA[r1..rm] to SA*[r1..rm]
    for (int q = r1; q <= rm; ++q) {
        sa_temp[q] = sa[q];
    }

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
