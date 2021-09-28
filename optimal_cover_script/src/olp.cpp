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

// COMPUTE ========================================================================================


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
