// INCLUDES =======================================================================================

#include "globals.h"
#include <cmath>
#include <tuple>
#include "RMQ/RMQ_succinct.hpp"
#include "olp.h"

using namespace std;

// Global Variables to compute OLP Definition ======================================================


vector<set<pair<int, int>>> runs;
map<int, std::vector<int>> runsHT;

int sorted_i = 0; 
int sorted_j = 0;  
int lcp_i = 0;
vector<int> SA_temp(0);
int top_index = 0;

// UTILS ==========================================================================================

int get_max(vector<vector<int>> &vec) {
    return (*max_element(
        vec.begin(), vec.end(),
        [](const vector<int> &a, const vector<int> &b) -> int {
            return a[3] < b[3];
        }
    ))[0];
}

void reset_SA_temp(){
    vector<int> SA_star(input_size, 0);
    SA_temp = SA_star;
}

void set_OLP_nlogn_variables(
    int &lcp_i_src, 
    int top_index_src
    ){
    
    lcp_i = lcp_i_src;
    top_index = top_index_src;
}

#define printvar(x) cout << (#x) << " = " << x << endl;

void PrintArrays(vector<int> arr){
    for(auto a : arr){
        cout << a << ",";
    } cout << endl; cout << endl;
}

// COMPUTE OLP_nlogn O(nlogn) Implementation ===========================================================

/* compute OLP_nlogn array of string */
/* O(nlogn) implementation */
vector<int> compute_OLP_nlogn(vector<set<pair<int, int>>> runs_src) {

    runs = runs_src;
    vector<int> OLP_nlogn = vector<int>(input_size, 0);
    OLP_nlogn[1] = 0;

    stack<pair<int,int>> st; // stack of pairs (index, r1)
    pair<int,int> top; // (index, r1)

    runs_for_exrun();

    int sorted_i = 0; // lower index of sorted range in SA_star
    int sorted_j = 0;  // upper index of sorted range in SA_star

    int i = 1;
    while (i < input_size) {
        if (st.empty()) { st.push( make_pair(i,i-1)); }
        else {
            top = st.top();
            if (LCP[top.first] < LCP[i]) { st.push( make_pair(i, i-1) );}
            else if (LCP[top.first] == LCP[i]) { OLP_nlogn[i] = 0; }
            else{
                if (LCP[i] <= 1){
                    runsHT.clear();  // clears the contents of hash table; original: runsHT.slots = 0
                    sorted_i = 0; 
                    sorted_j = 0;
                    reset_SA_temp(); cout << "We reset_SA_temp" << endl;
                }
                while ( !st.empty() && LCP[top.first] > LCP[i]) {
                    top = st.top();
                    OLP_nlogn[top.first] = 0;
                    if(LCP[top.first] != 1) {
                        set_OLP_nlogn_variables(LCP[top.first], top.first);
                        eligible_run(top.second, i-1, sorted_i, sorted_j); 
                        OLP_nlogn[top.first] = compute_OLP_nlogn_at_index();
                    }
                    st.pop();
                }
                if (!st.empty() && top.first == LCP[i]) { OLP_nlogn[i] = 0; }
                else { st.push( make_pair(i, top.first-1) ); }
            }
        }
        i++;
    }
    while (!st.empty()) {
        if (LCP[top.first] == 0) { OLP_nlogn[top.first] = 0; }
        else {
            eligible_run(top.second, i-1, sorted_i, sorted_j); 
            OLP_nlogn[top.first] = compute_OLP_nlogn_at_index();
        }
        st.pop();
    }
    return OLP_nlogn;

}

void remove_non_eligible_runs() { 

    // Remove non-eligible runs of u from runsHT
    if (runsHT.size() != 0) { 
        for (auto it = runsHT.begin(); it != runsHT.end();) { // r = (i, j, p)
            if (it->first > lcp_i) { 
                it = runsHT.erase(it); // |v| is actually len of u; Delete all slots with period > |v| in runsHT
            }
            else {
                it++;
            }
        }
    }
}


void eligible_run(  
    int i, // (i,j) is the range pair (r1 .. RM) for u
    int j,
    int &sorted_i,
    int &sorted_j
    ) {

    // (sorted_i,sorted_j) is the sorted range pair of the previously computed (r1 .. RM) for u;  
 
    cout << "\neligible_run\n" << endl;

    printvar(i);printvar(j);printvar(sorted_i);printvar(sorted_j);

    remove_non_eligible_runs();

    PrintArrays(SA_temp);

    cout << endl;

    if(sorted_i == 0 || sorted_j == 0){ // either sorted_i or sorted_j are 0
        cout << "Hi 1" << endl;
        copy_sort_compute(i,j);     
        sorted_i = i;
        sorted_j = j;
    }
    else{
        if(i < sorted_i){ // [i..sorted_i..sorted_j..j]
            if(sorted_j < j){
                cout << "Hi 2" << endl;

                copy_sort_compute(i, sorted_i-1); // Compute left: [i..sorted_i)
                copy_sort_compute(sorted_j+1, j); // Compute right: (sorted_j ..j]

                PrintArrays(SA_temp);
                
                merge_compute_eruns(i, sorted_i-1, sorted_i, sorted_j); // Merge left to mid: [i..sorted_j)
                merge_compute_eruns(i, sorted_j, sorted_j+1, j); // Merge left-mid to right: [i..j]        
                sorted_i = i;
                sorted_j = j; // [i..j]
                printvar(sorted_i);printvar(sorted_j);
            }
            else{
                if(sorted_i <= j){ // [i..sorted_i..j..sorted_j]
                    cout << "Hi 3" << endl;
                    copy_sort_compute(i, sorted_i -1); 

                PrintArrays(SA_temp);
                
                    merge_compute_eruns(i, sorted_i-1,sorted_i,sorted_j); 
                    sorted_i = i; // [i..sorted_j]
                    printvar(sorted_i);printvar(sorted_j);
                }
                else{ // [i..j..sorted_i..sorted_j]
                    cout << "Hi 4" << endl;
                    copy_sort_compute(i, j);  
                    
                    PrintArrays(SA_temp);
                    
                     if ( (sorted_i - j) > 1 ) {printvar(sorted_i-j);sorted_j = j;} // if gap b/t [j..sorted_i] then range is [i..j]
                    // else the range is [i..sorted_j]
                    sorted_i = i;
                    printvar(sorted_i);printvar(sorted_j);
                }
            }
        }
        else{
            if(sorted_j < i){ // [sorted_i..sorted_j..i..j]
                cout << "Hi 5" << endl;
                copy_sort_compute(i,j);   

                PrintArrays(SA_temp);
                
                 if ( (i - sorted_j) > 1) {printvar(i-sorted_j); sorted_i = i;} // if gap b/t [sorted_j..i] then range is [i..j]
                // else the range is [sorted_i..j]
                sorted_j = j;
                printvar(sorted_i);printvar(sorted_j);
            }   
            else{
                if(sorted_j < j){ // [sorted_i..i..sorted_j..j]
                    cout << "Hi 6" << endl;
                    copy_sort_compute(sorted_j+1,j); 

                PrintArrays(SA_temp);
                
                    merge_compute_eruns(sorted_i,sorted_j,sorted_j+1,j); 
                    sorted_j = j;
                    printvar(sorted_i);printvar(sorted_j);
                }
                else{ // [sorted_i..i..j..sorted_j]
                    // Do Nothing
                    cout << "Hi 7" << endl;
                }
            } 
        }
    }
}

void copy_sort_compute(int x, int y){
    // Copy (unsorted) elements SA[x..y] to SA_temp[x..y]
    for (int i = x; i <= y; i++){ 
        SA_temp[i] = SA[i]; 
    } 

    sort(SA_temp.begin() + x, SA_temp.begin() + y+1); // Sort elements in SA* in ascending order

    compute_eruns(x, y);    // Compute eligible runs     
}

void compute_eruns(
    int i, 
    int j
    ) {
    
    cout << endl; cout << "Compute_eruns"<<endl;

    int k = i;
    while (k < j) {
        if (SA_temp[k+1] - SA_temp[k] < lcp_i) {    
            cout << "exrun input = " << SA_temp[k] << ", " << (SA_temp[k+1] + lcp_i-1) << ", " <<lcp_i << endl;     
            vector<int> r = exrun(SA_temp[k], SA_temp[k+1] + lcp_i - 1, lcp_i); 

            cout << "run = " << r[0] << "," << r[1] << "," << r[2] << endl;
            int fru = compute_frequency(r, SA_temp[k], SA_temp[k] + lcp_i-1); 
            r.push_back(SA_temp[top_index]);
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
    int j2
    ) { 

    cout << endl; cout << "merge_compute_eruns start" << endl;
    cout << "i1,j1,i2,j2 = " << i1 <<"," << j1 <<"," << i2 << "," << j2 << endl;

    int len_1 = j1 - i1 + 1; // length of sublist SA_temp[i1..j1]
    int len_2 = j2 - i2 + 1; // length of sublist SA_temp[i2..j2]

    int i = 0; 
    int j = 0; 
    
    // std::vector<int> L = std::vector<int>(len_1, 0) ;
    // std::vector<int> R = std::vector<int> (len_1, 0);
    
    // for (int i = 1; i <= len_1; i++){ 
    //     L[i] = SA_temp[i1 + i - 1];
    // }
    // for (int j = 1; j <= len_2; j++){ 
    //     R[j] = SA_temp[i2 + j - 1];
    // }

    sort(SA_temp.begin() +i1, SA_temp.begin() +j2+1); // TODO: Ensure this is a linear time operation O(M) and not MlogM where M = j2-i1+1

    PrintArrays(SA_temp);

    for (int k = i1; k <= j2; k++) { 
        // if( i <= len_1 && j <= len_2) {
        //     if (L[i] <= R[j]) {
        //         SA_temp[k] = L[i];
        //         i++;
        //     }
        //     else {
        //         SA_temp[k] = R[j]; 
        //         j++;
        //     }
        // }
        // else {
        //     if (i = len_1 + 1) {
        //         SA_temp[k] = R[j];
        //         j++;
        //     }
        //     if (j = len_2 + 1) { 
        //         SA_temp[k] = L[i];
        //         i++; 
        //     }
        // }

        if( (k > i1) && (SA_temp[k] - SA_temp[k-1] < lcp_i)) { 

            cout << "exrun inputs: " << endl;
            printvar(SA_temp[k-1]); printvar(SA_temp[k] + lcp_i - 1); printvar(lcp_i);
            cout << endl;

            vector<int> r = exrun(SA_temp[k-1], SA_temp[k] + lcp_i - 1, lcp_i); 

            r.push_back(SA_temp[top_index]);
            
            runsHT.insert({r[2], r}); // Hash function to add r to hashtable runsHT
        }
    }

}

/* Compute_OLP at index functions*/

int compute_OLP_nlogn_at_index() {

    int OLPi = 0;

    for (auto const& r : runsHT ) { // for each r in hash table r = (i, j, p); 
        if (r.first < lcp_i) { // if p < l;

            vector<int> run = {r.second[0], r.second[1], r.second[2]};

            int fru = compute_frequency(run, r.second[3], r.second[3] + lcp_i- 1); 

            OLPi = OLPi + (lcp_i - run[2]) * (fru - 1);

        }
    }

    return OLPi;

}

// void compute_sorted_range(
//     int r1, 
//     int rm
//     ) {
    
//     if (sorted_i == 0 && sorted_j == 0) {
//         sorted_i = r1; sorted_j = rm;
//     }
//     else {
//         if (sorted_j < r1) { sorted_j = rm; } // Sorted range is before an unsorted range
//         else if (rm < sorted_i) { sorted_i = r1; } // Unsorted range is before the sorted range
//         else { sorted_i = r1; sorted_j = rm; }
//     }
// }

// // TODO: Check that exrun runs correctly for O(nlogn) implementation;
// // Answer: This is optimized version. TODO: Update to O(n) version
// vector<int> exrun_nlogn(
//     const int &i,
//     const int &j,
//     const int &lcp_i,
//     vector<set<pair<int, int>>> &runs // TODO: change to hashtable
//     ) {
//     vector<int> r = vector<int>(3);
//     if (p = (j-i+1)/2) {
//         if (runs[p].size() > 0) {
//             for (auto run : runs[p]) {
//                 if (get<0>(run) <= i && j <= get<1>(run)) {
//                     return r;
//                 }
//             }
//         }
//     }
// }

// COMPUTE OLP Quadratic Implementation ===========================================================


/* compute OLP array of string */
/* QUADRATIC implementation */
vector<int> compute_olp(vector<set<pair<int, int>>> runs_src) {

    runs = runs_src;
    vector<int> olp(input_size,0);
    reset_SA_temp();
    runs_for_exrun();
    OLP.reserve(input_size);
    for (int i = 0; i < input_size; i++) {
        if (RSF[i] != 0) {
            olp[i] = compute_olpi(i, R1[i], RM[i]);
        }
    }

    return olp;
}

vector<int> compute_rank(vector<int> &sa) {
    vector<int> rank(sa.size());

    for (int i = 0; i < input_size; ++i){
        rank[sa[i]] = i;
    }

    return rank;
}

int lc(
    RMQ_succinct &RMq,
    const vector<int> &lc_arr,
    const vector<int> &r_arr,
    const int &a,
    const int &b
    ) {

    return RMq.query(
        min(r_arr[a], r_arr[b]) + 1,
        max(r_arr[a], r_arr[b])
    );
}

vector<set<pair<int, int>>> compute_runs(
    vector<int> &lcp,
    vector<int> &lcs,
    vector<int> &rank,
    vector<int> &rev_rank
    ) {

    // create runs list
    vector<vector<int>> runs(input_size, vector<int>(3));

    if (input_size < 120) {
        vector<int> tmp(120, 0);
        lcp.insert(lcp.end(), tmp.begin(), tmp.end());
        lcs.insert(lcs.end(), tmp.begin(), tmp.end());
    }

    RMQ_succinct RMq_lcp = RMQ_succinct(LCP.data(), lcp.size());
    RMQ_succinct RMq_lcs = RMQ_succinct(lcs.data(), lcs.size());

    vector<set<pair<int, int>>> unique_runs(input_size); // TODO: This should be a hashtable instead of a vector b/c resizing using .insert() is inefficient, which requires reallocation of vector size. Hashtable addresses this with chaining.

    int top = 0;
    for (int per = 1; per <= floor(input_size / 2); per++) {
        int pos = per - 1;
        while (pos + per < input_size) {
            int right = lcp[lc(RMq_lcp, lcp, rank, pos, pos + per)];
            int left = lcs[lc(RMq_lcs, lcs, rev_rank, input_size - pos - 1, input_size - pos - per - 1)];

            if (left + right > per) {
                pair<int, int> run_pair = make_pair(
                    pos - left + 1,
                    pos + per + right - 1
                );
                unique_runs[per].insert(run_pair); // TODO: This should be a hashtable instead of a vector b/c resizing using .insert() is inefficient, which requires reallocation of vector size. Hashtable addresses this with chaining.
            }

            pos = pos + per;

        }
    }

    if (input_size < 120) {
        lcp.resize(input_size);
        lcs.resize(input_size);
    }
    return unique_runs;
}

void compute_R1() {
   R1.reserve(input_size);
   stack<tuple<int, int>> st;
   int lastr1 = 0;

   R1[0] = -1;
   st.push(make_tuple(0, 0));

    for (int i = 0; i < input_size; i++) {
        if ((LCP[i] != 0) && (LCP[i] == LCP[get<0>(st.top())]) ) {
            R1[i] = get<1>(st.top());
        }
        else if (LCP[i] > LCP[get<0>(st.top())]) {
            st.push(make_tuple(i, i-1));
            R1[i] = get<1>(st.top());
        }
        else {
            while (LCP[i] < LCP[get<0>(st.top())]) {
                lastr1 = get<1>(st.top());
                st.pop();
            }
            if ((LCP[i] != 0) && (LCP[i] == LCP[get<0>(st.top())]) ) {
                R1[i] = get<1>(st.top());
            }
            else if (LCP[i] == 0) {
                st.push(make_tuple(0, -1));
                R1[i] = -1;
            }
            else if (LCP[i] > LCP[get<0>(st.top())]) {
                st.push(make_tuple(i, lastr1));
                R1[i] = lastr1;
            }
        }
    }
}

void compute_RM() {
    RM.reserve(input_size);
    for (int i = 0; i < input_size; ++i) {
        if (RSF[i] != 0) {
            RM[i] = R1[i] + RSF[i] - 1;
        }
    }
}

vector<int> exrun(
    const int &i,
    const int &j,
    const int &lcp_i
    ) {

    vector<int> r = vector<int>(3);    

    for (int p = 0; p <= ((j - i + 1) / 2); ++p) { // don't need this, if (p = (j-i+1)/2) then ...
        if (runs[p].size() > 0) {
            for (auto run : runs[p]) {
                r = { get<0>(run), get<1>(run), p };
                if (r[0] <= i && j <= r[1]) {
                    return r;
                }
            }
        }
    }

    cout << "run = {0,0,0}" << endl; 
    return {0, 0, 0}; 

}

int compute_frequency(
    vector<int> r, 
    int i_prime, 
    int j_prime
    ) {

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
    int &rm
    ) { 

    int olpi = 0;
    lcp_i = LCP[index];

    vector<int> SA_temp(input.size());



    // copy elements of SA[r1..rm] to SA*[r1..rm]
    for (int q = r1; q <= rm; ++q) {
        SA_temp[q] = SA[q];
    }

    // sort SA_temp in ascending order using mergesort O(nlogn)
    if (r1 < rm) {
        sort(SA_temp.begin() + r1, SA_temp.begin() + rm + 1);
    }

    int k = r1;
    while (k < rm) {
        // determines if there exists an overlap

        if ((SA_temp[k + 1] - SA_temp[k]) < lcp_i) {
            // if the difference between r[j] and r[j+1] is less than the length of the substring 

            vector<int> r = exrun(SA_temp[k], SA_temp[k + 1] + lcp_i - 1, lcp_i);
            // the frequency of the substring occuring in the run

            int fru = compute_frequency(r, SA_temp[k], SA_temp[k] + lcp_i - 1);

            // compute the overlap which is the freq of the substring minus one, multipled by the amount of overlap (which is length of the substring minus the period)
            olpi = olpi + (lcp_i - r[2]) * (fru-1);
            k = k + fru-1;
        }
        else {
            k++;
        }
    }

    return olpi;
}
