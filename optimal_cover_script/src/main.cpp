//
// Created by vi.
//
#include <execution>
#include <time.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <stack>
#include <string>
#include <iomanip>
#include <map>
#include "olp.h"
#include "sais/sais.h"
#include "lcpdc/dc_lcp_oracle.h"


/* compute SA array of string */
std::vector<int> compute_sa(const std::string &input) {
    int *sa = (int*)malloc(input.size() * sizeof(int));

    sais((unsigned char*)input.c_str(), sa, input.size());
    std::vector<int> sa_vec(sa, sa + input.size());

    free(sa);
    
    return sa_vec;
}

/* compute LCP array of string */
std::vector<int> compute_lcp(std::string &input, std::vector<int> &sa) {

    int *sa_arr = (int*)calloc(sa.size(), sizeof(int));
    for (int i = 0; i < sa.size(); ++i) {
        sa_arr[i] = sa[i];
    }

    //dc_lcp_initialize((unsigned char*)input.c_str(), lcp.data(), lcp.size(), 2);
    dc_lcp_initialize((unsigned char*)input.c_str(), sa_arr, sa.size(), 2);
    dc_lcp_construct();

    std::vector<int> lcp({ 0 });
    lcp.insert(lcp.end(), sa_arr, sa_arr + sa.size() - 1);

    dc_lcp_free();
    free(sa_arr);

    return lcp;
}

/* compute RSF array of string */
std::vector<int> compute_rsf(
        std::string &input,
        std::vector<int> &SA,
        std::vector<int> &LCP //,
        // std::vector<int> &R1,
        // std::vector<int> &RM
    ) {

    std::vector<int> rsf(SA.size());
    std::stack<int> st;
    int i = 0;

    while (i < SA.size()) {
        if (st.empty()) {
            st.push(i);
        } else {
            if (LCP[st.top()] < LCP[i]) {
                st.push(i);
            } else if (LCP[st.top()] == LCP[i]){
                rsf[i] = 0;
            } else {
                while (!st.empty() && LCP[st.top()] > LCP[i]) {
                    rsf[st.top()] = i - st.top();
                    // RM[st.top()] = i - 1;
                    st.pop();
                }
                if (!st.empty() && LCP[st.top()] == LCP[i] && LCP[i] != 0) {
                    rsf[i] = 0;
                } else {
                    st.push(i);
                }
            }
        }
        ++i;
    }

    while (!st.empty()) {
        if (LCP[st.top()] != 0) {
            rsf[st.top()] = SA.size() - st.top();
        } else {
            rsf[st.top()] = 0;
        }
        st.pop();
    }

    i = SA.size() - 1;
    while(i >= 0) {
        if(st.empty()) {
            st.push(i);
        } else {
            if(LCP[st.top()] <= LCP[i]) {
                st.push(i);
            } else {
                while(!st.empty() && LCP[st.top()] > LCP[i]){
                    if(rsf[st.top()] != 0) {
                        rsf[st.top()] = rsf[st.top()] + st.top() - i;
                    }
                    st.pop();
                }
                st.push(i);
            }
        }
        i = i - 1;
    }
    
    /*
    for (auto i : R1.length()){
        if (rsf[i] != 0){ 
            R1[i] = RM[i] - RSF[i] + 1;
        }
    }
    */
    
    return rsf;
}

/* compute RSF_all array of string */
std::vector<int> compute_rsf_all(
        std::string &input,
        std::vector<int> &SA,
        std::vector<int> &LCP
    ) {
    
    std::vector<int> RSF_all(SA.size());
    std::stack<int> st;
    int i = 0;

    while (i < SA.size()) {
        if (st.empty()) {
            st.push(i);
        }
        else {
            if (LCP[st.top()] <= LCP[i]) {
                st.push(i);}
            else {
                while (!st.empty() && LCP[st.top()] > LCP[i]) {
                    RSF_all[st.top()] = i - st.top();
                    st.pop();
                }
                    st.push(i);
            }
        }
        ++i;
    }

    while (!st.empty()) {
        if (LCP[st.top()] != 0) {
            RSF_all[st.top()] = SA.size() - st.top();
        } else {
            RSF_all[st.top()] = 0;
        }
        
        st.pop();
    }

    i = SA.size() - 1;
    while(i >= 0) {
        if(st.empty()) {
            st.push(i);
        } else {
            if(LCP[st.top()] <= LCP[i]) {
                st.push(i);
            } else {
                while(!st.empty() && LCP[st.top()] > LCP[i]){
                    if(RSF_all[st.top()] != 0) {
                        RSF_all[st.top()] = RSF_all[st.top()] + st.top() - i;
                    }
                    st.pop();
                }
                st.push(i);
            }
        }
        i = i - 1;
    }
       
    /*
    std::cout << "RSF_all: ";
    for (auto i : RSF_all){std::cout << i << ", "; }
    std::cout << std::endl;
    */
    
    return RSF_all;
}


/* compute OLP array of string */
/* QUADRATIC implementation */
std::vector<int> compute_olp(
    std::string &input,
    std::vector<int> &SA,
    std::vector<int> &LCP,
    std::vector<int> &RSF,
    std::vector<int> &R1,
    std::vector<int> &RM
    ) {

    std::vector<int> RANK = compute_rank(SA);

    std::string reverse(input);
    std::reverse(reverse.begin(), reverse.end());
    std::vector<int> SA_rev = compute_sa(reverse);

    std::vector<int> LCS = compute_lcp(reverse, SA_rev);
    // std::vector<int> LCS = {0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9}; 

    std::vector<int> revRANK = compute_rank(SA_rev);

    std::vector<std::set<std::pair<int, int>>> runs = compute_runs(LCP, LCS, RANK, revRANK);

    std::vector<int> OLP = std::vector<int>(input.size(), 0);

    for (int i = 0; i < OLP.size(); i++) {
    // std::cout << i << " out of " << OLP.size() << " " << OLP.size() << "left..." << std::endl;
        if (RSF[i] != 0) {
            OLP[i] = compute_olpi(i, R1[i], RM[i], LCP, SA, runs, input);
        }
    }

    return OLP;
}




/* compute OLP_nlogn array of string */
/* O(nlogn) implementation */
// TODO: Update header file to include OLP_nlogn functions
// TODO: Update makefile to ensure it compiles with new functions
std::vector<int> compute_OLP_nlogn(
    std::string &input,
    std::vector<int> &SA,
    std::vector<int> &LCP,
    std::vector<int> &RSF,
    std::vector<int> &R1,
    std::vector<int> &RM
    ) {
    
    int i = 2;
    std::vector<int> OLP_nlogn = std::vector<int>(input.size(), 0);
    OLP_nlogn[1] = 0;
    std::stack<pair<int,int>> stack; // stack of pairs (index, r1)
    std::pair<int,int> top; // (index, r1)
    // int runsHT.slots = 0; // # of slots filled in hashtable runsHT
    std::map<int, std::array<int, 3>> runsHT; // Hashtable int Key, tuple r = (i, j, p) Value 
    int Sorted_LI = 0; // lower index of sorted range in SA
    int Sorted_UI = 0;  // upper index of sorted range in SA
    
    while (i <= input.size()) {
        if (stack.empty()){
            stack.push( std::make_pair(i,i-1) );
        }
        else {
            top = stack.top();
            if (LCP[top.first] < LCP[i]) { stack.push( std::make_pair(i, i-1) ); }
            else if (LCP[top.first] = LCP[i]) { OLP_nlogn[i] = 0; }
            else{
                if (LCP[i] <= 1){
                    runsHT.clear();  // clears the contents of hash table; original: runsHT.slots = 0
                    Sorted_LI = 0; 
                    Sorted_UI = 0;
                }
                while (LCP[top.first] > LCP[i]) {
                    OLP_nlogn[top.first] = 0;
                    if(LCP[top.first] != 1) {
                        compute_Ru(top.first, top.second, i-1, Sorted_LI, Sorted_UI, SA, LCP, runsHT); 
                        OLP_nlogn[top.first] = compute_OLP_nlogn_at_index(top.first, LCP, runsHT);
                        compute_sorted_range(Sorted_LI, Sorted_UI, top.second, i-1); 
                    }
                    stack.pop();
                    top = stack.top();
                }
                if (top.first = LCP[i]) { OLP_nlogn[i] = 0; }
                else { stack.push( std::make_pair(i, top.first-1) ); }
            }
        }
        i++;
    }
    compute_OLP_nlogn_stack(i, Sorted_LI, Sorted_UI, top, SA, LCP, OLP_nlogn, stack, runsHT); 
    return OLP_nlogn;
}



/* compute RSPC array of string */
std::vector<int> compute_rspc(
    const std::string &input,
    std::vector<int> &LCP,
    std::vector<int> &RSF,
    std::vector<int> &OLP
) {

    std::vector<int> RSPC(input.size(), 0);

    for (int i = 0; i < input.size(); i++) {
        RSPC[i] = RSF[i] * LCP[i] - OLP[i];
    }

    return RSPC;
}


/* compute all optimal covers of a string */
std::vector<int> compute_optimal_covers(
    const std::vector<int> &LCP,
    std::vector<int> &RSPC
) {
    int max_pc = 1;
    int max_cover_length = 0;
    std::vector<int> OCList(RSPC.size()); // allocate memory to OCList array
    int list_length = 0;                                          // position to place elements in array list

    for (int i = 0; i < RSPC.size(); i++) {
        if (RSPC[i] > 1) {
            if (RSPC[i] > max_pc) {
                max_pc = RSPC[i];
                max_cover_length = LCP[i];
                list_length = 0;            //reset list_length position counter
                OCList[list_length] = i; // add i to the list
                list_length++;
            } else if (RSPC[i] == max_pc) {
                if (max_cover_length == LCP[i]) {
                    OCList[list_length] = i; // add i to the list
                    list_length++;
                } else if (max_cover_length < LCP[i]) {
                    max_cover_length = LCP[i];
                    list_length = 0;            //reset list_length position counter
                    OCList[list_length] = i; // add i to the list
                    list_length++;
                }
            }
        }
    }

    /* Set new OCListarr to the values computed above */
    std::vector<int> OCListarr(list_length, 0);
    for (int j = 0; j < list_length; j++) {
        OCListarr[j] = OCList[j];
    }

    return OCListarr;
}

std::vector<int> compute_top_ten_covers(
        const std::vector<int> &LCP,
        std::vector<int> &RSPC
) {
    std::vector<int> out;
//    int max_val = LCP[0];
    int min_vi = 0;
    for (int i = 0; i < LCP.size(); ++i) {
        if (out.size() < 10) {
            out.push_back(i);
            min_vi = (LCP[out[min_vi]] > LCP[i]) ? (out.size() - 1) : min_vi;
        } else if (LCP[i] > LCP[out[min_vi]]) {
            out[min_vi] = i;
            min_vi = std::min_element(out.begin(), out.end(), [&LCP](const int &a, const int &b) {
                return LCP[a] < LCP[b];
            }) - out.begin();
        }
    }
    return out;
}

/* compute frequency (%) of optimal cover */
int compute_oc_freq(int index, std::vector<int> &RSPC) {
    int frequency = (int)(((float)RSPC[index] / RSPC.size()) * 100);
    return frequency;
}

/* Driver function to sort the vector elements by second element of pairs in descending order */
bool sortbysec(const pair<int,int> &a,
              const pair<int,int> &b)
{
    return (a.second > b.second);
}

/* compute non-extendible repeating substrings in a string longer than repeat_size*/
std::vector<std::pair<int, int>> compute_ne(
        std::vector<int> &SA, 
        std::vector<int> &LCP, 
        std::vector<int> &RSF_all, 
        std::vector<int> &R1,
        std::vector<int> &RM, 
        int repeat_size)
    {
    std::vector<int> ISA = std::vector<int>(SA.size());
    std::vector<int> ILCP = std::vector<int>(LCP.size());
    std::vector<int> IRSF = std::vector<int>(RSF_all.size());
        
    for (int i = 0; i < SA.size(); i++){ ISA[SA[i]] = i;}
    for (int i = 0; i < LCP.size(); i++){ ILCP[i] = LCP[ISA[i]];}
    for (int i = 0; i < RSF_all.size(); i++){ IRSF[i] = RSF_all[ISA[i]];}
    
    /*
    std::cout << "ISA: ";
    for (auto i : ISA){std::cout << i << ", "; }
    std::cout << std::endl;

    std::cout << "ILCP: ";
    for (auto i : ILCP){std::cout << i << ", "; }
    std::cout << std::endl;
        
    std::cout << "IRSF: ";
    for (auto i : IRSF){std::cout << i << ", "; }
    std::cout << std::endl;
    */
    
    std::vector<std::pair<int, int>> NE;
    
    if ((IRSF[0] != 0) && (ILCP[0] >= repeat_size)){
        //Output NE repeating substring pair (1, ILCP[1])
        NE.push_back(std::make_pair(ISA[0], ILCP[0]));
    }
    int i = 1;
    while (i < SA.size()){
        if ((IRSF[i] != 0) && (ILCP[i] >= repeat_size) && (RM[ISA[i]] != 0)){
            if (not (IRSF[i]==IRSF[i-1] && ILCP[i]+1 == ILCP[i-1])){
                //Output NE repeating substring pair (i, i+ILCP[i]-1)
                NE.push_back(std::make_pair(ISA[i], ILCP[i]));
            }
        }
        i++;
    }
    
    sort(NE.begin(), NE.end(), sortbysec);
    
    /*
    std::cout << "NE: ";
    for (int i = 0; i < NE.size(); i++){
        std::cout << "(" << NE[i].first << "," << NE[i].second << ");";
    }
    std::cout << std::endl;
    */
    
    //std::vector<std::pair<int, int>> NE_Matches;
    
    /*
    i = 0;
    //set<int> NE_indices;
    while (i < NE.size()){
        //NE_indices.insert(NE[i].first);
        for(int r = R1[ISA[NE[i].first]]; r <= RM[ISA[NE[i].first]]; r++){ 
            //if(not (NE_indices.find(SA[r]) != NE_indices.end())){
                //NE_indices.insert(SA[r]);
                NE_Matches.push_back(std::make_pair(SA[r], NE[i].second));
            //}
        }
        i++;
    }
    */
        
    //sort(NE_Matches.begin(), NE_Matches.end(), sortbysec);
    
    /*std::cout << "NE: ";
    for (int i = 0; i < NE.size(); i++){
        std::cout << "(" << NE[i].first << "," << NE[i].second << ");";
    }
    std::cout << std::endl;
    */
    //return NE_Matches;
    return NE;
}



/* compute optimal covers, the main script */
void Compute_OC_main(std::string input, int top_ten){
    auto t_start = clock();
    std::vector<int> SAarr = compute_sa(input);
    std::vector<int> LCParr = compute_lcp(input, SAarr);
    std::vector<int> RSFarr = compute_rsf(input, SAarr, LCParr);
    std::vector<int> R1 = compute_r1(LCParr, RSFarr);
    std::vector<int> RM = compute_rm(R1, RSFarr);    
    std::vector<int> OLParr =  compute_olp(input, SAarr, LCParr, RSFarr, R1, RM);
    std::vector<int> RSPCarr = compute_rspc(input, LCParr, RSFarr, OLParr);
    std::vector<int> OCListarr = compute_optimal_covers(LCParr, RSPCarr);
    std::vector<int> top_ten_covers;
    if(top_ten == 1){
        top_ten_covers = compute_top_ten_covers(LCParr, RSPCarr);
    }
    auto end = clock();


    /* print all Optimal Covers and frequency */
    std::cout << "string of length " << input.size() << std::endl;
    std::cout << "\nOptimal Covers:" << std::endl;
    for (auto a : OCListarr) {
        int start = SAarr[a];
        int length = LCParr[a];
        std::cout << "\t" << input.substr(start, length) << std::endl; 
        std::cout << "\tOC length=" << length << std::endl;

        int OC_freq = compute_oc_freq(a, RSPCarr);
        std::cout << "\tnumber of positions covered=" << RSPCarr[a]
                  << "\n\t% covered = " << OC_freq << "%"
                  << "\n\ttime elapsed " << ((float)(end - t_start)) / CLOCKS_PER_SEC << " seconds\n" << std::endl;
    }
    
    if(top_ten == 1){
        std::cout << "Top Ten Covers: " << std::endl;
        for (auto a : top_ten_covers){
            std::cout << "out[min_vi] = " << a << " Position: " << SAarr[a] << " of length: " << LCParr[a] << std::endl;
        }
    }
    
}


/* compute repeat matches, the main script, using repeat_size (int) to determine minimum repeat length */
void Compute_RepeatMatches(std::string input, int repeat_size){
    auto t_start = clock();
    std::vector<int> SAarr = compute_sa(input);
    std::vector<int> LCParr = compute_lcp(input, SAarr);
    std::vector<int> RSFarr = compute_rsf(input, SAarr, LCParr);
    std::vector<int> RSF_all = compute_rsf_all(input, SAarr, LCParr);
    std::vector<int> R1 = compute_r1(LCParr, RSFarr);
    std::vector<int> RM = compute_rm(R1, RSFarr);
    std::vector<std::pair<int, int>> NE = compute_ne(SAarr, LCParr, RSF_all, R1, RM, repeat_size);
    auto end = clock();
    
    /*
    std::cout << "SA: ";
    for (auto i : SAarr){std::cout << i << ", "; }
    std::cout << std::endl;

    std::cout << "LCP: ";
    for (auto i : LCParr){std::cout << i << ", "; }
    std::cout << std::endl;
    
    std::cout << "RSF: ";
    for (auto i : RSFarr){std::cout << i << ", "; }
    std::cout << std::endl;
    
    std::cout << "R1: ";
    for (auto i : R1){std::cout << i << ", "; }
    std::cout << std::endl;
    
    std::cout << "RM: ";
    for (auto i : RM){std::cout << i << ", "; }
    std::cout << std::endl;
    */
    
    /* print the repeat matches */
    std::cout << "Long Exact Matches:" << std::endl;
    std::cout << std::right << std::setw(10) << "Start1" << std::setw(10) << "Start2" << std::setw(10) << "Length" << std::endl;
    
    /*
    for (int i = 0; i < input.length(); i++){
        if ((RM[i] != 0) && (LCParr[i] >= repeat_size)){
            for(int r = R1[i]+1; r <= RM[i]; r++){
                int Start1 = SAarr[R1[i]]+1; // SAarr[R1]
                int Start2 = SAarr[r]+1; // SAarr[RM] for (all R1+1 to RM) loop
                int Length = LCParr[i]; 
                std::cout << std::right << std::setw(10) << Start1 << std::setw(10) << Start2 << std::setw(10) << Length << std::endl;            
            }
        }
    }
    */

    std::vector<int> Reported;
    int start1, start2;
    
    for (int i = 0; i < NE.size(); i++){
        //std::cout << "NE[" << i << "]: " << NE[i].first << "," << NE[i].second << std::endl;
        for(int r1 = R1[NE[i].first]; r1 <= RM[NE[i].first]; r1++){
            start1 = SAarr[r1];
            //std::cout << "Start1: " << start1 << std::endl;
            if (std::find(Reported.begin(), Reported.end(), start1) == Reported.end()){
                for(int r2 = r1; r2 <= RM[NE[i].first]; r2++){
                    start2 = SAarr[r2];
                    if(start1 != start2){
                        std::cout << std::right;
                        std::cout << std::setw(10) << start1+1;
                        std::cout << std::setw(10) << start2+1; 
                        std::cout << std::setw(10) << NE[i].second;
                        std::cout << std::endl;
                    }
                }
            }
            Reported.push_back(start1);
        }
    }
    

  /*    
    int Start1 = NE[0].first;
    int Start2;
    int Length = NE[0].second;

    int i = 1;    
    while (i < NE.size()){
        if(NE[i].second == Length){
            Start2 = NE[i].first;            
            std::cout << std::right << std::setw(10) << Start1+1 << std::setw(10) << Start2+1 << std::setw(10) << Length << std::endl;
        }
        else{
            Start1 = NE[i].first;
            Length = NE[i].second;
        }
        i++;        
    }
  */  
    
    // std::cout << "Repeating Substring" << i << "of Length" << NE[i].second << "=" << for loop {11, 1, 7} << std::endl;
    
}

int main(int argc, char* argv[]) {
    
    /*std::string input = "aababcaba";
    Compute_OC_main(input);
    Compute_RepeatMatches(input, 1);

    input = "abacababacabacaba";
    Compute_OC_main(input);
    Compute_RepeatMatches(input, 1);
    */

    int Compute_OC = std::atoi(argv[1]); // 1 = True; 0 = False;
    int top_ten = std::atoi(argv[2]); // 1 = True; 0 = False;
    int RepeatMatches = std::atoi(argv[3]); // 1 = True; 0 = False;
    int repeat_size = std::atoi(argv[4]); // int, 1+
    
    /*
    std::cout << "Compute_OC: " << Compute_OC << std::endl;
    std::cout << "RepeatMatches: " << RepeatMatches << std::endl;
    std::cout << "repeat_size: " << repeat_size << std::endl;
    */
        

        
    /* initialize input string*/
    std::ifstream ifs(argv[5]);
    std::string input_EOL(
        (std::istreambuf_iterator<char>(ifs)),
        (std::istreambuf_iterator<char>())
    );
    std::string input = input_EOL.substr(0, input_EOL.size()-1);
    //std::cout << "Input String Length: " << input.size() << std::endl;
    
    if (Compute_OC == 1) {Compute_OC_main(input, top_ten);}
    if (RepeatMatches == 1){Compute_RepeatMatches(input, repeat_size);}
    
    //std::cout << "Input String: \n\t" << input << std::endl;
    
    return 0;
    
}
