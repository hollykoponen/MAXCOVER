// INCLUDES =======================================================================================

#include "globals.h"
#include <execution>
#include <time.h>
#include <fstream>
#include <iomanip>
#include "olp.h"
#include "sais/sais.h"
#include "lcpdc/dc_lcp_oracle.h"

using namespace std;

/* -------- Global Variables Declaration -------- */
string input = "";
int input_size = 0;
vector<int> SA(0);
vector<int> LCP(0);
vector<int> RSF(0);
vector<int> RSF_all(0);
vector<int> R1(0);
vector<int> RM(0);
vector<int> ISA(0);
vector<int> ILCP(0);
vector<int> IRSF(0);
vector<pair<int, int>> NE(0);

string reversed_input = "";
vector<int> RANK(0);
vector<int> SA_rev(0);
vector<int> LCS(0);
vector<int> revRANK(0);

vector<int> OLP(0);
vector<int> RSPC(0);
vector<int> OCList(0);

// UTILS ==========================================================================================

void PrintStack(stack<pair<int,int>> s)
{
    if (s.empty())
        return;
    pair<int,int> x = s.top();
    s.pop();
    PrintStack(s);
    cout << "(" << x.first << "," << x.second << ") " << endl;;
    s.push(x);
}

/* Driver function to sort the vector elements by second element of pairs in descending order */
bool sortbysec(const pair<int,int> &a, const pair<int,int> &b) {return (a.second > b.second);}

/* compute frequency (%) of optimal cover */
int compute_oc_freq(int index) {
    int frequency = (int)(((float)RSPC[index] / input_size) * 100);
    return frequency;
}

#define watch(x) cout << (#x) << " = "

#define printvar(x) cout << (#x) << " = " << x << endl;

void PrintArray(vector<int> arr){
    for(auto a : arr){
        cout << a << ",";
    } cout << endl; cout << endl;
}

/* -------- Compute Array Functions -------- */

/* compute SA array of string */
vector<int> compute_sa(string input) {
    int *sa = (int*)malloc(input_size * sizeof(int));

    sais((unsigned char*)input.c_str(), sa, input_size);
    vector<int> sa_vec(sa, sa + input_size);

    free(sa);
    
    return sa_vec;
}

/* compute LCP array of string */
vector<int> compute_lcp(string input, vector<int> sa) {

    int *sa_arr = (int*)calloc(input_size, sizeof(int));
    for (int i = 0; i < input_size; ++i) {
        sa_arr[i] = sa[i];
    }

    dc_lcp_initialize((unsigned char*)input.c_str(), sa_arr, input_size, 2);
    dc_lcp_construct();

    vector<int> lcp({ 0 });
    lcp.insert(lcp.end(), sa_arr, sa_arr + input_size - 1);

    dc_lcp_free();
    free(sa_arr);

    return lcp;
}

void runs_for_exrun(){
    reversed_input = input;
    reverse(reversed_input.begin(), reversed_input.end());
    RANK = compute_rank(SA);
    SA_rev = compute_sa(reversed_input);
    LCS = compute_lcp(reversed_input, SA_rev);
    revRANK = compute_rank(SA_rev);
    runs = compute_runs(LCP, LCS, RANK, revRANK);
}

/* compute RSF array of string */
vector<int> compute_rsf() {

    vector<int> rsf(input_size);
    stack<int> st;
    int i = 0;

    while (i < input_size) {
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
            rsf[st.top()] = input_size - st.top();
        } else {
            rsf[st.top()] = 0;
        }
        st.pop();
    }

    i = input_size - 1;
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
        i--;
    }
    
    return rsf;
}

/* compute RSF_all array of string */
vector<int> compute_rsf_all() {
    
    vector<int> RSF_all(input_size);
    stack<int> st;

    int i = 0;
    while (i < input_size) {
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
        i++;
    }

    while (!st.empty()) {
        if (LCP[st.top()] != 0) {
            RSF_all[st.top()] = input_size - st.top();
        } else {
            RSF_all[st.top()] = 0;
        }
        
        st.pop();
    }

    i = input_size - 1;
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
        i--;
    }

    return RSF_all;
}

/* compute RSPC array of string */
vector<int> compute_rspc() {

    vector<int> rspc(input_size);

    for (int i = 0; i < input_size; i++) {
        rspc[i] = RSF[i] * LCP[i] - OLP[i];
    }

    return rspc;
}


/* compute all optimal covers of a string */
vector<int> compute_optimal_covers() {
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

/* -------- Flags Functions -------- */

vector<int> compute_top_ten_covers() {
    vector<int> out;
    int min_vi = 0;

    for (int i = 0; i < input_size; ++i) {
        if (out.size() < 10) {
            out.push_back(i);
            min_vi = (LCP[out[min_vi]] > LCP[i]) ? (out.size() - 1) : min_vi;
        } else if (LCP[i] > LCP[out[min_vi]]) {
            out[min_vi] = i;
            min_vi = min_element(out.begin(), out.end(), [](const int &a, const int &b) {
                return LCP[a] < LCP[b];
            }) - out.begin();
        }
    }
    return out;
}

/* compute non-extendible repeating substrings in a string longer than repeat_size*/
vector<pair<int, int>> compute_ne(int repeat_size) {
    
    ISA.reserve(input_size);
    ILCP.reserve(input_size);
    IRSF.reserve(input_size);
        
    for (int i = 0; i < input_size; i++){ 
        ISA[SA[i]] = i;
        ILCP[i] = LCP[ISA[i]];
        IRSF[i] = RSF_all[ISA[i]];
    }
        
    if ((IRSF[0] != 0) && (ILCP[0] >= repeat_size)){
        //Output NE repeating substring pair (1, ILCP[1])
        NE.push_back(make_pair(ISA[0], ILCP[0]));
    }

    int i = 1;
    while (i < input_size){
        if ((IRSF[i] != 0) && (ILCP[i] >= repeat_size) && (RM[ISA[i]] != 0)){
            if (not (IRSF[i]==IRSF[i-1] && ILCP[i]+1 == ILCP[i-1])){
                //Output NE repeating substring pair (i, i+ILCP[i]-1)
                NE.push_back(make_pair(ISA[i], ILCP[i]));
            }
        }
        i++;
    }
    
    sort(NE.begin(), NE.end(), sortbysec);
    
    return NE;
}


/* -------- Procedures -------- */

/* compute optimal covers, the main script */
void Compute_OC_main(int top_ten){
    
    auto t_start = clock();
    
    SA = compute_sa(input);
    LCP = compute_lcp(input, SA);
    RSF = compute_rsf();
    compute_R1();
    compute_RM();    

    PrintArray(SA); PrintArray(LCP); PrintArray(RSF);


    //OLP =  compute_olp(runs);
    OLP =  compute_OLP_nlogn(runs);
    RSPC = compute_rspc();
    OCList = compute_optimal_covers();

    vector<int> top_ten_covers;
    
    if(top_ten == 1){
        top_ten_covers = compute_top_ten_covers();
    }
    
    auto end = clock();


    /* print all Optimal Covers and frequency */
    cout << "string of length " << input_size << endl;
    cout << "\nOptimal Covers:" << endl;
    for (auto a : OCList) {
        int start = SA[a];
        int length = LCP[a];
        cout << "\t" << input.substr(start, length) << endl; 
        cout << "\tOC length=" << length << endl;

        int OC_freq = compute_oc_freq(a);
        cout << "\tnumber of positions covered=" << RSPC[a]
                  << "\n\t% covered = " << OC_freq << "%"
                  << "\n\ttime elapsed " << ((float)(end - t_start)) / CLOCKS_PER_SEC << " seconds\n" << endl;
    }
    
    if(top_ten == 1){
        cout << "Top Ten Covers: " << endl;
        for (auto a : top_ten_covers){
            cout << "out[min_vi] = " << a << " Position: " << SA[a] << " of length: " << LCP[a] << endl;
        }
    }
    
}

/* compute repeat matches, the main script, using repeat_size (int) to determine minimum repeat length */
void Compute_RepeatMatches(int repeat_size){
    auto t_start = clock();
    SA = compute_sa(input);
    LCP = compute_lcp(input, SA);
    RSF = compute_rsf();
    RSF_all = compute_rsf_all();
    compute_R1();
    compute_RM();
    NE = compute_ne(repeat_size);
    auto end = clock();
    
    /* print the repeat matches */
    cout << "Long Exact Matches:" << endl;
    cout << right << setw(10) << "Start1" << setw(10) << "Start2" << setw(10) << "Length" << endl;
    
    vector<int> Reported;
    int start1, start2;
    
    for (int i = 0; i < NE.size(); i++){
        //cout << "NE[" << i << "]: " << NE[i].first << "," << NE[i].second << endl;
        for(int r1 = R1[NE[i].first]; r1 <= RM[NE[i].first]; r1++){
            start1 = SA[r1];
            //cout << "Start1: " << start1 << endl;
            if (find(Reported.begin(), Reported.end(), start1) == Reported.end()){
                for(int r2 = r1; r2 <= RM[NE[i].first]; r2++){
                    start2 = SA[r2];
                    if(start1 != start2){
                        cout << right;
                        cout << setw(10) << start1+1;
                        cout << setw(10) << start2+1; 
                        cout << setw(10) << NE[i].second;
                        cout << endl;
                    }
                }
            }
            Reported.push_back(start1);
        }
    }
    
}

int main(int argc, char* argv[]) {

    int Compute_OC = atoi(argv[1]); // 1 = True; 0 = False;
    int top_ten = atoi(argv[2]); // 1 = True; 0 = False;
    int RepeatMatches = atoi(argv[3]); // 1 = True; 0 = False;
    int repeat_size = atoi(argv[4]); // int, 1+
    
        
    /* initialize input string*/
    ifstream ifs(argv[5]);
    string input_EOL(
        (istreambuf_iterator<char>(ifs)),
        (istreambuf_iterator<char>())
    );
    
    input = input_EOL.substr(0, input_EOL.size()-1);
    input_size = input.size();
    
    if (Compute_OC == 1) {Compute_OC_main(top_ten);}
    if (RepeatMatches == 1){Compute_RepeatMatches(repeat_size);}
    
    return 0;
    
}
