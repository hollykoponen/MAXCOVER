#ifndef OPTIMAL_COVERS_GLOBALS_H
#define OPTIMAL_COVERS_GLOBALS_H

// INCLUDES ====================================================================

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <set>
#include <stack>
#include <map>

using namespace std;

/* -------- Global Variables Declaration -------- */
extern string input;
extern int input_size;
extern int sorted_i; 
extern int sorted_j;  
extern vector<int> SA;
extern vector<int> LCP;
extern vector<int> RSF;
extern vector<int> RSF_all;
extern vector<int> R1;
extern vector<int> RM;
extern vector<int> ISA;
extern vector<int> ILCP;
extern vector<int> IRSF;
extern vector<pair<int, int>> NE;
extern string reversed_input;
extern vector<int> RANK;
extern vector<int> SA_rev;
extern vector<int> LCS;
extern vector<int> revRANK;
extern vector<set<pair<int, int>>> runs;
extern map<int, std::vector<int>> runsHT;
extern vector<int> OLP;
extern vector<int> RSPC;
extern vector<int> OCList;

#endif // OPTIMAL_COVERS_GLOBALS_H
