#include "readData.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <chrono>
#include <limits>
#include <cmath>

using namespace std;

enum NL{N1, N2, N3, N4, N5};

double **costM;
int dimension;

typedef struct{
    double C;
    double T;
    int W;
} ReoptData;

void printCostM() 
{
    std::cout << "dimension: " << dimension << endl;
    for (size_t i = 1; i <= dimension; i++) {
        for (size_t j = 1; j <= dimension; j++) {
            std::cout << costM[i][j] << " ";
        }
        std::cout << endl;
    }
    cout << endl << endl;
}

vector<int> construction()
{
    vector<int> solution = {1};
    vector<int> candidateList;

    //Fills candidates list according to dimension
    for(int i = 2; i <= dimension; i++){
        candidateList.push_back(i);
    }

    int r = 0;
    while(!candidateList.empty())
    {
        int closestNode;
        double closestNodeDistance = std::numeric_limits<double>::max();

        for(int k : candidateList)
        {
            if(costM[solution[r]][k] < closestNodeDistance)
            {
                closestNodeDistance = costM[solution[r]][k];
                closestNode = k;
            }    
        }

        solution.push_back(closestNode);
        candidateList.erase(std::remove(candidateList.begin(), candidateList.end(), closestNode), candidateList.end());
        r++;
    }

    solution.push_back(1);

    return solution;
}

double solutionCost (vector <int> s)
{
    ReoptData cost;
    const int i = 0, j = s.size();

    cost.C = 0;
    cost.T = 0;

    for(int k = i; k < j; k++)
    {
        cost.T += costM[s[k]][s[k+1]];
        cost.C += costM[s[k]][s[k+1]] + cost.T;
    }

    return cost.C;
}

double solutionCost2(vector<int> s) //DEBUG
{
    vector<vector<ReoptData>> reopt = calcReopt(s);
    return reopt[0][0].C;
}

void calcReopt(vector<int> s, vector<vector<ReoptData>> *reopt, int ii, int jj)
{
    int begin = (ii <= jj)? ii : jj;
    
    for(int j = begin-1; j < s.size(); j++)
    {
        for(int i = begin-1; i <= j; i++)
        {
            (*reopt)[s[i]][s[j]].C = 0;
            (*reopt)[s[i]][s[j]].T = 0;
            (*reopt)[s[i]][s[j]].W = (s[i] == 1) ? j-i : j-i+1;
            //(*reopt)[s[i]][s[j]].W = (i == j) ? 1 : j-i;

            for(int k = i; k < j; k++)
            {
                (*reopt)[s[i]][s[j]].T += costM[s[k]][s[k+1]];
                (*reopt)[s[i]][s[j]].C += costM[s[k]][s[k+1]] + (*reopt)[s[i]][s[k]].T;
            }
        }
    }

    //Reverse (for 2opt)
    for(int i = begin; i < s.size()-1; i++)
    {
        for(int j = begin; j < i; j++)
        {
            (*reopt)[s[i]][s[j]].C = 0;
            (*reopt)[s[i]][s[j]].T = 0;
            (*reopt)[s[i]][s[j]].W = (s[i] == 1) ? i-j : i-j+1;

            for(int k = i, l = j; l < i; k--, l++)
            {
                (*reopt)[s[i]][s[j]].T += costM[s[k]][s[k-1]];
            }
        }
        for(int j = begin; j < i; j++)
        {
            for(int k = i, l = j; l < i; k--, l++){
                (*reopt)[s[i]][s[j]].C += costM[s[k]][s[k-1]] + (*reopt)[s[i]][s[j]].T;
            }
        }
    }

    (*reopt)[0][0].W = (*reopt)[s[0]][s[s.size()-3]].W + (*reopt)[s[s.size()-2]][s[s.size()-1]].W;
    (*reopt)[0][0].T = (*reopt)[s[0]][s[s.size()-3]].T + (*reopt)[s[s.size()-2]][s[s.size()-1]].T;
    (*reopt)[0][0].C = (*reopt)[s[0]][s[s.size()-3]].C + (*reopt)[s[s.size()-2]][s[s.size()-1]].W*((*reopt)[s[0]][s[s.size()-3]].T + (*reopt)[s[s.size()-2]][s[s.size()-1]].T) + (*reopt)[s[s.size()-2]][s[s.size()-1]].C;
}

vector<vector<ReoptData>> calcReopt(vector<int> s)
{
    vector<vector<ReoptData>> reopt(s.size(), vector<ReoptData>(s.size()));
    calcReopt(s, &reopt, 1, 1);
    return reopt;
}

ReoptData concatCost(vector<vector<ReoptData>> reopt, int u, int v, int w, int x)
{ // O = (Ou, ... , Ov) and O' = (O'w, ... , O'x) are our subsequences
    ReoptData cost;

    cost.T = reopt[u][v].T + reopt[v][w].T + reopt[w][x].T;
    cost.C = reopt[u][v].C + reopt[w][x].W * (reopt[u][v].T + reopt[v][w].T) + reopt[w][x].C;
    cost.W = reopt[u][v].W + reopt[w][x].W;

    return cost;  
}

ReoptData concatCost(vector<vector<ReoptData>> reopt, ReoptData sq1, int v, int w, int x)
{ // O = (Ou, ... , Ov) = sq1 and O' = (O'w, ... , O'x) are our subsequences
    ReoptData cost;
    int reopt_w_x_W = reopt[w][x].W;
    if(w == 1 && x == 1)
        reopt_w_x_W = 1; // Force reopt[1][1] = 1, as here 1 always comes from s[last]

    cost.T = sq1.T + reopt[v][w].T + reopt[w][x].T;
    cost.C = sq1.C + reopt_w_x_W * (sq1.T + reopt[v][w].T) + reopt[w][x].C;
    cost.W = sq1.W + reopt_w_x_W;

    return cost;  
}

vector<int> swap (vector<int> s, vector<vector<ReoptData>> *reopt, bool *improved){
    int    best_i = 0, best_j = 0;
    double bestDelta = 0, delta = 0;
    double preMoveCost = (*reopt)[0][0].C;
    ReoptData sq[3];
    ReoptData cost;

    for(int j = 2; j < s.size() - 1; j++)
    {
        for(int i = 1; i < j; i++)
        {
            cost.C = 0;
            cost.T = 0;
            cost.W = 0;

            if(i == j - 1){ // If they're adjacent
                sq[0] = concatCost(*reopt, s[0], s[i-1], s[j], s[j]);
                sq[1] = concatCost(*reopt, sq[0], s[j], s[i], s[i]);
                cost  = concatCost(*reopt, sq[1], s[i], s[j+1], s[s.size()-1]);
                //reopt[0][i-1] o reopt[j][j] o reopt[i][i] o reopt[j+1][dimension]
            }else{
                sq[0] = concatCost(*reopt, s[0], s[i-1], s[j], s[j]);
                sq[1] = concatCost(*reopt, sq[0], s[j], s[i+1], s[j-1]);
                sq[2] = concatCost(*reopt, sq[1], s[j-1], s[i], s[i]);
                cost  = concatCost(*reopt, sq[2], s[i], s[j+1], s[s.size()-1]);
                //reopt[0][i-1] o reopt[j][j] o reopt[i+1][j-1] o reopt[i][i] o reopt[j+1][dimension]
            }

            if(i == 47 && j == 51) //DEBUG
            {
                int a = i + j;
            }

            delta = cost.C - preMoveCost;
            if(delta < bestDelta){
                bestDelta = delta;
                *improved = true;

                best_i = i;
                best_j = j;
            }
        }
    }

    if(*improved)
    {
        std::swap(s[best_j], s[best_i]);
        calcReopt(s, reopt, best_i, best_j);
    }

    return s;
}

vector<int> flip (vector<int> s, vector<vector<ReoptData>> *reopt, bool *improved){
    int best_i = 0, best_j = 0;
    double bestDelta = 0, delta = 0;
    double preMoveCost = (*reopt)[0][0].C;
    ReoptData sq;
    ReoptData cost;

    for(int j = 3; j < s.size() - 1; j++)
    {
        for(int i = 1; i < j - 1; i++)
        {
            sq = concatCost(*reopt, s[0], s[i-1], s[j], s[i]);
            cost = concatCost(*reopt, sq, s[i], s[j+1], s[s.size()-1]);
            //reopt[0][i-1] o reopt[j][i] o reopt[j+1][dimension]
            
            delta = cost.C - preMoveCost;
            if(delta < bestDelta){
                bestDelta = delta;
                *improved = true;

                best_i = i;
                best_j = j;
            }
        }
    }

    if(*improved)
    {
        std::reverse(s.begin() + best_i, s.begin() + best_j + 1);
        calcReopt(s, reopt, best_i, best_j);
    }

    return s;
}

vector<int> reinsertion (vector<int> s, vector<vector<ReoptData>> *reopt, bool *improved, int subsegSize){
    int best_i = 0, best_j = 0;
    double bestDelta = 0, delta = 0;
    double preMoveCost = (*reopt)[0][0].C;
    ReoptData sq[2];
    ReoptData cost;

    for(int i = 1; i < s.size() - subsegSize; i++)
    {
        for(int j = 1; j < s.size(); j++)
        {
            if(i <= j && j <= i + subsegSize) continue;

            if(i < j){
                sq[0] = concatCost(*reopt, s[0], s[i-1], s[i+subsegSize], s[j-1]);
                sq[1] = concatCost(*reopt, sq[0], s[j-1], s[i], s[i+subsegSize-1]);
                cost = concatCost(*reopt, sq[1], s[i+subsegSize-1], s[j], s[s.size()-1]);

                //reopt[0][i-1] o reopt[i+subsegSize][j-1] o reopt[i][i+subsegSize-1] o reopt[j][dimension]
            }else{
                sq[0] = concatCost(*reopt, s[0], s[j-1], s[i], s[i+subsegSize-1]);  
                sq[1] = concatCost(*reopt, sq[0], s[i+subsegSize-1], s[j], s[i-1]);
                cost = concatCost(*reopt, sq[1], s[i-1], s[i+subsegSize], s[s.size()-1]);

                //reopt[0][j-1] o reopt[i][i+subsegSize-1] o reopt[j][i-1] o reopt[i+subsegSize][dimension]
            }
            
            delta = cost.C - preMoveCost;
            if(delta < bestDelta){
                bestDelta = delta;
                *improved = true;

                best_i = i;
                best_j = j;
            }
        }
    }

    if(*improved)
    {
        vector<int> subseg(s.begin() + best_i, s.begin() + best_i + subsegSize);

        if(best_i < best_j){
            s.insert(s.begin() + best_j, subseg.begin(), subseg.end());
            s.erase(s.begin() + best_i, s.begin() + best_i + subsegSize);
        }else{
            s.erase(s.begin() + best_i, s.begin() + best_i + subsegSize);
            s.insert(s.begin() + best_j, subseg.begin(), subseg.end());
        }

        calcReopt(s, reopt, best_i, best_j);
    }

    return s;
}

vector<int> RVND (vector<int> s, double *mainCost){
    vector<int> ngbhList = {N1/*, N2, N3, N4, N5*/};
    int ngbh_n;

    vector<int> neighbour_s = s;
    bool improved;

    vector<vector<ReoptData>> costs = calcReopt(s);

    int i = 0; //DEBUG

    while(!ngbhList.empty())
    {
        ngbh_n = ngbhList[rand() % ngbhList.size()];

        switch(ngbh_n){
            case N1:
                neighbour_s = swap(s, &costs, &improved);
                break;
            case N2:
                neighbour_s = flip(s, &costs, &improved);
                break;
            case N3:
                neighbour_s = reinsertion(s, &costs, &improved, 1);
                break;
            case N4:
                neighbour_s = reinsertion(s, &costs, &improved, 2);
                break;
            case N5:
                neighbour_s = reinsertion(s, &costs, &improved, 3);
                break;
        }

        if(improved){
            s = neighbour_s;
            *mainCost = costs[0][0].C;
            improved = false;

            cout << "loop #" << i << '\t' << "IMPROVED " << (*mainCost) << endl; //DEBUG
            cout << "movement: " << ngbh_n << endl; // DEBUG

            ngbhList = {N1/*, N2, N3, N4, N5*/};
            //Reopt update is done in movement functions
        }else{
            cout << "loop #" << i << '\t' << "!improved" << endl; //DEBUG
            cout << "movement: " << ngbh_n << endl; // DEBUG

            ngbhList.erase(std::remove(ngbhList.begin(), ngbhList.end(), ngbh_n), ngbhList.end());
        }

        i++; //DEBUG
    }

    return s;
}

int randomRange(int min, int max)
{   // Inclusive
    return min + (rand() % (max - min + 1));
}

vector<int> perturb (vector<int> s){
    int subseg1Start = 1, subseg1End = 1;
    int subseg2Start = 1, subseg2End = 1;
    // If s.size()/10 >= 2 -> max = s.size()/10, else -> max = 2
    const int maxSubsegSize = (std::ceil(s.size()/10) >= 2) ? std::ceil(s.size()/10) : 2; 
    const int minSubsegSize = 2;

    while( subseg1Start <= subseg2Start && subseg2Start <= subseg1End
        || subseg2Start <= subseg1Start && subseg1Start <= subseg2End )
    {
        subseg1Start = randomRange(1, s.size() - 1 - maxSubsegSize);
        subseg1End = subseg1Start + randomRange(minSubsegSize, maxSubsegSize);

        subseg2Start = randomRange(1, s.size() - 1 - maxSubsegSize);
        subseg2End = subseg2Start + randomRange(minSubsegSize, maxSubsegSize);
    }

    const int subseg1Size = subseg1End - subseg1Start;
    const int subseg2Size = subseg2End - subseg2Start;

    vector<int> subseg1(s.begin() + subseg1Start, s.begin() + subseg1End);
    vector<int> subseg2(s.begin() + subseg2Start, s.begin() + subseg2End);

    std::reverse(subseg1.begin(), subseg1.end());
    std::reverse(subseg2.begin(), subseg2.end());

    if(subseg1End < subseg2Start){
        s.insert(s.begin() + subseg2Start, subseg1.begin(), subseg1.end());
        s.insert(s.begin() + subseg1Start, subseg2.begin(), subseg2.end());
        s.erase(s.begin() + subseg2Start + subseg1Size + subseg2Size, s.begin() + subseg2End + subseg1Size + subseg2Size);
        s.erase(s.begin() + subseg1Start + subseg2Size, s.begin() + subseg1End + subseg2Size);
    }else{
        s.insert(s.begin() + subseg1Start, subseg2.begin(), subseg2.end());
        s.insert(s.begin() + subseg2Start, subseg1.begin(), subseg1.end());
        s.erase(s.begin() + subseg1Start + subseg1Size + subseg2Size, s.begin() + subseg1End + subseg1Size + subseg2Size);
        s.erase(s.begin() + subseg2Start + subseg1Size, s.begin() + subseg2End + subseg1Size);
    }

    return s;
}

int main(int argc, char** argv) {
    srand(time(NULL));

    readData(argc, argv, &dimension, &costM);
    //std::cout << "\tCOST MATRIX: \n";
    //printCostM();

    auto timerStart = chrono::system_clock::now();

    const int I_MAX = 50;
    const int I_ILS = (dimension >= 150) ? (dimension/2) : (dimension);

    vector<int> solutionAlpha, solutionBeta, solutionOmega;
    double costAlpha, costBeta, costOmega = numeric_limits<double>::infinity();

    double costAlpha2, firstCost; //DEBUG

    for(int i = 0; i < I_MAX; i++)
    {
        solutionAlpha = construction();
        solutionBeta = solutionAlpha;

        costAlpha2 = solutionCost(solutionAlpha);
        costAlpha = solutionCost2(solutionAlpha); //DEBUG
        
        if(i == 0){ //DEBUG
            firstCost = costAlpha;
        }

        costBeta = costAlpha;

        for(int iterILS = 0; iterILS < I_ILS; iterILS++){

            cout << "Iter " << i << " of " << I_MAX << endl; //DEBUG
            cout << "IterILS " << iterILS << " of " << I_ILS << endl; //DEBUG

            solutionAlpha = RVND(solutionAlpha, &costAlpha);

            if(costAlpha < costBeta){
                solutionBeta = solutionAlpha;
                costBeta = costAlpha;
                iterILS = 0;
            }

            solutionAlpha = perturb(solutionBeta);
        }

        if(costBeta < costOmega){
            solutionOmega = solutionBeta;
            costOmega = costBeta;
        }
    }

    auto timerEnd = chrono::system_clock::now();
    chrono::duration<double> elapsedSeconds = timerEnd - timerStart;

    //PRINT COST AND SOLUTION
    std::cout << "SOLUTION:\n";
    for(auto k : solutionOmega){
        std::cout << k << ' ';
    }
    std::cout << "\n\n" << "COST: " << costOmega << "\n";
    std::cout << "TIME: " << elapsedSeconds.count() << "\n";
    if(costOmega == firstCost){
        cout << "Solution never changed" << endl;
    }

    return 0;
}