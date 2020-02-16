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

    for(int k = i; k < j; k++){
        cost.T += costM[s[k]][s[k+1]];
        cost.C += costM[s[k]][s[k+1]] + cost.T; //Checar se esse calculo ta certo
    }

    return cost.C;
}

void calcReopt(vector<int> s, vector<vector<ReoptData>> *reopt, int ii, int jj)
{
    int small = (ii <= jj)? ii : jj;
    
    for(int j = small; j < s.size(); j++)
    {
        for(int i = small; i <= j; i++)
        {
            (*reopt)[s[i]][s[j]].C = 0;
            (*reopt)[s[i]][s[j]].T = 0;
            (*reopt)[s[i]][s[j]].W = j - i + 1;

            if(i == j)
                continue;
            

            for(int k = i; k < j; k++)
            {
                (*reopt)[s[i]][s[j]].T += costM[s[k]][s[k+1]];
                (*reopt)[s[i]][s[j]].C += costM[s[k]][s[k+1]] + (*reopt)[s[i]][s[k]].T;
            }
        }
    }

    for(int i = small; i < s.size(); i++)
    {
        for(int j = small; j < i; j++)
        {
            (*reopt)[s[i]][s[j]].C = 0;
            (*reopt)[s[i]][s[j]].T = 0;
            (*reopt)[s[i]][s[j]].W = i - j + 1;

            for(int k = i, l = j; l < i; k--, l++)
            {
                (*reopt)[s[i]][s[j]].T += costM[s[k]][s[k-1]];
            }
        }
        for(int j = small; j < i; j++)
        {
            for(int k = i, l = j; l < i; k--, l++){
                (*reopt)[s[i]][s[j]].C += costM[s[k]][s[k-1]] + (*reopt)[s[i]][s[j]].T;
            }
        }
    }

    (*reopt)[s[0]][s[s.size()-1]].W = s.size();
    (*reopt)[s[0]][s[s.size()-1]].T = (*reopt)[s[0]][s[s.size()-2]].T + (*reopt)[s[s.size()-2]][s[s.size()-1]].T;
    (*reopt)[s[0]][s[s.size()-1]].C = (*reopt)[0][s.size()-2].C + (*reopt)[s[0]][s.size()-2].T + (*reopt)[s.size()-2][s.size()-1].T;
}

vector<vector<ReoptData>> calcReopt(vector<int> s){
    vector<vector<ReoptData>> reopt(s.size(), vector<ReoptData>(s.size()));
    
    for(int j = 1; j < s.size(); j++)
    {
        for(int i = 0; i <= j; i++)
        {
            reopt[s[i]][s[j]].C = 0;
            reopt[s[i]][s[j]].T = 0;
            reopt[s[i]][s[j]].W = j - i + 1;

            if(i == j)
                continue;
            

            for(int k = i; k < j; k++)
            {
                reopt[s[i]][s[j]].T += costM[s[k]][s[k+1]];
                reopt[s[i]][s[j]].C += costM[s[k]][s[k+1]] + reopt[s[i]][s[k]].T;
            }
        }
    }

    //Reverse (for 2opt)
    for(int i = 1; i < s.size(); i++)
    {
        for(int j = 1; j < i; j++)
        {
            reopt[s[i]][s[j]].C = 0;
            reopt[s[i]][s[j]].T = 0;
            reopt[s[i]][s[j]].W = i - j + 1;

            for(int k = i; k > j; k--)
            {
                reopt[s[i]][s[j]].T += costM[s[k]][s[k-1]];
            }
        }
        for(int j = 1; j < i; j++)
        {
            for(int k = i; k > j; k--)
            {
                reopt[s[i]][s[j]].C += costM[s[k]][s[k-1]] + reopt[s[i]][s[k]].T;
            }
        }
    }

    reopt[s[0]][s[s.size()-1]].W = s.size();
    reopt[s[0]][s[s.size()-1]].T = reopt[s[0]][s[s.size()-2]].T + reopt[s[s.size()-2]][s[s.size()-1]].T;
    reopt[s[0]][s[s.size()-1]].C = reopt[0][s.size()-2].C + reopt[s[0]][s.size()-2].T + reopt[s.size()-2][s.size()-1].T;

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
{ // O = (Ou, ... , Ov) = sq1 and O' = (O'w, ... , O'x) = sq2 are our subsequences
    ReoptData cost;

    cost.T = sq1.T + reopt[v][w].T + reopt[w][x].T;
    cost.C = sq1.C + reopt[w][x].W * (sq1.T + reopt[v][w].T) + reopt[w][x].C;
    cost.W = sq1.W + reopt[w][x].W;

    return cost;  
}

vector<int> swap (vector<int> s, vector<vector<ReoptData>> *reopt, bool *improved){
    int best_i = 0, best_j = 0;
    double bestDelta = 0, delta = 0;
    double preMoveCost = (*reopt)[s[0]][s[s.size()-1]].C;
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
                cost = concatCost(*reopt, sq[1], s[i], s[j+1], s[s.size()-1]);
                //reopt[0][i-1] o reopt[j][j] o reopt[i][i] o reopt[j+1][dimension]
            }else{
                sq[0] = concatCost(*reopt, s[0], s[i-1], s[j], s[j]);
                sq[1] = concatCost(*reopt, sq[0], s[j], s[i+1], s[j-1]);
                sq[2] = concatCost(*reopt, sq[1], s[j-1], s[i], s[i]);
                cost = concatCost(*reopt, sq[2], s[i], s[j+1], s[s.size()-1]);
                //reopt[0][i-1] o reopt[j][j] o reopt[i+1][j-1] o reopt[i][i] o reopt[j+1][dimension]
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
    double preMoveCost = (*reopt)[s[0]][s[s.size()-1]].C;
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
    double preMoveCost = (*reopt)[s[0]][s[s.size()-1]].C;
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
            *mainCost = costs[s[0]][s[s.size()-1]].C;
            improved = false;

            ngbhList = {N1/*, N2, N3, N4, N5*/};
            //Reopt update is done in movement functions
        }else{
            ngbhList.erase(std::remove(ngbhList.begin(), ngbhList.end(), ngbh_n), ngbhList.end());
        }
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

    for(int i = 0; i < I_MAX; i++)
    {
        solutionAlpha = construction();
        solutionBeta = solutionAlpha;

        costAlpha = solutionCost(solutionAlpha);
        costBeta = costAlpha;

        for(int iterILS = 0; iterILS < I_ILS; iterILS++){
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

    return 0;
}
