#include "readData.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <chrono>
#include <limits>
#include <cmath>

#define MAX_COST 17976 // DEBUG
#define LAST s.size() - 1

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

double calculaCustoAcumulado(std::vector<int> &s)
{
  double custo = 0;
  for (int u = 0; u < dimension; u++)
  {
    custo += (dimension - u) * costM[s[u]][s[u + 1]];
  }

  return custo;
}

double solutionCost (vector <int> s)
{
    vector<vector<ReoptData>> reopt(s.size(), vector<ReoptData>(s.size()));
    
    // First assign value of individual elements
    for(int i = 0; i < s.size(); i++)
    {
        reopt[s[i]][s[i]].W = 1;
        reopt[s[i]][s[i]].T = 0;
        reopt[s[i]][s[i]].C = 0;    
    }
    reopt[s[0]][s[0]].W = 0; // Deposit has W = 0

    // Perform concatenation of [0][i-1] + [i][i] from 0 till end of s
    for(int i = 0; i < LAST; i++)
    {
        reopt[s[0]][s[i]].W = reopt[s[0]][s[i-1]].W + reopt[s[i]][s[i]].W;
        reopt[s[0]][s[i]].T = reopt[s[0]][s[i-1]].T + costM[s[i-1]][s[i]] /*+reopt[s[i]][s[i]].T(==0)*/;
        reopt[s[0]][s[i]].C = reopt[s[0]][s[i-1]].C + reopt[s[i]][s[i]].W*(reopt[s[0]][s[i-1]].T+costM[s[i-1]][s[i]]) /*+reopt[s[i]][s[i]].C(==0)*/; 
    }

    const int v = s.size()-2, w = LAST;

    reopt[0][0].W = reopt[s[0]][s[v]].W + reopt[s[w]][s[w]].W;
    reopt[0][0].T = reopt[s[0]][s[v]].T + costM[s[v]][s[w]] /*+reopt[s[w]][s[w]].T(==0)*/;
    reopt[0][0].C = reopt[s[0]][s[v]].C + reopt[s[w]][s[w]].W * (reopt[s[0]][s[v]].T+costM[s[v]][s[w]]) /*+reopt[s[w]][s[w]].C(==0)*/; 

    return reopt[0][0].C;
}

void calcReopt(vector<int> &s, vector<vector<ReoptData>> &reopt)
{
    for (int i = 0; i < dimension; i++)
    {
        reopt[i][i].T = 0;
        reopt[i][i].C = 0;
        reopt[i][i].W = 1;
    }
    reopt[0][0].W = 0;
    reopt[LAST][LAST].W = 0;

    for (int t = 2; t <= s.size(); t++)
    {
        for (int i = 0, j; i < dimension - t + 2; i++)
        {
            j = i + t - 1;

            reopt[i][j].W = reopt[i][j - 1].W + reopt[j][j].W;
            reopt[i][j].T = reopt[i][j - 1].T + costM[s[j - 1]][s[j]];
            reopt[i][j].C = reopt[i][j - 1].C + reopt[j][j].W * (reopt[i][j - 1].T + costM[s[j - 1]][s[j]]) + reopt[j][j].C;

            reopt[j][i].W = reopt[i][j].W;
            reopt[j][i].T = reopt[i][j].T;
            reopt[j][i].C = reopt[j - 1][i].C + reopt[j - 1][i].W * (reopt[j][j].T + costM[s[j]][s[j - 1]]) + reopt[j][j].C;
        }
    }

    if(reopt[0][LAST].C < MAX_COST){
        cout << "DEU RUIM" << endl;
    }
}

vector<int> swap (vector<int> s, vector<vector<ReoptData>> &reopt, bool *improved){
    int    best_i = 0, best_j = 0;
    double bestDelta = 0, delta = 0;
    double preMoveCost = reopt[0][LAST].C;
    ReoptData sq[4];
    double cost;

    for(int i = 1; i < LAST-1; i++)
    {
        for(int j = i+2; j < LAST; j++)
        {
            if(i == j - 1){ // If they're adjacent
                cout << "adjacent" << endl;

                sq[0].W = reopt[0][i-1].W + reopt[j][j].W;
                sq[1].W = sq[0].W + reopt[i][i].W;
                sq[2].W = sq[1].W + reopt[j+1][LAST].W;
                
                sq[0].T = reopt[0][i-1].T + costM[s[i-1]][s[j]] + reopt[j][j].T;
                sq[1].T = sq[0].T + costM[s[j]][s[i]] + reopt[i][i].T;
                sq[2].T = sq[1].T + costM[s[i]][s[j+1]] + reopt[j+1][LAST].T;

                sq[0].C = reopt[0][i-1].C + reopt[j][j].W * ( reopt[0][i-1].T + costM[s[i-1]][s[j]] ) + reopt[j][j].C;
                sq[1].C = sq[0].C + reopt[i][i].W * ( sq[0].T + costM[s[j]][s[i]] ) + reopt[i][i].C;
                sq[2].C = sq[1].C + reopt[j+1][LAST].W * ( sq[1].T + costM[s[i]][s[j+1]] ) + reopt[j+1][LAST].C;

                cost = sq[2].C + sq[2].T;
                //reopt[0][i-1] o reopt[j][j] o reopt[i][i] o reopt[j+1][dimension]
            }else{ 
                sq[0].W = reopt[0][i-1].W + reopt[j][j].W;
                sq[1].W = sq[0].W + reopt[i+1][j-1].W;
                sq[2].W = sq[1].W + reopt[i][i].W;
                sq[3].W = sq[2].W + reopt[j+1][LAST].W;
                
                sq[0].T = reopt[0][i-1].T + costM[s[i-1]][s[j]] + reopt[j][j].T;
                sq[1].T = sq[0].T + costM[s[j]][s[i+1]] + reopt[i+1][j-1].T;
                sq[2].T = sq[1].T + costM[s[j-1]][s[i]] + reopt[i][i].T;
                sq[3].T = sq[2].T + costM[s[i]][s[j+1]] + reopt[j+1][LAST].T;

                sq[0].C = reopt[0][i-1].C + reopt[j][j].W * ( reopt[0][i-1].T + costM[s[i-1]][s[j]] ) + reopt[j][j].C;
                sq[1].C = sq[0].C + reopt[i+1][j-1].W * ( sq[0].T + costM[s[j]][s[i+1]] ) + reopt[i+1][j-1].C;
                sq[2].C = sq[1].C + reopt[i][i].W * ( sq[1].T + costM[s[j-1]][s[i]] ) + reopt[i][j].C;
                sq[3].C = sq[2].C + reopt[j+1][LAST].W * ( sq[2].T + costM[s[i]][s[j+1]] ) + reopt[j+1][LAST].C;

                cost = sq[3].C + sq[3].T;
                //reopt[0][i-1] o reopt[j][j] o reopt[i+1][j-1] o reopt[i][i] o reopt[j+1][dimension]
            } 

            delta = cost - preMoveCost;
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
        calcReopt(s, reopt);
    }

    return s;
}

vector<int> flip (vector<int> s, vector<vector<ReoptData>> &reopt, bool *improved){
    int best_i = 0, best_j = 0;
    double bestDelta = 0, delta = 0;
    double preMoveCost = reopt[0][0].C;
    ReoptData sq[2];
    double cost;

    for(int j = 3; j < LAST; j++)
    {
        for(int i = 1; i < j - 1; i++)
        {
            sq[0].W = reopt[s[0]][s[i-1]].W + reopt[s[j]][s[i]].W;
            sq[1].W = sq[0].W + reopt[s[j+1]][s[LAST]].W;
            
            sq[0].T = reopt[s[0]][s[i-1]].T + costM[s[i-1]][s[j]] + reopt[s[j]][s[i]].T;
            sq[1].T = sq[0].T + costM[s[i]][s[j+1]] + reopt[s[j+1]][s[LAST]].T;

            sq[0].C = reopt[s[0]][s[i-1]].C + reopt[s[j]][s[i]].W * ( reopt[s[0]][s[i-1]].T + costM[s[i-1]][s[j]] ) + reopt[s[j]][s[i]].C;
            sq[1].C = sq[0].C + reopt[s[j+1]][s[LAST]].W * ( sq[0].T + costM[s[i]][s[j+1]] ) + reopt[s[j+1]][s[LAST]].C;


            cost = sq[1].C + sq[1].T;
            //reopt[0][i-1] o reopt[j][i] o reopt[j+1][dimension]
            
            delta = cost - preMoveCost;
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
        calcReopt(s, reopt);
    }

    return s;
}

vector<int> reinsertion (vector<int> s, vector<vector<ReoptData>> &reopt, bool *improved, int subsegSize){
    int best_i = 0, best_j = 0;
    double bestDelta = 0, delta = 0;
    double preMoveCost = reopt[0][0].C;
    ReoptData sq[3];
    double cost;

    for(int i = 1; i < s.size() - subsegSize; i++)
    {
        for(int j = 1; j < s.size(); j++)
        {
            if(i <= j && j <= i + subsegSize) continue;

            if(i < j){
                sq[0].W = reopt[s[0]][s[i-1]].W + reopt[s[i+subsegSize]][s[j-1]].W;
                sq[1].W = sq[0].W + reopt[s[i]][s[i+subsegSize-1]].W;
                sq[2].W = sq[1].W + reopt[s[j]][s[LAST]].W;
                
                sq[0].T = reopt[s[0]][s[i-1]].T + costM[s[i-1]][s[i+subsegSize]] + reopt[s[i+subsegSize]][s[j-1]].T;
                sq[1].T = sq[0].T + costM[s[j-1]][s[i]] + reopt[s[i]][s[i+subsegSize-1]].T;
                sq[2].T = sq[1].T + costM[s[i+subsegSize-1]][s[j]] + reopt[s[j]][s[LAST]].T;

                sq[0].C = reopt[s[0]][s[i-1]].C + reopt[s[i+subsegSize]][s[j-1]].W * ( reopt[s[0]][s[i-1]].T + costM[s[i-1]][s[i+subsegSize]] ) + reopt[s[i+subsegSize]][s[j-1]].C;
                sq[1].C = sq[0].C + reopt[s[i]][s[i+subsegSize-1]].W * ( sq[0].T + costM[s[j-1]][s[i]] ) + reopt[s[i]][s[i+subsegSize-1]].C;
                sq[2].C = sq[1].C + reopt[s[j]][s[LAST]].W * ( sq[1].T + costM[s[i+subsegSize-1]][s[j]] ) + reopt[s[j]][s[LAST]].C;

                cost = sq[2].C + sq[2].T;
                //reopt[0][i-1] o reopt[i+subsegSize][j-1] o reopt[i][i+subsegSize-1] o reopt[j][dimension]
            }else{
                sq[0].W = reopt[s[0]][s[j-1]].W + reopt[s[i]][s[i+subsegSize-1]].W;
                sq[1].W = sq[0].W + reopt[s[j]][s[i-1]].W;
                sq[2].W = sq[1].W + reopt[s[i+subsegSize]][s[LAST]].W;
                
                sq[0].T = reopt[s[0]][s[j-1]].T + costM[s[j-1]][s[i]] + reopt[s[i]][s[i+subsegSize-1]].T;
                sq[1].T = sq[0].T + costM[s[i+subsegSize-1]][s[j]] + reopt[s[j]][s[i-1]].T;
                sq[2].T = sq[1].T + costM[s[i-1]][s[i+subsegSize]] + reopt[s[i+subsegSize]][s[LAST]].T;

                sq[0].C = reopt[s[0]][s[j-1]].C + reopt[s[i]][s[i+subsegSize-1]].W * ( reopt[s[0]][s[j-1]].T + costM[s[j-1]][s[i]] ) + reopt[s[i]][s[i+subsegSize-1]].C;
                sq[1].C = sq[0].C + reopt[s[j]][s[i-1]].W * ( sq[0].T + costM[s[i+subsegSize-1]][s[j]] ) + reopt[s[j]][s[i-1]].C;
                sq[2].C = sq[1].C + reopt[s[i+subsegSize]][s[LAST]].W * ( sq[1].T + costM[s[i-1]][s[i+subsegSize]] ) + reopt[s[i+subsegSize]][s[LAST]].C;

                cost = sq[2].C + sq[2].T;
                //reopt[0][j-1] o reopt[i][i+subsegSize-1] o reopt[j][i-1] o reopt[i+subsegSize][dimension]
            }
            
            delta = cost - preMoveCost;
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

        calcReopt(s, reopt);
    }

    return s;
} 

void RVND (vector<int> &s, double *mainCost){
    vector<int> ngbhList = {N1, N2, N3, N4, N5};
    int ngbh_n;

    vector<int> neighbour_s = s;
    bool improved;

    vector<vector<ReoptData>> costs(s.size(), vector<ReoptData>(s.size()));
    calcReopt(s, costs);

    int i = 0; //DEBUG

    while(!ngbhList.empty())
    {
        ngbh_n = ngbhList[rand() % ngbhList.size()];

        switch(ngbh_n){
            case N1:
                neighbour_s = swap(s, costs, &improved);
                break;
            case N2:
                neighbour_s = flip(s, costs, &improved);
                break;
            case N3:
                neighbour_s = reinsertion(s, costs, &improved, 1);
                break;
            case N4:
                neighbour_s = reinsertion(s, costs, &improved, 2);
                break;
            case N5:
                neighbour_s = reinsertion(s, costs, &improved, 3);
                break; 
        }

        if(improved){
            s = neighbour_s;
            *mainCost = costs[0][LAST].C;
            improved = false;

            /* if(*mainCost < MAX_COST) // DEBUG
            {
                cout << "movement: " << ngbh_n << endl;
            } */

            //cout << "loop #" << i << '\t' << "IMPROVED " << (*mainCost) << endl; //DEBUG

            ngbhList = {N1, N2, N3, N4, N5};
            //Reopt update is done in movement functions
        }else{
            //cout << "loop #" << i << '\t' << "!improved" << endl; //DEBUG

            ngbhList.erase(std::remove(ngbhList.begin(), ngbhList.end(), ngbh_n), ngbhList.end());
        }

        i++; //DEBUG
    }
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
        subseg1Start = randomRange(1, LAST - maxSubsegSize);
        subseg1End = subseg1Start + randomRange(minSubsegSize, maxSubsegSize);

        subseg2Start = randomRange(1, LAST - maxSubsegSize);
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
        costAlpha = calculaCustoAcumulado(solutionAlpha);

        /* if(costAlpha != costAlpha2) // DEBUG
        {
            cout << "SOLUTIONCOST ESTA ERRADO" << endl;
        } */

        if(i == 0){ //DEBUG
            firstCost = costAlpha;
        }

        costBeta = costAlpha;

        for(int iterILS = 0; iterILS < I_ILS; iterILS++){

            //cout << "Iter " << i << " of " << I_MAX << endl; //DEBUG
            //cout << "IterILS " << iterILS << " of " << I_ILS << endl; //DEBUG

            RVND(solutionAlpha, &costAlpha);

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