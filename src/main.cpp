#include "readData.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <chrono>
#include <limits>
#include <cmath>

#define LAST s.size() - 1

using namespace std;

enum NL{N1, N2, N3, N4, N5};

double **costM;
int dimension;

chrono::duration<double> timerReopt = chrono::system_clock::now() - chrono::system_clock::now();
chrono::duration<double> timerSwap = chrono::system_clock::now() - chrono::system_clock::now();
chrono::duration<double> timerRein = chrono::system_clock::now() - chrono::system_clock::now();
chrono::duration<double> timer2Opt = chrono::system_clock::now() - chrono::system_clock::now();

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

void calcReopt(vector<int> &s, vector<vector<ReoptData>> &reopt)
{
    auto timerStart = chrono::system_clock::now();

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

    auto timerEnd = chrono::system_clock::now();
    timerReopt += timerEnd - timerStart;
}

vector<int> swap (vector<int> s, vector<vector<ReoptData>> &reopt, bool *improved)
{
    auto timerStart = chrono::system_clock::now();

    int    best_i = 0, best_j = 0;
    double bestCost = reopt[0][LAST].C + reopt[0][LAST].T;
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

            if(cost + 1e-6 < bestCost){
                bestCost = cost;
                *improved = true;

                best_i = i;
                best_j = j;
            }
        }
    }

    //auto timerEnd = chrono::system_clock::now();

    if(*improved)
    {
        std::swap(s[best_j], s[best_i]);
        //timerEnd = chrono::system_clock::now();
        //calcReopt(s, reopt);

        for (int j = best_i; j <= best_j; j++)
        {
            for (int i = j - 1; i >= 0; i--)
            {
                reopt[i][j].W = reopt[i][j-1].W + reopt[j][j].W;
                reopt[i][j].T = reopt[i][j - 1].T + costM[s[j - 1]][s[j]];
                reopt[i][j].C = reopt[i][j-1].C + reopt[j][j].W * (reopt[i][j-1].T + costM[s[j - 1]][s[j]]) + reopt[j][j].C;

                reopt[j][i].W = reopt[i][j].W;
                reopt[j][i].T = reopt[i][j].T;

                reopt[j][i].C = reopt[j - 1][i].C + reopt[j - 1][i].W * (reopt[j][j].T + costM[s[j]][s[j - 1]]) + reopt[j][j].C;
            }
        }

        for (int i = best_i; i <= best_j; i++)
        {
            for (int j = best_j + 1; j <= dimension; j++)
            {
                reopt[i][j].W = reopt[i][j-1].W + reopt[j][j].W;
                reopt[i][j].T = reopt[i][j - 1].T + costM[s[j - 1]][s[j]];
                reopt[i][j].C = reopt[i][j-1].C + reopt[j][j].W * (reopt[i][j-1].T + costM[s[j - 1]][s[j]]) + reopt[j][j].C;

                reopt[j][i].W = reopt[i][j].W;
                reopt[j][i].T = reopt[i][j].T;

                reopt[j][i].C = reopt[j - 1][i].C + reopt[j - 1][i].W * (reopt[j][j].T + costM[s[j]][s[j - 1]]) + reopt[j][j].C;
            }
        }

        for (int i = best_i - 1, k = best_i, l; i >= 0; i--)
        {
            for (int j = best_j + 1; j <= dimension; j++)
            {

                reopt[i][j].W = reopt[i][k].W + reopt[k + 1][j].W;
                reopt[i][j].T = reopt[i][k].T + costM[s[k]][s[k + 1]] + reopt[k + 1][j].T;
                reopt[i][j].C = reopt[i][k].C + reopt[k + 1][j].W * (reopt[i][k].T + costM[s[k]][s[k + 1]]) + reopt[k + 1][j].C;

                reopt[j][i].W = reopt[i][j].W;
                reopt[j][i].T = reopt[i][j].T;

                reopt[j][i].C = reopt[k][i].C + reopt[k][i].W * (reopt[j][k + 1].T + costM[s[k + 1]][s[k]]) + reopt[j][k + 1].C;
            }
        }
    }

    auto timerEnd = chrono::system_clock::now();
    timerSwap += timerEnd - timerStart;


    return s;
}

vector<int> flip (vector<int> s, vector<vector<ReoptData>> &reopt, bool *improved)
{
    auto timerStart = chrono::system_clock::now();

    int best_i = 0, best_j = 0;
    double bestCost = reopt[0][LAST].C + reopt[0][LAST].T;
    ReoptData sq[2];
    double cost;

    for(int j = 3; j < LAST; j++)
    {
        for(int i = 1; i < j - 1; i++)
        {
            sq[0].W = reopt[0][i-1].W + reopt[j][i].W;
            sq[1].W = sq[0].W + reopt[j+1][LAST].W;
            
            sq[0].T = reopt[0][i-1].T + costM[s[i-1]][s[j]] + reopt[j][i].T;
            sq[1].T = sq[0].T + costM[s[i]][s[j+1]] + reopt[j+1][LAST].T;

            sq[0].C = reopt[0][i-1].C + reopt[j][i].W * ( reopt[0][i-1].T + costM[s[i-1]][s[j]] ) + reopt[j][i].C;
            sq[1].C = sq[0].C + reopt[j+1][LAST].W * ( sq[0].T + costM[s[i]][s[j+1]] ) + reopt[j+1][LAST].C;

            cost = sq[1].C + sq[1].T;
            //reopt[0][i-1] o reopt[j][i] o reopt[j+1][dimension]
            
            if(cost + 1e-6 < bestCost)
            {
                bestCost = cost;
                *improved = true;

                best_i = i;
                best_j = j;
            }
        }
    }

    auto timerEnd = chrono::system_clock::now();

    if(*improved)
    {
        std::reverse(s.begin() + best_i, s.begin() + best_j + 1);
        timerEnd = chrono::system_clock::now();
        calcReopt(s, reopt);
    }

    timer2Opt += timerEnd - timerStart;

    return s;
}

vector<int> reinsertion (vector<int> s, vector<vector<ReoptData>> &reopt, bool *improved, int subsegSize)
{
    auto timerStart = chrono::system_clock::now();

    int best_i = 0, best_j = 0;
    double bestCost = reopt[0][LAST].C + reopt[0][LAST].T;
    ReoptData sq[3];
    double cost;

    for(int i = 1; i < s.size() - subsegSize; i++)
    {
        for(int j = 1; j < s.size(); j++)
        {
            if(i <= j && j <= i + subsegSize) continue;

            if(i < j){
                sq[0].W = reopt[0][i-1].W + reopt[i+subsegSize][j-1].W;
                sq[1].W = sq[0].W + reopt[i][i+subsegSize-1].W;
                sq[2].W = sq[1].W + reopt[j][LAST].W;
                
                sq[0].T = reopt[0][i-1].T + costM[s[i-1]][s[i+subsegSize]] + reopt[i+subsegSize][j-1].T;
                sq[1].T = sq[0].T + costM[s[j-1]][s[i]] + reopt[i][i+subsegSize-1].T;
                sq[2].T = sq[1].T + costM[s[i+subsegSize-1]][s[j]] + reopt[j][LAST].T;

                sq[0].C = reopt[0][i-1].C + reopt[i+subsegSize][j-1].W * ( reopt[0][i-1].T + costM[s[i-1]][s[i+subsegSize]] ) + reopt[i+subsegSize][j-1].C;
                sq[1].C = sq[0].C + reopt[i][i+subsegSize-1].W * ( sq[0].T + costM[s[j-1]][s[i]] ) + reopt[i][i+subsegSize-1].C;
                sq[2].C = sq[1].C + reopt[j][LAST].W * ( sq[1].T + costM[s[i+subsegSize-1]][s[j]] ) + reopt[j][LAST].C;

                cost = sq[2].C + sq[2].T;
                //reopt[0][i-1] o reopt[i+subsegSize][j-1] o reopt[i][i+subsegSize-1] o reopt[j][dimension]
            }else{
                sq[0].W = reopt[0][j-1].W + reopt[i][i+subsegSize-1].W;
                sq[1].W = sq[0].W + reopt[j][i-1].W;
                sq[2].W = sq[1].W + reopt[i+subsegSize][LAST].W;
                
                sq[0].T = reopt[0][j-1].T + costM[s[j-1]][s[i]] + reopt[i][i+subsegSize-1].T;
                sq[1].T = sq[0].T + costM[s[i+subsegSize-1]][s[j]] + reopt[j][i-1].T;
                sq[2].T = sq[1].T + costM[s[i-1]][s[i+subsegSize]] + reopt[i+subsegSize][LAST].T;

                sq[0].C = reopt[0][j-1].C + reopt[i][i+subsegSize-1].W * ( reopt[0][j-1].T + costM[s[j-1]][s[i]] ) + reopt[i][i+subsegSize-1].C;
                sq[1].C = sq[0].C + reopt[j][i-1].W * ( sq[0].T + costM[s[i+subsegSize-1]][s[j]] ) + reopt[j][i-1].C;
                sq[2].C = sq[1].C + reopt[i+subsegSize][LAST].W * ( sq[1].T + costM[s[i-1]][s[i+subsegSize]] ) + reopt[i+subsegSize][LAST].C;

                cost = sq[2].C + sq[2].T;
                //reopt[0][j-1] o reopt[i][i+subsegSize-1] o reopt[j][i-1] o reopt[i+subsegSize][dimension]
            }
            
            if(cost + 1e-6 < bestCost){
                bestCost = cost;
                *improved = true;

                best_i = i;
                best_j = j;
            }
        }
    }

    auto timerEnd = chrono::system_clock::now();

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

        timerEnd = chrono::system_clock::now();

        calcReopt(s, reopt);
    }

    timerRein += timerEnd - timerStart;

    return s;
} 

void RVND (vector<int> &s, vector<vector<ReoptData>> &reopt, double *mainCost){
    vector<int> ngbhList = {N1, N2, N3, N4, N5};
    int ngbh_n;

    vector<int> neighbour_s = s;
    bool improved;

    while(!ngbhList.empty())
    {
        ngbh_n = ngbhList[rand() % ngbhList.size()];

        switch(ngbh_n){
            case N1:
                neighbour_s = swap(s, reopt, &improved);
                break;
            case N2:
                neighbour_s = flip(s, reopt, &improved);
                break;
            case N3:
                neighbour_s = reinsertion(s, reopt, &improved, 1);
                break;
            case N4:
                neighbour_s = reinsertion(s, reopt, &improved, 2);
                break;
            case N5:
                neighbour_s = reinsertion(s, reopt, &improved, 3);
                break; 
        }

        if(improved){
            s = neighbour_s;
            *mainCost = reopt[0][LAST].C + reopt[0][LAST].T;
            improved = false;
            ngbhList = {N1, N2 , N3, N4, N5};
            //Reopt update is done in movement functions
        }else{
            ngbhList.erase(std::remove(ngbhList.begin(), ngbhList.end(), ngbh_n), ngbhList.end());
        }
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
    //cout << "\tCOST MATRIX: \n";
    //printCostM();

    auto timerStart = chrono::system_clock::now();

    const int I_MAX = 10;
    const int I_ILS = (dimension >= 100) ? 100 : (dimension);

    vector<int> solutionAlpha, solutionBeta, solutionOmega;
    double costAlpha, costBeta, costOmega = numeric_limits<double>::infinity();

    vector<vector<ReoptData>> reopt(solutionAlpha.size(), vector<ReoptData>(solutionAlpha.size()));

    for(int i = 0; i < I_MAX; i++)
    {
        solutionAlpha = construction();
        solutionBeta = solutionAlpha;

        vector<vector<ReoptData>> reopt(solutionAlpha.size(), vector<ReoptData>(solutionAlpha.size()));
        calcReopt(solutionAlpha, reopt);

        costAlpha = reopt[0][dimension].C + reopt[0][dimension].T;
        costBeta = costAlpha;

        for(int iterILS = 0; iterILS < I_ILS; iterILS++)
        {
            RVND(solutionAlpha, reopt, &costAlpha);

            if(costAlpha < costBeta){
                solutionBeta = solutionAlpha;
                costBeta = costAlpha;
                iterILS = 0;
            }

            solutionAlpha = perturb(solutionBeta);
            calcReopt(solutionAlpha, reopt);
            costAlpha = reopt[0][dimension].C + reopt[0][dimension].T;
        }

        if(costBeta < costOmega){
            solutionOmega = solutionBeta;
            costOmega = costBeta;
        }
    }

    auto timerEnd = chrono::system_clock::now();
    chrono::duration<double> timerGilsRVND = timerEnd - timerStart;

    //PRINT COST AND SOLUTION
    
    /* std::cout << "SOLUTION:\n";
    for(auto k : solutionOmega){
        std::cout << k << ' ';
    } */
    std::cout << "\n\n" << "COST: " << costOmega << "\n";
    std::cout << "TIME: " << timerGilsRVND.count() << "\n";

    /* std::cout << "TIME REOPT: " << timerReopt.count() << "\n";
    std::cout << "TIME SWAP: " << timerSwap.count() << "\n";
    std::cout << "TIME 2OPT: " << timer2Opt.count() << "\n";
    std::cout << "TIME REIN: " << timerRein.count() << "\n"; */

    return 0;
}