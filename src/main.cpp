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

double ** costM;
int dimension;

typedef struct{
    vector<vector<double>> C;
    vector<vector<double>> T;
    vector<vector<int>> W;
} ReoptData;

void printCostM() {
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

    int j = 0;
    while(!candidateList.empty())
    {
        int closestNode;
        double closestNodeDistance = std::numeric_limits<double>::max();

        for(int k : candidateList)
        {
            if(costM[solution[j]][k] < closestNodeDistance)
            {
                closestNodeDistance = costM[solution[j]][k];
                closestNode = k;
            }    
        }

        solution.push_back(closestNode);
        candidateList.erase(std::remove(candidateList.begin(), candidateList.end(), closestNode), candidateList.end());

        j++;
    }

    solution.push_back(1);

    return solution;
}

double solutionCost (vector <int> solution){
    double cost = 0;
    for(int i = 0; i < solution.size() - 1; i++){
        cost += costM[solution[i]][solution[i+1]];
    }
    return cost;
}

ReoptData calcReopt(ReoptData reopt, vector<int> s, int ii, int jj){
    int small = (ii <= jj)? ii : jj;
    
    for(int j = small; j < dimension; j++){
        for(int i = small; i <= j; i++){
            reopt.C[i][j] = 0;
            reopt.T[i][j] = 0;
            reopt.W[i][j] = j - i + 1;

            if(i == j){
                continue;
            }

            for(int k = i; k < j; k++){
                reopt.C[s[i]][s[j]] += costM[s[k]][s[k+1]];
                reopt.T[s[i]][s[j]] += costM[s[k]][s[k+1]] + reopt.T[s[i]][s[j]];
            }
        }
    }

    for(int i = small; i < dimension; i++){
        for(int j = small; j < i; j++){
            reopt.C[i][j] = 0;
            reopt.T[i][j] = 0;
            reopt.W[i][j] = i - j + 1;

            for(int k = i, l = j; l < i; k--, l++){
                reopt.C[s[i]][s[j]] += costM[s[k]][s[k-1]];
                reopt.T[s[i]][s[j]] += costM[s[k]][s[k-1]] + reopt.T[s[i]][s[j]];
            }
        }
    }

    return reopt;
}

ReoptData calcReopt(ReoptData reopt, vector<int> s){
    reopt.W[0][0] = 0;
    reopt.T[0][0] = 0;
    reopt.C[0][0] = 0;
    
    for(int j = 1; j < dimension; j++){
        for(int i = 0; i <= j; i++){
            reopt.C[i][j] = 0;
            reopt.T[i][j] = 0;
            reopt.W[i][j] = j - i + 1;

            if(i == j){
                continue;
            }

            for(int k = i; k < j; k++){
                reopt.C[s[i]][s[j]] += costM[s[k]][s[k+1]];
                reopt.T[s[i]][s[j]] += costM[s[k]][s[k+1]] + reopt.T[s[i]][s[j]]; //Checar se esse calculo ta certo
            }
        }
    }

    //Reverse (for 2opt)
    for(int i = 1; i < dimension; i++){
        for(int j = 1; j < i; j++){
            reopt.C[i][j] = 0;
            reopt.T[i][j] = 0;
            reopt.W[i][j] = i - j + 1;

            for(int k = i, l = j; l < i; k--, l++){
                reopt.C[s[i]][s[j]] += costM[s[k]][s[k-1]];
                reopt.T[s[i]][s[j]] += costM[s[k]][s[k-1]] + reopt.T[s[i]][s[j]]; //Checar se esse calculo ta certo (se o reverso esta realmente sendo contabilizado)
            }
        }
    }

    return reopt;
}

/*swap:
    i    j
 01|2|34|5|60
 01|5|34|2|60
reopt[0][i-1] o reopt[j][j] o reopt[i+1][j-1] o reopt[i][i] o reopt[j+1][dimension]

flip:
     i j
 012|345|60
 012|543|60
reopt[0][i-1] o reopt[j][i] o reopt[j+1][dimension]

reinsertion:
     i       j
 012|345|678|90
 012|678|345|90
reopt[0][i-1] o reopt[i+subsegSize][j-1] o reopt[i][i+subsegSize-1] o reopt[j][dimension]

double bridge:
   1s  1e 2s 2e
 01|234|56|78|90
 01|87|56|432|90
reopt[0][subseg1Start-1] o reopt[subseg2End-1][subseg2Start] o reopt[subseg1End][subseg2Start-1] o reopt[subseg1End-1][subseg1Start] o reopt[subseg2End][dimension]

*/

vector<int> swap (vector<int> s, double *bestDelta){
    int best_i = 0, best_j = 0;
    *bestDelta = 0; //0 = no change
    double delta = 0, semiDelta = 0;

    for(int j = 2; j < s.size() - 1; j++)
    {
        semiDelta =  -costM[s[j]][s[j-1]] -costM[s[j]][s[j+1]]; //relies only on j

        for(int i = 1; i < j; i++)
        {
            if(i == j - 1){ // If they're adjacent
                delta = costM[s[j]][s[i - 1]] + costM[s[i]][s[j + 1]] - costM[s[i]][s[i - 1]] - costM[s[j]][s[j + 1]];
            }else{
                delta = costM[s[i]][s[j-1]] + costM[s[i]][s[j+1]] + costM[s[j]][s[i-1]] + costM[s[j]][s[i+1]]
                - costM[s[i]][s[i-1]] - costM[s[i]][s[i+1]] + semiDelta;
            }

            if(delta < *bestDelta){
                *bestDelta = delta;

                best_i = i;
                best_j = j;
            }
        }
    }

    if(*bestDelta < 0){
        std::swap(s[best_j], s[best_i]);
    }else{
        *bestDelta = 0;
    }

    return s;
}

vector<int> flip (vector<int> s, double *bestDelta){
    int best_i = 0, best_j = 0;
    *bestDelta = 0; //0 = no change
    double delta = 0, semiDelta = 0;

    for(int j = 3; j < s.size() - 1; j++)
    {
        semiDelta = -costM[s[j]][s[j+1]]; //relies only on j

        for(int i = 1; i < j - 1; i++)
        {
            delta = -costM[s[i-1]][s[i]] + semiDelta + costM[s[i]][s[j+1]] +costM[s[j]][s[i-1]];

            if(delta < *bestDelta - std::numeric_limits<double>::epsilon())
            {
                *bestDelta = delta;

                best_i = i;
                best_j = j;
            }
        }
    }

    if(*bestDelta < std::numeric_limits<double>::epsilon()){
        std::reverse(s.begin() + best_i, s.begin() + best_j + 1);
    }else{
        *bestDelta = 0;
    }

    return s;
}

vector<int> reinsertion (vector<int> s, double *bestDelta, int subsegSize){
    int best_i = 0, best_j = 0;
    *bestDelta = 0; //0 = no change
    double delta = 0, semiDelta = 0;

    for(int i = 1; i < s.size() - subsegSize; i++)
    {
        semiDelta = costM[s[i - 1]][s[i + subsegSize]] - costM[s[i]][s[i - 1]] - costM[s[i + subsegSize - 1]][s[i + subsegSize]];

        for(int j = 1; j < s.size(); j++)
        {
            if(i <= j && j <= i + subsegSize) continue;

            delta = costM[s[j - 1]][s[i]] + costM[s[i + subsegSize - 1]][s[j]] - costM[s[j]][s[j - 1]] + semiDelta;

            if(delta < *bestDelta - std::numeric_limits<double>::epsilon()){
                *bestDelta = delta;

                best_i = i;
                best_j = j;
            }

        }
    }

    if(*bestDelta < std::numeric_limits<double>::epsilon())
    {
        vector<int> subseg(s.begin() + best_i, s.begin() + best_i + subsegSize);

        if(best_i < best_j){
            s.insert(s.begin() + best_j, subseg.begin(), subseg.end());
            s.erase(s.begin() + best_i, s.begin() + best_i + subsegSize);
        }else{
            s.erase(s.begin() + best_i, s.begin() + best_i + subsegSize);
            s.insert(s.begin() + best_j, subseg.begin(), subseg.end());
        }

    }else{
        *bestDelta = 0;
    }

    return s;
}

vector<int> RVND (vector<int> s, double *mainCost){
    vector<int> ngbhList = {N1, N2, N3, N4, N5};
    int ngbh_n;

    vector<int> neighbour_s = s;
    double neighbourCost = numeric_limits<double>::max();
    double delta = 0; //Receives delta from neighbouthood movements

    while(!ngbhList.empty())
    {
        ngbh_n = ngbhList[rand() % ngbhList.size()];

        switch(ngbh_n){
            case N1:
                neighbour_s = swap(s, &delta);
                break;
            case N2:
                neighbour_s = flip(s, &delta);
                break;
            case N3:
                neighbour_s = reinsertion(s, &delta, 1);
                break;
            case N4:
                neighbour_s = reinsertion(s, &delta, 2);
                break;
            case N5:
                neighbour_s = reinsertion(s, &delta, 3);
                break;
        }

        neighbourCost = *mainCost + delta;

        if(neighbourCost < *mainCost){
            s = neighbour_s;
            *mainCost = neighbourCost;

            ngbhList = {N1, N2, N3, N4, N5};
        }else{
            ngbhList.erase(std::remove(ngbhList.begin(), ngbhList.end(), ngbh_n), ngbhList.end());
        }
    }

    return s;
}

int randomRange(int min, int max){ // Inclusive
    return min + (rand() % (max - min + 1));
}

vector<int> perturb (vector<int> s, double costIn, double *costOut){
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

    double delta = -costM[s[subseg1Start - 1]][s[subseg1Start]] -costM[s[subseg1End - 1]][s[subseg1End]]
                   -costM[s[subseg2Start - 1]][s[subseg2Start]] -costM[s[subseg2End - 1]][s[subseg2End]]
                   +costM[s[subseg1Start - 1]][s[subseg2End - 1]] +costM[s[subseg2Start]][s[subseg1End]]
                   +costM[s[subseg2Start - 1]][s[subseg1End - 1]] +costM[s[subseg1Start]][s[subseg2End]];

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

    *costOut = costIn + delta;

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

            solutionAlpha = perturb(solutionBeta, costBeta, &costAlpha);
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
