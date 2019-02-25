#ifndef FORCESOLVER_H
#define FORCESOLVER_H

#include "statevector.h"
#include "vectoralgebra.h"
#include "solvers.h"
using namespace std;

class FORCESolver
{
public:
    FORCESolver();

    static vector<double> computeLaxFriedrichsFlux(StateVector leftStateVector, StateVector rightStateVector, double cellSpacing, double timeStep);
    static vector<double> computeRichtmyerFlux(StateVector leftStateVector, StateVector rightStateVector, double cellSpacing, double timeStep);

    static vector<double> computeFORCEFlux(StateVector leftStateVector, StateVector rightStateVector, double cellSpacing, double timeStep);

    static void computeFORCETimeStep(vector<StateVector> & newCells, vector<StateVector> & currentCells, double cellSpacing, double timeStep);
    static vector<StateVector> solve(vector<StateVector> & cells, double cellSpacing, double CFLCoefficient, double finalTime);
};

#endif // FORCESOLVER_H
