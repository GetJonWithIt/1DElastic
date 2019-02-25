#ifndef SOLVERS_H
#define SOLVERS_H

#include "statevector.h"
#include "vectoralgebra.h"
#include "slopelimiters.h"
using namespace std;

class Solvers
{
public:
    Solvers();

    static vector<StateVector> insertBoundaryCells(vector<StateVector> & cells, int boundaryCells);

    static double computeMaximumWaveSpeed(vector<StateVector> & cells);
    static double computeStableTimeStep(vector<StateVector> & cells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime, int currentIteration);

    static StateVector evolveForHalfTimeStep(StateVector leftStateVector, StateVector middleStateVector, StateVector rightStateVector, double cellSpacing, double timeStep, double bias,
                                             int slopeLimiter, int side);
};

#endif // SOLVERS_H
