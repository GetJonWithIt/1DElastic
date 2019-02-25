#ifndef SLICSOLVER_H
#define SLICSOLVER_H

#include "forcesolver.h"

class SLICSolver
{
public:
    SLICSolver();

    static vector<double> computeSLICFlux(StateVector leftLeftStateVector, StateVector leftStateVector, StateVector rightStateVector, StateVector rightRightStateVector, double cellSpacing,
                                   double timeStep, double bias, int slopeLimiter);

    static void computeSLICTimeStep(vector<StateVector> & newCells, vector<StateVector> & currentCells, double cellSpacing, double timeStep, double bias, int slopeLimiter);
    static vector<StateVector> solve(vector<StateVector> & cells, double cellSpacing, double CFLCoefficient, double finalTime, double bias, int slopeLimiter);
};

#endif // SLICSOLVER_H
