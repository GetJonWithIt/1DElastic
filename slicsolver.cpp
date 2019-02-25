#include "slicsolver.h"

// This class encapsulates the second-order Slope-Limiter Centred (SLIC) solver for nonlinear elastic systems, as detailed in Toro.
SLICSolver::SLICSolver()
{
}

// Computes the SLIC flux, using the interfaces between four neighbouring cells, with states given by leftLeftStateVector, leftStateVector, rightStateVector and rightRightStateVector,
// respectively.
vector<double> SLICSolver::computeSLICFlux(StateVector leftLeftStateVector, StateVector leftStateVector, StateVector rightStateVector, StateVector rightRightStateVector,
                                           double cellSpacing, double timeStep, double bias, int slopeLimiter)
{
    StateVector evolvedRightStateVector = Solvers::evolveForHalfTimeStep(leftLeftStateVector, leftStateVector, rightStateVector, cellSpacing, timeStep, bias, slopeLimiter, 1);
    StateVector evolvedLeftStateVector = Solvers::evolveForHalfTimeStep(leftStateVector, rightStateVector, rightRightStateVector, cellSpacing, timeStep, bias, slopeLimiter, 0);

    return FORCESolver::computeFORCEFlux(evolvedRightStateVector, evolvedLeftStateVector, cellSpacing, timeStep);
}

// Evolves the computational domain by one timestep, using the SLIC scheme.
void SLICSolver::computeSLICTimeStep(vector<StateVector> & newCells, vector<StateVector> & currentCells, double cellSpacing, double timeStep, double bias, int slopeLimiter)
{
    for (int i = 0; i < newCells.size(); i++)
    {
        vector<double> conservedVariableVector = newCells[i].computeConservedVariableVector();

        vector<double> leftFluxVector = computeSLICFlux(currentCells[i], currentCells[i + 1], currentCells[i + 2] , currentCells[i + 3], cellSpacing, timeStep, bias, slopeLimiter);
        vector<double> rightFluxVector = computeSLICFlux(currentCells[i + 1], currentCells[i + 2], currentCells[i + 3], currentCells[i + 4], cellSpacing, timeStep, bias, slopeLimiter);

        newCells[i].setConservedVariableVector(VectorAlgebra::addVectors(conservedVariableVector, VectorAlgebra::multiplyVector(
                                                                             (timeStep / cellSpacing), VectorAlgebra::subtractVectors(leftFluxVector, rightFluxVector))));
    }
}

// Evolves the computational domain until finalTime, using the SLIC scheme.
vector<StateVector> SLICSolver::solve(vector<StateVector> & cells, double cellSpacing, double CFLCoefficient, double finalTime, double bias, int slopeLimiter)
{
    double currentTime = 0;
    int currentIteration = 0;
    vector<StateVector> currentCells = cells;

    while (currentTime < finalTime)
    {
        vector<StateVector> currentCellsWithBoundary = Solvers::insertBoundaryCells(currentCells, 2);
        double timeStep = Solvers::computeStableTimeStep(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration);

        computeSLICTimeStep(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter);
        currentTime += timeStep;

        currentIteration += 1;
    }

    return currentCells;
}
