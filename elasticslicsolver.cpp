#include "elasticslicsolver.h"

// This class encapsulates the second-order Slope-Limiter Centred (SLIC) scheme for nonlinear elastic systems, as detailed in Toro.
ElasticSLICSolver::ElasticSLICSolver()
{
}

// Computes the SLIC flux, using the interfaces between four neighbouring cells, with states given by leftLeftStateVector, leftStateVector, rightStateVector and rightRightStateVector,
// respectively.
vector<double> ElasticSLICSolver::computeSLICFlux(ElasticStateVector leftLeftStateVector, ElasticStateVector leftStateVector, ElasticStateVector rightStateVector,
                                                  ElasticStateVector rightRightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                  HyperelasticVariables hyperelasticVariables)
{
    ElasticStateVector evolvedRightStateVector = ElasticSolvers::evolveForHalfTimeStep(leftLeftStateVector, leftStateVector, rightStateVector, cellSpacing, timeStep, bias,
                                                                                       slopeLimiter, 1, hyperelasticVariables);
    ElasticStateVector evolvedLeftStateVector = ElasticSolvers::evolveForHalfTimeStep(leftStateVector, rightStateVector, rightRightStateVector, cellSpacing, timeStep, bias,
                                                                                      slopeLimiter, 0, hyperelasticVariables);

    return ElasticFORCESolver::computeFORCEFlux(evolvedRightStateVector, evolvedLeftStateVector, cellSpacing, timeStep, hyperelasticVariables);
}

// Evolves the computational domain by one timestep, using the SLIC scheme.
void ElasticSLICSolver::computeSLICTimeStep(vector<ElasticStateVector> & newCells, vector<ElasticStateVector> & currentCells, double cellSpacing, double timeStep, double bias,
                                            int slopeLimiter, HyperelasticVariables hyperelasticVariables)
{
    for (int i = 0; i < newCells.size(); i++)
    {
        vector<double> conservedVariableVector = newCells[i].computeConservedVariableVector(hyperelasticVariables);

        vector<double> leftFluxVector = computeSLICFlux(currentCells[i], currentCells[i + 1], currentCells[i + 2], currentCells[i + 3], cellSpacing, timeStep, bias, slopeLimiter,
                hyperelasticVariables);
        vector<double> rightFluxVector = computeSLICFlux(currentCells[i + 1], currentCells[i + 2], currentCells[i + 3], currentCells[i + 4], cellSpacing, timeStep, bias, slopeLimiter,
                hyperelasticVariables);

        newCells[i].setConservedVariableVector(VectorAlgebra::addVectors(conservedVariableVector, VectorAlgebra::multiplyVector(
                                                                             (timeStep / cellSpacing), VectorAlgebra::subtractVectors(leftFluxVector, rightFluxVector))), hyperelasticVariables);
    }
}

// Evolves the computational domain until finalTime, using the SLIC scheme.
vector<ElasticStateVector> ElasticSLICSolver::solve(vector<ElasticStateVector> & cells, double cellSpacing, double CFLCoefficient, double finalTime, double bias, int slopeLimiter,
                                                    HyperelasticVariables hyperelasticVariables)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<ElasticStateVector> currentCells = cells;

    while (currentTime < finalTime)
    {
        vector<ElasticStateVector> currentCellsWithBoundary = ElasticSolvers::insertBoundaryCells(currentCells, 2);
        double timeStep = ElasticSolvers::computeStableTimeStep(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, hyperelasticVariables);

        computeSLICTimeStep(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, hyperelasticVariables);

        currentTime += timeStep;
        currentIteration += 1;
    }

    return currentCells;
}
