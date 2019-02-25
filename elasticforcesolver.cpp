#include "elasticforcesolver.h"

// This class encapsulates the First-Order Centred Scheme (FORCE) solver for nonlinear elastic systems, as detailed in Toro.
ElasticFORCESolver::ElasticFORCESolver()
{
}

// Computes the Lax-Friedrichs flux at the interface between two neighbouring cells, with states given by leftStateVector and rightStateVector respectively.
vector<double> ElasticFORCESolver::computeLaxFriedrichsFlux(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, double cellSpacing, double timeStep,
                                                            HyperelasticVariables hyperelasticVariables)
{
    vector<double> leftFluxVector = leftStateVector.computeFluxVector(hyperelasticVariables);
    vector<double> rightFluxVector = rightStateVector.computeFluxVector(hyperelasticVariables);

    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(hyperelasticVariables);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(hyperelasticVariables);

    return VectorAlgebra::multiplyVector(
                0.5, VectorAlgebra::addVectors(
                    VectorAlgebra::addVectors(
                        leftFluxVector, rightFluxVector), VectorAlgebra::multiplyVector(
                        (cellSpacing / timeStep), VectorAlgebra::subtractVectors(leftConservedVariableVector, rightConservedVariableVector))));
}

// Computes the Richtmyer flux at the interface between two neighbouring cells, with states given by leftStateVector and rightStateVector respectively.
vector<double> ElasticFORCESolver::computeRichtmyerFlux(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, double cellSpacing, double timeStep,
                                                        HyperelasticVariables hyperelasticVariables)
{
    vector<double> leftFluxVector = leftStateVector.computeFluxVector(hyperelasticVariables);
    vector<double> rightFluxVector = rightStateVector.computeFluxVector(hyperelasticVariables);

    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(hyperelasticVariables);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(hyperelasticVariables);

    vector<double> intermediateConservedVariableVector = VectorAlgebra::multiplyVector(
                0.5, VectorAlgebra::addVectors(
                    VectorAlgebra::addVectors(
                        leftConservedVariableVector, rightConservedVariableVector), VectorAlgebra::multiplyVector(
                        (timeStep / cellSpacing), VectorAlgebra::subtractVectors(leftFluxVector, rightFluxVector))));
    return ElasticStateVector::computeFluxVector(intermediateConservedVariableVector, hyperelasticVariables);
}

// Computes the FORCE flux at the interface between two neighbouring cells, with states given by leftStateVector and rightStateVector respectively.
vector<double> ElasticFORCESolver::computeFORCEFlux(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, double cellSpacing, double timeStep,
                                                    HyperelasticVariables hyperelasticVariables)
{
    vector<double> laxFriedrichsFlux = computeLaxFriedrichsFlux(leftStateVector, rightStateVector, cellSpacing, timeStep, hyperelasticVariables);
    vector<double> richtmyerFlux = computeRichtmyerFlux(leftStateVector, rightStateVector, cellSpacing, timeStep, hyperelasticVariables);

    return VectorAlgebra::multiplyVector(0.5, VectorAlgebra::addVectors(laxFriedrichsFlux, richtmyerFlux));
}

// Evolves the computational domain by one timestep, using the FORCE scheme.
void ElasticFORCESolver::computeFORCETimeStep(vector<ElasticStateVector> & newCells, vector<ElasticStateVector> & currentCells, double cellSpacing, double timeStep,
                                              HyperelasticVariables hyperelasticVariables)
{
    for (int i = 0; i < newCells.size(); i++)
    {
        vector<double> conservedVariableVector = newCells[i].computeConservedVariableVector(hyperelasticVariables);

        vector<double> leftFluxVector = computeFORCEFlux(currentCells[i], currentCells[i + 1], cellSpacing, timeStep, hyperelasticVariables);
        vector<double> rightFluxVector = computeFORCEFlux(currentCells[i + 1], currentCells[i + 2], cellSpacing, timeStep, hyperelasticVariables);

        newCells[i].setConservedVariableVector(
                    VectorAlgebra::addVectors(conservedVariableVector, VectorAlgebra::multiplyVector((timeStep / cellSpacing), VectorAlgebra::subtractVectors(leftFluxVector, rightFluxVector))),
                    hyperelasticVariables);
    }
}

// Evolves the computational domain until finalTime, using the FORCE scheme.
vector<ElasticStateVector> ElasticFORCESolver::solve(vector<ElasticStateVector> & cells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                     HyperelasticVariables hyperelasticVariables)
{
    double currentTime = 0;
    int currentIteration = 0;
    vector<ElasticStateVector> currentCells = cells;

    while (currentTime < finalTime)
    {
        vector<ElasticStateVector> currentCellsWithBoundary = ElasticSolvers::insertBoundaryCells(currentCells, 1);
        double timeStep = ElasticSolvers::computeStableTimeStep(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, hyperelasticVariables);

        computeFORCETimeStep(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, hyperelasticVariables);

        currentTime += timeStep;
        currentIteration += 1;
    }

    return currentCells;
}
