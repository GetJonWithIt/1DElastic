#include "forcesolver.h"

// This class encapsulates the First-Order Centred Scheme (FORCE) solver for the (stiffened gas) Euler equations, as detailed in Toro.
FORCESolver::FORCESolver()
{
}

// Computes the Lax-Friedrichs flux at the interface between two neighbouring cells, with states given by leftStateVector and rightStateVector respectively.
vector<double> FORCESolver::computeLaxFriedrichsFlux(StateVector leftStateVector, StateVector rightStateVector, double cellSpacing, double timeStep)
{
    vector<double> leftFluxVector = leftStateVector.computeFluxVector();
    vector<double> rightFluxVector = rightStateVector.computeFluxVector();

    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector();
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector();

    return VectorAlgebra::multiplyVector(
                0.5, VectorAlgebra::addVectors(
                    VectorAlgebra::addVectors(
                        leftFluxVector, rightFluxVector), VectorAlgebra::multiplyVector(
                        (cellSpacing / timeStep), VectorAlgebra::subtractVectors(leftConservedVariableVector, rightConservedVariableVector))));
}

// Computes the Richtmyer flux at the interface between two neighbouring cells, with states given by leftStateVector and rightStateVector respectively.
vector<double> FORCESolver::computeRichtmyerFlux(StateVector leftStateVector, StateVector rightStateVector, double cellSpacing, double timeStep)
{
    vector<double> leftFluxVector = leftStateVector.computeFluxVector();
    vector<double> rightFluxVector = rightStateVector.computeFluxVector();

    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector();
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector();

    vector<double> intermediateConservedVariableVector = VectorAlgebra::multiplyVector(
                0.5, VectorAlgebra::addVectors(
                    VectorAlgebra::addVectors(
                        leftConservedVariableVector, rightConservedVariableVector), VectorAlgebra::multiplyVector(
                        (timeStep / cellSpacing), VectorAlgebra::subtractVectors(leftFluxVector, rightFluxVector))));
    return StateVector::computeFluxVector(intermediateConservedVariableVector);
}

// Computes the FORCE flux at the interface between two neighbouring cells, with states given by leftStateVector and rightStateVector respectively.
vector<double> FORCESolver::computeFORCEFlux(StateVector leftStateVector, StateVector rightStateVector, double cellSpacing, double timeStep)
{
    vector<double> laxFriedrichsFlux = computeLaxFriedrichsFlux(leftStateVector, rightStateVector, cellSpacing, timeStep);
    vector<double> richtmyerFlux = computeRichtmyerFlux(leftStateVector, rightStateVector, cellSpacing, timeStep);

    return VectorAlgebra::multiplyVector(0.5, VectorAlgebra::addVectors(laxFriedrichsFlux, richtmyerFlux));
}

// Evolves the computation domain by one timestep, using the FORCE scheme.
void FORCESolver::computeFORCETimeStep(vector<StateVector> & newCells, vector<StateVector> & currentCells, double cellSpacing, double timeStep)
{
    for (int i = 0; i < newCells.size(); i++)
    {
        vector<double> conservedVariableVector = newCells[i].computeConservedVariableVector();

        vector<double> leftFluxVector = computeFORCEFlux(currentCells[i], currentCells[i + 1], cellSpacing, timeStep);
        vector<double> rightFluxVector = computeFORCEFlux(currentCells[i + 1], currentCells[i + 2], cellSpacing, timeStep);

        newCells[i].setConservedVariableVector(
                    VectorAlgebra::addVectors(conservedVariableVector, VectorAlgebra::multiplyVector((timeStep / cellSpacing), VectorAlgebra::subtractVectors(leftFluxVector, rightFluxVector))));
    }
}

// Evolves the computational domain until finalTime, using the FORCE scheme.
vector<StateVector> FORCESolver::solve(vector<StateVector> & cells, double cellSpacing, double CFLCoefficient, double finalTime)
{
    double currentTime = 0;
    int currentIteration = 0;
    vector<StateVector> currentCells = cells;

    while (currentTime < finalTime)
    {
        vector<StateVector> currentCellsWithBoundary = Solvers::insertBoundaryCells(currentCells, 1);
        double timeStep = Solvers::computeStableTimeStep(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration);

        computeFORCETimeStep(currentCells, currentCellsWithBoundary, cellSpacing, timeStep);
        currentTime += timeStep;

        currentIteration += 1;
    }

    return currentCells;
}
