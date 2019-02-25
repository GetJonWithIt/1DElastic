#include "solvers.h"

// This class encapsulates the methods that are common to all solvers for the (stiffened gas) Euler equations.
Solvers::Solvers()
{
}

// Returns the computational domain, but with a given number of boundary cells inserted on each side (reflexive boundary conditions).
vector<StateVector> Solvers::insertBoundaryCells(vector<StateVector> & cells, int boundaryCells)
{
    int cellCount = cells.size();
    vector<StateVector> cellsWithBoundary(cellCount + (2 * boundaryCells));

    if (boundaryCells == 1)
    {
        cellsWithBoundary[0] = cells[0];
        cellsWithBoundary[cellCount + 1] = cells[cellCount - 1];

        cellsWithBoundary[0].setXVelocity(-cellsWithBoundary[0].getXVelocity());
        cellsWithBoundary[cellCount + 1].setXVelocity(-cellsWithBoundary[cellCount + 1].getXVelocity());
    }
    else if (boundaryCells == 2)
    {
        cellsWithBoundary[0] = cells[1];
        cellsWithBoundary[1] = cells[0];
        cellsWithBoundary[cellCount + 2] = cells[cellCount - 1];
        cellsWithBoundary[cellCount + 3] = cells[cellCount - 2];

        cellsWithBoundary[0].setXVelocity(-cellsWithBoundary[0].getXVelocity());
        cellsWithBoundary[1].setXVelocity(-cellsWithBoundary[1].getXVelocity());
        cellsWithBoundary[cellCount + 2].setXVelocity(-cellsWithBoundary[cellCount + 2].getXVelocity());
        cellsWithBoundary[cellCount + 3].setXVelocity(-cellsWithBoundary[cellCount + 3].getXVelocity());
    }

    for (int i = 0; i < cellCount; i++)
    {
        cellsWithBoundary[i + boundaryCells] = cells[i];
    }

    return cellsWithBoundary;
}

// Computes the estimated maximum wave speed across the entire computational domain, for the purposes of timestep calculation.
double Solvers::computeMaximumWaveSpeed(vector<StateVector> & cells)
{
    double maximumWaveSpeed = 0;

    for (int i = 0; i < cells.size(); i++)
    {
        double waveSpeed = abs(cells[i].getXVelocity()) + cells[i].computeSoundSpeed();
        if (waveSpeed > maximumWaveSpeed)
        {
            maximumWaveSpeed = waveSpeed;
        }
    }

    return maximumWaveSpeed;
}

// Computes the maximum stable timestep, in accordance with the CFL condition.
double Solvers::computeStableTimeStep(vector<StateVector> & cells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime, int currentIteration)
{
    double timeStep = CFLCoefficient * (cellSpacing / computeMaximumWaveSpeed(cells));

    if (currentIteration <= 5)
    {
        timeStep *= 0.2;
    }

    if (currentTime + timeStep > finalTime)
    {
        timeStep = finalTime - currentTime;
    }

    return timeStep;
}

// Evolves the state given by middleStateVector for half a timestep, using the slope across the neighbouring cells, given by leftStateVector, middleStateVector and rightStateVector,
// respectively.
StateVector Solvers::evolveForHalfTimeStep(StateVector leftStateVector, StateVector middleStateVector, StateVector rightStateVector, double cellSpacing, double timeStep, double bias,
                                           int slopeLimiter, int side)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVector(leftStateVector, middleStateVector, rightStateVector, bias, slopeLimiter);
    vector<double> leftExtrapolatedValue = VectorAlgebra::subtractVectors(middleStateVector.computeConservedVariableVector(), VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> rightExtrapolatedValue = VectorAlgebra::addVectors(middleStateVector.computeConservedVariableVector(), VectorAlgebra::multiplyVector(0.5, slopeVector));

    vector<double> fluxEvolution = VectorAlgebra::multiplyVector(0.5 * (timeStep / cellSpacing), VectorAlgebra::subtractVectors(StateVector::computeFluxVector(leftExtrapolatedValue),
                                                                                                                                StateVector::computeFluxVector(rightExtrapolatedValue)));

    StateVector evolvedStateVector;
    if (side == 0)
    {
        evolvedStateVector.setConservedVariableVector(VectorAlgebra::addVectors(leftExtrapolatedValue, fluxEvolution));
    }
    else
    {
        evolvedStateVector.setConservedVariableVector(VectorAlgebra::addVectors(rightExtrapolatedValue, fluxEvolution));
    }

    return evolvedStateVector;
}
