#include "statevector.h"

// This class encapsulates the complete state vector for a (stiffened gas) Euler system, as detailed in Toro.

// Constructs a state vector with default values for the density, x/y/z velocities, pressure, adiabatic exponent, and stiffening parameter.
StateVector::StateVector()
{
    density = 1;
    xVelocity = 0;
    yVelocity = 0;
    zVelocity = 0;
    pressure = 1;
    adiabaticIndex = 1;
    stiffeningParameter = 0;
}

// Constructs a state vector with specified values for the density, x/y/z velocities, pressure, adiabatic exponent, and stiffening parameter.
StateVector::StateVector(double newDensity, double newXVelocity, double newYVelocity, double newZVelocity, double newPressure, double newAdiabaticIndex, double newStiffeningParameter)
{
    density = newDensity;
    xVelocity = newXVelocity;
    yVelocity = newYVelocity;
    zVelocity = newZVelocity;
    pressure = newPressure;
    adiabaticIndex = newAdiabaticIndex;
    stiffeningParameter = newStiffeningParameter;
}

// Computes the specific internal energy for a stiffened gas.
double StateVector::computeSpecificInternalEnergy()
{
    return EquationOfState::computeSpecificInternalEnergy(density, pressure, adiabaticIndex, stiffeningParameter);
}

// Computes the total energy (sum of internal and kinetic energies) for a stiffened gas.
double StateVector::computeTotalEnergy()
{
    double velocitySquared = (xVelocity * xVelocity) + (yVelocity * yVelocity) + (zVelocity * zVelocity);
    return density * ((velocitySquared / 2) + computeSpecificInternalEnergy());
}

// Computes the sound speed for a stiffened gas.
double StateVector::computeSoundSpeed()
{
    return EquationOfState::computeSoundSpeed(density, pressure, adiabaticIndex, stiffeningParameter);
}

// Computes the vector of primitive variables (i.e. density, x/y/z velocities, pressure, adiabatic exponent, and stiffening parameter).
vector<double> StateVector::computePrimitiveVariableVector()
{
    vector<double> primitiveVariableVector(7);

    primitiveVariableVector[0] = density;
    primitiveVariableVector[1] = xVelocity;
    primitiveVariableVector[2] = yVelocity;
    primitiveVariableVector[3] = zVelocity;
    primitiveVariableVector[4] = pressure;
    primitiveVariableVector[5] = adiabaticIndex;
    primitiveVariableVector[6] = stiffeningParameter;

    return primitiveVariableVector;
}

// Computes the vector of conserved variables (i.e. density, x/y/z momenta, total energy, and conserved forms of the adiabatic exponent and stiffening parameter).
vector<double> StateVector::computeConservedVariableVector()
{
    vector<double> conservedVariableVector(7);

    conservedVariableVector[0] = density;
    conservedVariableVector[1] = density * xVelocity;
    conservedVariableVector[2] = density * yVelocity;
    conservedVariableVector[3] = density * zVelocity;
    conservedVariableVector[4] = computeTotalEnergy();
    conservedVariableVector[5] = density * adiabaticIndex;
    conservedVariableVector[6] = density * stiffeningParameter;

    return conservedVariableVector;
}

// Computes the flux vector in the x-direction, using the vector of conserved variables.
vector<double> StateVector::computeFluxVector(vector<double> conservedVariableVector)
{
    vector<double> fluxVector(7);

    double computedDensity = conservedVariableVector[0];
    double computedXVelocity = conservedVariableVector[1] / computedDensity;
    double computedYVelocity = conservedVariableVector[2] / computedDensity;
    double computedZVelocity = conservedVariableVector[3] / computedDensity;
    double computedAdiabaticIndex = conservedVariableVector[5] / computedDensity;
    double computedStiffeningParameter = conservedVariableVector[6] / computedDensity;

    double velocitySquared = (computedXVelocity * computedXVelocity) + (computedYVelocity * computedYVelocity) + (computedZVelocity * computedZVelocity);
    double computedPressure = EquationOfState::computePressure(computedDensity, (conservedVariableVector[4] / computedDensity) - (velocitySquared / 2), computedAdiabaticIndex,
            computedStiffeningParameter);

    fluxVector[0] = conservedVariableVector[1];
    fluxVector[1] = (conservedVariableVector[1] * computedXVelocity) + computedPressure;
    fluxVector[2] = conservedVariableVector[1] * computedYVelocity;
    fluxVector[3] = conservedVariableVector[1] * computedZVelocity;
    fluxVector[4] = computedXVelocity * (conservedVariableVector[4] + computedPressure);
    fluxVector[5] = conservedVariableVector[1] * computedAdiabaticIndex;
    fluxVector[6] = conservedVariableVector[1] * computedStiffeningParameter;

    return fluxVector;
}

// Computes the flux vector in the x-direction.
vector<double> StateVector::computeFluxVector()
{
    return computeFluxVector(computeConservedVariableVector());
}

// Sets the vector of primitive variables (i.e. density, x/y/z velocities, pressure, adiabatic exponent, and stiffening parameter).
void StateVector::setPrimitiveVariableVector(vector<double> primitiveVariableVector)
{
    density = primitiveVariableVector[0];
    xVelocity = primitiveVariableVector[1];
    yVelocity = primitiveVariableVector[2];
    zVelocity = primitiveVariableVector[3];
    pressure = primitiveVariableVector[4];
    adiabaticIndex = primitiveVariableVector[5];
    stiffeningParameter = primitiveVariableVector[6];
}

// Sets the vector of conserved variables (i.e. density, x/y/z momenta, and conserved forms of the adiabatic exponent and stiffening parameter).
void StateVector::setConservedVariableVector(vector<double> conservedVariableVector)
{
    density = conservedVariableVector[0];
    xVelocity = conservedVariableVector[1] / density;
    yVelocity = conservedVariableVector[2] / density;
    zVelocity = conservedVariableVector[3] / density;
    adiabaticIndex = conservedVariableVector[5] / density;
    stiffeningParameter = conservedVariableVector[6] / density;

    double velocitySquared = (xVelocity * xVelocity) + (yVelocity * yVelocity) + (zVelocity * zVelocity);
    pressure = EquationOfState::computePressure(density, (conservedVariableVector[4] / density) - (velocitySquared / 2), adiabaticIndex, stiffeningParameter);
}

// Sets the density.
void StateVector::setDensity(double newDensity)
{
    density = newDensity;
}

// Sets the x-velocity.
void StateVector::setXVelocity(double newXVelocity)
{
    xVelocity = newXVelocity;
}

// Sets the y-velocity.
void StateVector::setYVelocity(double newYVelocity)
{
    yVelocity = newYVelocity;
}

// Sets the z-velocity.
void StateVector::setZVelocity(double newZVelocity)
{
    zVelocity = newZVelocity;
}

// Sets the pressure.
void StateVector::setPressure(double newPressure)
{
    pressure = newPressure;
}

// Sets the adiabatic exponent.
void StateVector::setAdiabaticIndex(double newAdiabaticIndex)
{
    adiabaticIndex = newAdiabaticIndex;
}

// Sets the stiffening parameter.
void StateVector::setStiffeningParameter(double newStiffeningParameter)
{
    stiffeningParameter = newStiffeningParameter;
}

// Retrieves the density.
double StateVector::getDensity()
{
    return density;
}

// Retrieves the x-velocity.
double StateVector::getXVelocity()
{
    return xVelocity;
}

// Retrieves the y-velocity.
double StateVector::getYVelocity()
{
    return yVelocity;
}

// Retrieves the z-velocity.
double StateVector::getZVelocity()
{
    return zVelocity;
}

// Retrieves the pressure.
double StateVector::getPressure()
{
    return pressure;
}

// Retrieves the adiabatic exponent.
double StateVector::getAdiabaticIndex()
{
    return adiabaticIndex;
}

// Retrieves the stiffening parameter.
double StateVector::getStiffeningParameter()
{
    return stiffeningParameter;
}
