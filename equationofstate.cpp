#include "equationofstate.h"

EquationOfState::EquationOfState()
{
}

// Computes the specific internal energy for a stiffened gas.
double EquationOfState::computeSpecificInternalEnergy(double density, double pressure, double adiabaticIndex, double stiffeningParameter)
{
    return (pressure + (adiabaticIndex * stiffeningParameter)) / ((adiabaticIndex - 1) * density);
}

// Computes the pressure for a stiffened gas.
double EquationOfState::computePressure(double density, double specificInternalEnergy, double adiabaticIndex, double stiffeningParameter)
{
    return (specificInternalEnergy * (adiabaticIndex - 1) * density) - (adiabaticIndex * stiffeningParameter);
}

// Computes the sound speed for a stiffened gas.
double EquationOfState::computeSoundSpeed(double density, double pressure, double adiabaticIndex, double stiffeningParameter)
{
    return sqrt((adiabaticIndex * (pressure + stiffeningParameter)) / density);
}

// Computes the entropy of an ideal gas, using its density, pressure and adiabatic exponent.
double EquationOfState::computeEntropy(double density, double pressure, double adiabaticIndex)
{
    return pressure / pow(density, adiabaticIndex);
}
