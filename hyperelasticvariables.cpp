#include "hyperelasticvariables.h"

// This class encapsulates the material properties used within the Romenski hyperelastic equation of state, as detailed in Titarev et al.

// Instantiates the hyperelastic variables, with the default properties (reference mass density, sound speed, shear wave speed, specific heat capacity, initial temperature, alpha/beta/gamma parameters) of copper.
HyperelasticVariables::HyperelasticVariables()
{
    referenceMassDensity = 8.9;
    soundSpeed = 4.6;
    shearWaveSpeed = 2.1;
    specificHeatCapacity = 0.0004;
    initialTemperature = 300;

    alpha = 1.0;
    beta = 3.0;
    gamma = 2.0;
}


// Instantiates the hyperelastic variables with specified values for the reference mass density, sound speed, shear wave speed, specific heat capacity, initial temperature, and alpha/beta/gamma parameters.
HyperelasticVariables::HyperelasticVariables(double newReferenceMassDensity, double newSoundSpeed, double newShearWaveSpeed, double newSpecificHeatCapacity, double newInitialTemperature,
                                             double newAlpha, double newBeta, double newGamma)
{
    referenceMassDensity = newReferenceMassDensity;
    soundSpeed = newSoundSpeed;
    shearWaveSpeed = newShearWaveSpeed;
    specificHeatCapacity = newSpecificHeatCapacity;
    initialTemperature = newInitialTemperature;

    alpha = newAlpha;
    beta = newBeta;
    gamma = newGamma;
}

// Computes the square of the bulk sound speed within the hyperelastic material.
double HyperelasticVariables::computeBulkSoundSpeedSquared()
{
    return (soundSpeed * soundSpeed) - ((4.0 / 3.0) * (shearWaveSpeed * shearWaveSpeed));
}

// Computes the square of the shear sound speed within the hyperelastic material.
double HyperelasticVariables::computeShearWaveSpeedSquared()
{
    return shearWaveSpeed * shearWaveSpeed;
}

// Sets the reference mass density.
void HyperelasticVariables::setReferenceMassDensity(double newReferenceMassDensity)
{
    referenceMassDensity = newReferenceMassDensity;
}

// Sets the sound speed.
void HyperelasticVariables::setSoundSpeed(double newSoundSpeed)
{
    soundSpeed = newSoundSpeed;
}

// Sets the shear wave speed.
void HyperelasticVariables::setShearWaveSpeed(double newShearWaveSpeed)
{
    shearWaveSpeed = newShearWaveSpeed;
}

// Sets the specific heat capacity.
void HyperelasticVariables::setSpecificHeatCapacity(double newSpecificHeatCapacity)
{
    specificHeatCapacity = newSpecificHeatCapacity;
}

// Sets the initial temperature.
void HyperelasticVariables::setInitialTemperature(double newInitialTemperature)
{
    initialTemperature = newInitialTemperature;
}

// Sets the alpha parameter.
void HyperelasticVariables::setAlpha(double newAlpha)
{
    alpha = newAlpha;
}

// Sets the beta parameter.
void HyperelasticVariables::setBeta(double newBeta)
{
    beta = newBeta;
}

// Sets the gamma parameter.
void HyperelasticVariables::setGamma(double newGamma)
{
    gamma = newGamma;
}

// Retrieves the reference mass density.
double HyperelasticVariables::getReferenceMassDensity()
{
    return referenceMassDensity;
}

// Retrieves the sound speed.
double HyperelasticVariables::getSoundSpeed()
{
    return soundSpeed;
}

// Retrieves the shear wave speed.
double HyperelasticVariables::getShearWaveSpeed()
{
    return shearWaveSpeed;
}

// Retrieves the specific heat capacity.
double HyperelasticVariables::getSpecificHeatCapacity()
{
    return specificHeatCapacity;
}

// Retrieves the initial temperature.
double HyperelasticVariables::getInitialTemperature()
{
    return initialTemperature;
}

// Retrieves the alpha parameter.
double HyperelasticVariables::getAlpha()
{
    return alpha;
}

// Retrieves the beta parameter.
double HyperelasticVariables::getBeta()
{
    return beta;
}

// Retrieves the gamma parameter.
double HyperelasticVariables::getGamma()
{
    return gamma;
}
