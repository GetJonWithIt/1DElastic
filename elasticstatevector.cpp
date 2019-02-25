#include "elasticstatevector.h"

// This class encapsulates the complete state vector for a nonlinear elastic system, as detailed in equations 5 - 23 of Miller and Colella.

// Constructs a state vector with default values for the x/y/z velocities, deformation gradient, and entropy.
ElasticStateVector::ElasticStateVector()
{
    xVelocity = 0.0;
    yVelocity = 0.0;
    zVelocity = 0.0;
    deformationGradient = vector<vector<double> >(3, vector<double>(3));
    entropy = 0.0;
}

// Constructs a state vector with specified values for the x/y/z velocities, deformation gradient and entropy.
ElasticStateVector::ElasticStateVector(double newXVelocity, double newYVelocity, double newZVelocity, vector<vector<double> > newDeformationGradient, double newEntropy)
{
    xVelocity = newXVelocity;
    yVelocity = newYVelocity;
    zVelocity = newZVelocity;
    deformationGradient = newDeformationGradient;
    entropy =  newEntropy;
}

// Computes the vector of primitive variables (i.e. x/y/z velocities, specific deformation gradient, and entropy).
vector<double> ElasticStateVector::computePrimitiveVariableVector()
{
    vector<double> primitiveVariableVector(13);

    primitiveVariableVector[0] = xVelocity;
    primitiveVariableVector[1] = yVelocity;
    primitiveVariableVector[2] = zVelocity;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            primitiveVariableVector[3 + (i * 3) + j] = deformationGradient[i][j];
        }
    }
    primitiveVariableVector[12] = entropy;

    return primitiveVariableVector;
}

// Computes the vector of "conserved" variables (i.e. x/y/z momenta, total deformation gradient, and total energy).
vector<double> ElasticStateVector::computeConservedVariableVector(HyperelasticVariables hyperelasticVariables)
{
    vector<double> conservedVariableVector(13);

    double density;
    if (MatrixAlgebra::computeDeterminant(deformationGradient) <= pow(10, -8))
    {
        density = hyperelasticVariables.getReferenceMassDensity() / pow(10, -8);
    }
    else
    {
        density = hyperelasticVariables.getReferenceMassDensity() / MatrixAlgebra::computeDeterminant(deformationGradient);
    }
    vector<vector<double> > inverseDeformationGradient = MatrixAlgebra::computeInverseMatrix(deformationGradient);
    double totalEnergy = ElasticEquationOfState::computeTotalEnergy(inverseDeformationGradient, entropy, xVelocity, yVelocity, zVelocity, hyperelasticVariables);

    conservedVariableVector[0] = density * xVelocity;
    conservedVariableVector[1] = density * yVelocity;
    conservedVariableVector[2] = density * zVelocity;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            conservedVariableVector[3 + (i * 3) + j] = density * deformationGradient[i][j];
        }
    }
    conservedVariableVector[12] = density * totalEnergy;

    return conservedVariableVector;
}

// Computes the flux vector in the x-direction, using the vector of "conserved" variables.
vector<double> ElasticStateVector::computeFluxVector(vector<double> conservedVariableVector, HyperelasticVariables hyperelasticVariables)
{
    vector<double> fluxVector(13);

    double conservedXVelocity = conservedVariableVector[0];
    double conservedYVelocity = conservedVariableVector[1];
    double conservedZVelocity = conservedVariableVector[2];
    vector<vector<double> > conservedDeformationGradient(3, vector<double>(3));
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            conservedDeformationGradient[i][j] = conservedVariableVector[3 + (i * 3) + j];
        }
    }
    double conservedEnergy = conservedVariableVector[12];

    double referenceMassDensity = hyperelasticVariables.getReferenceMassDensity();
    double density = sqrt(MatrixAlgebra::computeDeterminant(conservedDeformationGradient) / referenceMassDensity);

    if (abs(density) <= pow(10, -8))
    {
        density = pow(10, -8);
    }

    double computedXVelocity = conservedXVelocity / density;
    double computedYVelocity = conservedYVelocity / density;
    double computedZVelocity = conservedZVelocity / density;
    vector<vector<double> > computedDeformationGradient(3, vector<double>(3));
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            computedDeformationGradient[i][j] = conservedDeformationGradient[i][j] / density;
        }
    }
    double computedEnergy = conservedEnergy / density;

    double computedEntropy = ElasticEquationOfState::computeEntropy(computedEnergy, computedDeformationGradient, computedXVelocity, computedYVelocity, computedZVelocity, hyperelasticVariables);
    vector<vector<double> > computedTotalStressTensor = ElasticEquationOfState::computeTotalStressTensor(density, computedDeformationGradient, computedEntropy, hyperelasticVariables);

    fluxVector[0] = (computedXVelocity * conservedXVelocity) - computedTotalStressTensor[0][0];
    fluxVector[1] = (computedXVelocity * conservedYVelocity) - computedTotalStressTensor[0][1];
    fluxVector[2] = (computedXVelocity * conservedZVelocity) - computedTotalStressTensor[0][2];
    for (int i = 0; i < 3; i++)
    {
        fluxVector[3 + i] = (computedXVelocity * conservedDeformationGradient[0][i]) - (computedXVelocity * conservedDeformationGradient[0][i]);
    }
    for (int i = 0; i < 3; i++)
    {
        fluxVector[6 + i] = (computedXVelocity * conservedDeformationGradient[1][i]) - (computedYVelocity * conservedDeformationGradient[0][i]);
    }
    for (int i = 0; i < 3; i++)
    {
        fluxVector[9 + i] = (computedXVelocity * conservedDeformationGradient[2][i]) - (computedZVelocity * conservedDeformationGradient[0][i]);
    }

    vector<double> velocityVector(3);
    velocityVector[0] = computedXVelocity;
    velocityVector[1] = computedYVelocity;
    velocityVector[2] = computedZVelocity;
    fluxVector[12] = (computedXVelocity * conservedEnergy) - VectorAlgebra::computeDotProduct(computedTotalStressTensor[0], velocityVector);

    return fluxVector;
}

// Computes the flux vector in the x-direction.
vector<double> ElasticStateVector::computeFluxVector(HyperelasticVariables hyperelasticVariables)
{
    return computeFluxVector(computeConservedVariableVector(hyperelasticVariables), hyperelasticVariables);
}

// Sets the vector of primitive variables (i.e. x/y/z velocities, specific deformation gradient, and entropy).
void ElasticStateVector::setPrimitiveVariableVector(vector<double> primitiveVariableVector)
{
    xVelocity = primitiveVariableVector[0];
    yVelocity = primitiveVariableVector[1];
    zVelocity = primitiveVariableVector[2];
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            deformationGradient[i][j] = primitiveVariableVector[3 + (i * 3) + j];
        }
    }
    entropy = primitiveVariableVector[12];
}

// Sets the vector of "conserved" variables (i.e. x/y/z momenta, total deformation gradient, and total energy).
void ElasticStateVector::setConservedVariableVector(vector<double> conservedVariableVector, HyperelasticVariables hyperelasticVariables)
{
    double conservedXVelocity = conservedVariableVector[0];
    double conservedYVelocity = conservedVariableVector[1];
    double conservedZVelocity = conservedVariableVector[2];
    vector<vector<double> > conservedDeformationGradient(3, vector<double>(3));
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            conservedDeformationGradient[i][j] = conservedVariableVector[3 + (i * 3) + j];
        }
    }
    double conservedEnergy = conservedVariableVector[12];

    double referenceMassDensity = hyperelasticVariables.getReferenceMassDensity();
    double density = sqrt(MatrixAlgebra::computeDeterminant(conservedDeformationGradient) / referenceMassDensity);

    if (abs(density) <= pow(10, -8))
    {
        density = pow(10, -8);
    }

    xVelocity = conservedXVelocity / density;
    yVelocity = conservedYVelocity / density;
    zVelocity = conservedZVelocity / density;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            deformationGradient[i][j] = conservedDeformationGradient[i][j] / density;
        }
    }
    double computedEnergy = conservedEnergy / density;

    entropy = ElasticEquationOfState::computeEntropy(computedEnergy, deformationGradient, xVelocity, yVelocity, zVelocity, hyperelasticVariables);
}

// Computes the mass density.
double ElasticStateVector::computeDensity(HyperelasticVariables hyperelasticVariables)
{
    if (MatrixAlgebra::computeDeterminant(deformationGradient) <= pow(10, -8))
    {
        return hyperelasticVariables.getReferenceMassDensity() / pow(10, -8);
    }
    else
    {
        return hyperelasticVariables.getReferenceMassDensity() / MatrixAlgebra::computeDeterminant(deformationGradient);
    }
}

// Sets the x-velocity.
void ElasticStateVector::setXVelocity(double newXVelocity)
{
    xVelocity = newXVelocity;
}

// Sets the y-velocity.
void ElasticStateVector::setYVelocity(double newYVelocity)
{
    yVelocity = newYVelocity;
}

// Sets the z-velocity.
void ElasticStateVector::setZVelocity(double newZVelocity)
{
    zVelocity = newZVelocity;
}

// Sets the deformation gradient.
void ElasticStateVector::setDeformationGradient(vector<vector<double> > newDeformationGradient)
{
    deformationGradient = newDeformationGradient;
}

// Sets the entropy.
void ElasticStateVector::setEntropy(double newEntropy)
{
    entropy = newEntropy;
}

// Retrieves the x-velocity.
double ElasticStateVector::getXVelocity()
{
    return xVelocity;
}

// Retrieves the y-velocity.
double ElasticStateVector::getYVelocity()
{
    return yVelocity;
}

// Retrieves the z-velocity.
double ElasticStateVector::getZVelocity()
{
    return zVelocity;
}

// Retrieves the deformation gradient.
vector<vector<double> > ElasticStateVector::getDeformationGradient()
{
    return deformationGradient;
}

// Retrieves the entropy.
double ElasticStateVector::getEntropy()
{
    return entropy;
}
