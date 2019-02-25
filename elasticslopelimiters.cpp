#include "elasticslopelimiters.h"

// This class encapsulates the methods for slope limiting and boundary value extrapolation, as required by the SLIC scheme for nonlinear elasticity.
ElasticSlopeLimiters::ElasticSlopeLimiters()
{
}

// Computes the value of the xi(R) parameter, as detailed in Toro.
double ElasticSlopeLimiters::computeR(double steepness, double bias)
{
    return 2.0 / ((1.0 - bias) + ((1.0 + bias) * steepness));
}

// Computes the value of the Super-Bee slope limiter.
double ElasticSlopeLimiters::computeSuperBeeLimiter(double steepness, double bias)
{
    if (steepness < 0.0)
    {
        return 0.0;
    }
    else if (steepness < 0.5)
    {
        return 2.0 * steepness;
    }
    else if (steepness < 1)
    {
        return 1.0;
    }
    else
    {
        return min(min(steepness, computeR(steepness, bias)), 2.0);
    }
}

// Computes the value of the Van-Leer slope limiter.
double ElasticSlopeLimiters::computeVanLeerLimiter(double steepness, double bias)
{
    if (steepness < 0.0)
    {
        return 0.0;
    }
    else
    {
        return min((2.0 * steepness) / (1.0 + steepness), computeR(steepness, bias));
    }
}

// Computes the value of the Min-Bee slope limiter.
double ElasticSlopeLimiters::computeMinBeeLimiter(double steepness, double bias)
{
    if (steepness < 0.0)
    {
        return 0.0;
    }
    else if (steepness < 1.0)
    {
        return steepness;
    }
    else
    {
        return min(1.0, computeR(steepness, bias));
    }
}

// Computes the value of the chosen slope limiter (0 = Super-Bee, 1 = Van-Leer, 2 = Min-Bee).
double ElasticSlopeLimiters::computeSlopeLimiter(double steepness, double bias, int slopeLimiter)
{
    if (slopeLimiter == 0)
    {
        return computeSuperBeeLimiter(steepness, bias);
    }
    else if (slopeLimiter == 1)
    {
        return computeVanLeerLimiter(steepness, bias);
    }
    else
    {
        return computeMinBeeLimiter(steepness, bias);
    }
}

// Computes the limited slope vector between three neighbouring cells, with states given by leftStateVector, middleStateVector and rightStateVector, respectively.
vector<double> ElasticSlopeLimiters::computeSlopeVector(ElasticStateVector leftStateVector, ElasticStateVector middleStateVector, ElasticStateVector rightStateVector, double bias,
                                                        int slopeLimiter, HyperelasticVariables hyperelasticVariables)
{
    int componentCount = leftStateVector.computeConservedVariableVector(hyperelasticVariables).size();

    vector<double> leftStateDifference = VectorAlgebra::subtractVectors(middleStateVector.computeConservedVariableVector(hyperelasticVariables),
                                                                        leftStateVector.computeConservedVariableVector(hyperelasticVariables));
    vector<double> rightStateDifference = VectorAlgebra::subtractVectors(rightStateVector.computeConservedVariableVector(hyperelasticVariables),
                                                                         middleStateVector.computeConservedVariableVector(hyperelasticVariables));

    vector<double> slopeVector = VectorAlgebra::addVectors(VectorAlgebra::multiplyVector(0.5 * (1.0 + bias), leftStateDifference),
                                                           VectorAlgebra::multiplyVector(0.5 * (1.0 - bias), rightStateDifference));

    for (int i = 0; i < componentCount; i++)
    {
        double numerator = pow(10, -8);
        double denominator = pow(10, -8);

        if (abs(leftStateDifference[i]) > numerator)
        {
            numerator = leftStateDifference[i];
        }
        if (abs(rightStateDifference[i]) > denominator)
        {
            denominator = rightStateDifference[i];
        }

        double steepness = numerator / denominator;
        slopeVector[i] *= computeSlopeLimiter(steepness, bias, slopeLimiter);
    }

    return slopeVector;
}
