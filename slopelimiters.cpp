#include "slopelimiters.h"

// This class encapsulated the methods for slope limiting and boundary value extrapolation, as required by the SLIC scheme for the (stiffened gas) Euler equations.
SlopeLimiters::SlopeLimiters()
{
}

// Computes the value of xi(R) parameter, as detailed in Toro.
double SlopeLimiters::computeR(double steepness, double bias)
{
    return 2 / ((1 - bias) + ((1 + bias) * steepness));
}

// Computes the value of the Super-Bee slope limiter.
double SlopeLimiters::computeSuperBeeLimiter(double steepness, double bias)
{
    if (steepness < 0)
    {
        return 0;
    }
    else if (steepness < 0.5)
    {
        return 2 * steepness;
    }
    else if (steepness < 1)
    {
        return 1;
    }
    else
    {
        return min(min(steepness, computeR(steepness, bias)), 2.0);
    }
}

// Computes the value of the Van-Leer slope limiter.
double SlopeLimiters::computeVanLeerLimiter(double steepness, double bias)
{
    if (steepness < 0)
    {
        return 0;
    }
    else
    {
        return min((2 * steepness) / (1 + steepness), computeR(steepness, bias));
    }
}

// Computes the value of the Min-Bee slope limiter.
double SlopeLimiters::computeMinBeeLimiter(double steepness, double bias)
{
    if (steepness < 0)
    {
        return 0;
    }
    else if (steepness < 1)
    {
        return steepness;
    }
    else
    {
        return min(1.0, computeR(steepness, bias));
    }
}

// Computes the value of the chosen slope limiter (0 = Super-Bee, 1 = Van-Leer, 2 = Min-Bee).
double SlopeLimiters::computeSlopeLimiter(double steepness, double bias, int slopeLimiter)
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
vector<double> SlopeLimiters::computeSlopeVector(StateVector leftStateVector, StateVector middleStateVector, StateVector rightStateVector, double bias, int slopeLimiter)
{
    int componentCount = leftStateVector.computeConservedVariableVector().size();
    vector<double> leftStateDifference = VectorAlgebra::subtractVectors(middleStateVector.computeConservedVariableVector(), leftStateVector.computeConservedVariableVector());
    vector<double> rightStateDifference = VectorAlgebra::subtractVectors(rightStateVector.computeConservedVariableVector(), middleStateVector.computeConservedVariableVector());

    vector<double> slopeVector = VectorAlgebra::addVectors(VectorAlgebra::multiplyVector(0.5 * (1 + bias), leftStateDifference),
                                                           VectorAlgebra::multiplyVector(0.5 * (1 - bias), rightStateDifference));

    for (int i = 0; i < componentCount; i++)
    {
        double numerator = pow(10, -5);
        double denominator = pow(10, -5);

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
