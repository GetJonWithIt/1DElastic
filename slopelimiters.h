#ifndef SLOPELIMITERS_H
#define SLOPELIMITERS_H

#include "statevector.h"
#include "vectoralgebra.h"
#include <vector>
using namespace std;

class SlopeLimiters
{
public:
    SlopeLimiters();

    static double computeR(double steepness, double bias);

    static double computeSuperBeeLimiter(double steepness, double bias);
    static double computeVanLeerLimiter(double steepness, double bias);
    static double computeMinBeeLimiter(double steepness, double bias);

    static double computeSlopeLimiter(double steepness, double bias, int slopeLimiter);

    static vector<double> computeSlopeVector(StateVector leftStateVector, StateVector middleStateVector, StateVector rightStateVector, double bias, int slopeLimiter);
};

#endif // SLOPELIMITERS_H
