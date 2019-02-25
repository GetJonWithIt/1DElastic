#ifndef ELASTICSLOPELIMITERS_H
#define ELASTICSLOPELIMITERS_H

#include "elasticstatevector.h"

class ElasticSlopeLimiters
{
public:
    ElasticSlopeLimiters();

    static double computeR(double steepness, double bias);

    static double computeSuperBeeLimiter(double steepness, double bias);
    static double computeVanLeerLimiter(double steepness, double bias);
    static double computeMinBeeLimiter(double steepness, double bias);

    static double computeSlopeLimiter(double steepness, double bias, int slopeLimiter);

    static vector<double> computeSlopeVector(ElasticStateVector leftStateVector, ElasticStateVector middleStateVector, ElasticStateVector rightStateVector, double bias, int slopeLimiter,
                                      HyperelasticVariables hyperelasticVariables);
};

#endif // ELASTICSLOPELIMITERS_H
