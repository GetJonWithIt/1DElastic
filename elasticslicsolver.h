#ifndef ELASTICSLICSOLVER_H
#define ELASTICSLICSOLVER_H

#include "elasticforcesolver.h"

class ElasticSLICSolver
{
public:
    ElasticSLICSolver();

    static vector<double> computeSLICFlux(ElasticStateVector leftLeftStateVector, ElasticStateVector leftStateVector, ElasticStateVector rightStateVector,
                                          ElasticStateVector rightRightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                          HyperelasticVariables hyperelasticVariables);

    static void computeSLICTimeStep(vector<ElasticStateVector> & newCells, vector<ElasticStateVector> & currentCells, double cellSpacing, double timeStep, double bias,
                               int slopeLimiter, HyperelasticVariables hyperelasticVariables);
    static vector<ElasticStateVector> solve(vector<ElasticStateVector> & cells, double cellSpacing, double CFLCoefficient, double finalTime, double bias, int slopeLimiter,
                                     HyperelasticVariables hyperelasticVariables);
};

#endif // ELASTICSLICSOLVER_H
