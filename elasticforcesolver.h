#ifndef ELASTICFORCESOLVER_H
#define ELASTICFORCESOLVER_H

#include "elasticsolvers.h"
using namespace std;

class ElasticFORCESolver
{
public:
    ElasticFORCESolver();

    static vector<double> computeLaxFriedrichsFlux(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, double cellSpacing, double timeStep,
                                                   HyperelasticVariables hyperelasticVariables);
    static vector<double> computeRichtmyerFlux(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, double cellSpacing, double timeStep,
                                               HyperelasticVariables hyperelasticVariables);

    static vector<double> computeFORCEFlux(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, double cellSpacing, double timeStep,
                                          HyperelasticVariables hyperelasticVariables);

    static void computeFORCETimeStep(vector<ElasticStateVector> & newCells, vector<ElasticStateVector> & currentCells, double cellSpacing, double timeStep,
                                     HyperelasticVariables hyperelasticVariables);
    static vector<ElasticStateVector> solve(vector<ElasticStateVector> & cells, double cellSpacing, double CFLCoefficient, double finalTime, HyperelasticVariables hyperelasticVariables);
};

#endif // ELASTICFORCESOLVER_H
