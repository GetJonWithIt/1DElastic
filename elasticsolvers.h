#ifndef ELASTICSOLVERS_H
#define ELASTICSOLVERS_H

#include "elasticstatevector.h"
#include "elasticslopelimiters.h"
using namespace std;

class ElasticSolvers
{
public:
    ElasticSolvers();

    static vector<ElasticStateVector> insertBoundaryCells(vector<ElasticStateVector> & cells, int boundaryCells);

    static vector<vector<double> > computeFluxJacobian(vector<double> conservedVariableVector, HyperelasticVariables hyperelasticVariables, double epsilon);
    static double computeFluxAcousticWaveSpeed(vector<double> conservedVariableVector, HyperelasticVariables hyperelasticVariables);

    static vector<vector<double> > computeStressTensorCentredDifference(vector<vector<double> > deformationGradient, double entropy, HyperelasticVariables hyperelasticVariables,
                                                                        double epsilon);
    static vector<vector<double> > computeAcousticTensor(vector<vector<double> > deformationGradient, double entropy, HyperelasticVariables hyperelasticVariables, double epsilon);
    static double computeAcousticTensorWaveSpeed(vector<vector<double> > deformationGradient, double entropy, HyperelasticVariables hyperelasticVariables);

    static double computeMaximumWaveSpeed(vector<ElasticStateVector> & cells, HyperelasticVariables hyperelasticVariables);
    static double computeStableTimeStep(vector<ElasticStateVector> & cells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime, int currentIteration,
                                        HyperelasticVariables hyperelasticVariables);

    static ElasticStateVector evolveForHalfTimeStep(ElasticStateVector leftStateVector, ElasticStateVector middleStateVector, ElasticStateVector rightStateVector, double cellSpacing,
                                                    double timeStep, double bias, int slopeLimiter, int side, HyperelasticVariables hyperelasticVariables);
};

#endif // ELASTICSOLVERS_H
