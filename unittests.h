#ifndef UNITTESTS_H
#define UNITTESTS_H

#include "elasticstatevector.h"
#include <iostream>
using namespace std;

class UnitTests
{
public:
    UnitTests();

    static vector<vector<double> > computeAnalyticHydrodynamicInternalEnergyFirstComponentDerivative(ElasticStateVector stateVector, HyperelasticVariables hyperelasticVariables);
    static vector<vector<double> > computeAnalyticHydrodynamicInternalEnergySecondComponentDerivative(ElasticStateVector stateVector, HyperelasticVariables hyperelasticVariables);
    static vector<vector<double> > computeAnalyticShearDeformationInternalEnergyDerivative(ElasticStateVector stateVector, HyperelasticVariables hyperelasticVariables);

    static vector<vector<double> > computeNumericalHydrodynamicInternalEnergyFirstComponentDerivative(ElasticStateVector stateVector, HyperelasticVariables hyperelasticVariables,
                                                                                                      double epsilon);
    static vector<vector<double> > computeNumericalHydrodynamicInternalEnergySecondComponentDerivative(ElasticStateVector stateVector, HyperelasticVariables hyperelasticVariables,
                                                                                                       double epsilon);
    static vector<vector<double> > computeNumericalShearDeformationInternalEnergyDerivative(ElasticStateVector stateVector, HyperelasticVariables hyperelasticVariables,
                                                                                            double epsilon);

    static void computeZhangTest1LeftMatrices();

    static void outputMatrix(vector<vector<double> > matrix);
};

#endif // UNITTESTS_H
