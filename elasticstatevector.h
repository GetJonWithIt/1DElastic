#ifndef ELASTICSTATEVECTOR_H
#define ELASTICSTATEVECTOR_H

#include "matrixalgebra.h"
#include "elasticequationofstate.h"
using namespace std;

class ElasticStateVector
{
public:
    ElasticStateVector();
    ElasticStateVector(double newXVelocity, double newYVelocity, double newZVelocity, vector<vector<double> > newDeformationGradient, double newEntropy);

    vector<double> computePrimitiveVariableVector();
    vector<double> computeConservedVariableVector(HyperelasticVariables hyperelasticVariables);

    static vector<double> computeFluxVector(vector<double> conservedVariableVector, HyperelasticVariables hyperelasticVariables);
    vector<double> computeFluxVector(HyperelasticVariables hyperelasticVariables);

    void setPrimitiveVariableVector(vector<double> primitiveVariableVector);
    void setConservedVariableVector(vector<double> conservedVariableVector, HyperelasticVariables hyperelasticVariables);

    double computeDensity(HyperelasticVariables hyperelasticVariables);

    void setXVelocity(double newXVelocity);
    void setYVelocity(double newYVelocity);
    void setZVelocity(double newZVelocity);
    void setDeformationGradient(vector<vector<double> > newDeformationGradient);
    void setEntropy(double newEntropy);

    double getXVelocity();
    double getYVelocity();
    double getZVelocity();
    vector<vector<double> > getDeformationGradient();
    double getEntropy();

private:
    double xVelocity;
    double yVelocity;
    double zVelocity;
    vector<vector<double> > deformationGradient;
    double entropy;
};

#endif // ELASTICSTATEVECTOR_H
