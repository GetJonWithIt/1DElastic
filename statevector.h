#ifndef STATEVECTOR_H
#define STATEVECTOR_H

#include "equationofstate.h"
#include <vector>
using namespace std;


class StateVector
{
public:
    StateVector();
    StateVector(double newDensity, double newXVelocity, double newYVelocity, double newZVelocity, double newPressure, double newAdiabaticIndex, double newStiffeningParameter);

    double computeSpecificInternalEnergy();
    double computeTotalEnergy();
    double computeSoundSpeed();

    vector<double> computePrimitiveVariableVector();
    vector<double> computeConservedVariableVector();

    static vector<double> computeFluxVector(vector<double> conservedVariableVector);
    vector<double> computeFluxVector();

    void setPrimitiveVariableVector(vector<double> primitiveVariableVector);
    void setConservedVariableVector(vector<double> conservedVariableVector);

    void setDensity(double newDensity);
    void setXVelocity(double newXVelocity);
    void setYVelocity(double newYVelocity);
    void setZVelocity(double newZVelocity);
    void setPressure(double newPressure);
    void setAdiabaticIndex(double newAdiabaticIndex);
    void setStiffeningParameter(double newStiffeningParameter);

    double getDensity();
    double getXVelocity();
    double getYVelocity();
    double getZVelocity();
    double getPressure();
    double getAdiabaticIndex();
    double getStiffeningParameter();

private:
    double density;
    double xVelocity;
    double yVelocity;
    double zVelocity;
    double pressure;
    double adiabaticIndex;
    double stiffeningParameter;
};

#endif // STATEVECTOR_H
