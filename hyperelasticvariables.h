#ifndef HYPERELASTICVARIABLES_H
#define HYPERELASTICVARIABLES_H

#include <cmath>
using namespace std;

class HyperelasticVariables
{
public:
    HyperelasticVariables();
    HyperelasticVariables(double newReferenceMassDensity, double newSoundSpeed, double newShearWaveSpeed, double newSpecificHeatCapacity, double newInitialTemperature,
                          double newAlpha, double newBeta, double newGamma);

    double computeBulkSoundSpeedSquared();
    double computeShearWaveSpeedSquared();

    void setReferenceMassDensity(double newReferenceMassDensity);
    void setSoundSpeed(double newSoundSpeed);
    void setShearWaveSpeed(double newShearWaveSpeed);
    void setSpecificHeatCapacity(double newSpecificHeatCapacity);
    void setInitialTemperature(double newInitialTemperature);

    void setAlpha(double newAlpha);
    void setBeta(double newBeta);
    void setGamma(double newGamma);

    double getReferenceMassDensity();
    double getSoundSpeed();
    double getShearWaveSpeed();
    double getSpecificHeatCapacity();
    double getInitialTemperature();

    double getAlpha();
    double getBeta();
    double getGamma();

private:
    double referenceMassDensity;
    double soundSpeed;
    double shearWaveSpeed;
    double specificHeatCapacity;
    double initialTemperature;

    double alpha;
    double beta;
    double gamma;
};

#endif // HYPERELASTICVARIABLES_H
