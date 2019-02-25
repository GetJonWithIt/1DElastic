#ifndef ELASTICEQUATIONOFSTATE_H
#define ELASTICEQUATIONOFSTATE_H

#include "tensoralgebra.h"
#include "hyperelasticvariables.h"
using namespace std;

class ElasticEquationOfState
{
public:
    ElasticEquationOfState();

    static double computeHydrodynamicInternalEnergyFirstComponent(double fingerTensorThirdInvariant, double bulkSoundSpeedSquared, double alpha);
    static double computeHydrodynamicInternalEnergySecondComponent(double fingerTensorThirdInvariant, double entropy, double specificHeatCapacity, double initialTemperature, double gamma);
    static double computeShearDeformationInternalEnergy(double fingerTensorFirstInvariant, double fingerTensorSecondInvariant, double fingerTensorThirdInvariant,
                                                        double shearWaveSpeedSquared, double beta);

    static vector<vector<double> > computeHydrodynamicInternalEnergyFirstComponentDerivative(double fingerTensorThirdInvariant, double bulkSoundSpeedSquared, double alpha);
    static vector<vector<double> > computeHydrodynamicInternalEnergySecondComponentDerivative(double fingerTensorThirdInvariant, double entropy, double specificHeatCapacity,
                                                                                              double initialTemperature, double gamma);
    static vector<vector<double> > computeShearDeformationInternalEnergyDerivative(vector<vector<double> > fingerTensor, double fingerTensorFirstInvariant, double fingerTensorSecondInvariant,
                                                                                   double fingerTensorThirdInvariant, double shearWaveSpeedSquared, double beta);

    static double computeTotalEnergy(vector<vector<double> > inverseDeformationGradient, double entropy, double xVelocity, double yVelocity, double zVelocity,
                                     HyperelasticVariables hyperelasticVariables);
    static double computeEntropy(double totalEnergy, vector<vector<double> > deformationGradient, double xVelocity, double yVelocity, double zVelocity,
                                 HyperelasticVariables hyperelasticVariables);
    static double computePressure(double density, vector<vector<double> > inverseDeformationGradient, double entropy, HyperelasticVariables hyperelasticVariables);

    static vector<vector<double> > computeTotalStressTensor(double density, vector<vector<double> > deformationGradient, double entropy, HyperelasticVariables hyperelasticVariables);

};

#endif // ELASTICEQUATIONOFSTATE_H
