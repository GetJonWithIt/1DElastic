#include "elasticequationofstate.h"

// This class encapsulates the hyperelastic equation of state, as detailed in equation (4) of Titarev, Romenski and Toro.

ElasticEquationOfState::ElasticEquationOfState()
{
}

// Computes the first component of the hydrodynamic (i.e. thermal) energy density, using the third invariant of the Finger tensor.

double ElasticEquationOfState::computeHydrodynamicInternalEnergyFirstComponent(double fingerTensorThirdInvariant, double bulkSoundSpeedSquared, double alpha)
{
    return (bulkSoundSpeedSquared / (2.0 * (alpha * alpha))) * (pow(fingerTensorThirdInvariant, (alpha / 2.0)) - 1.0) * (pow(fingerTensorThirdInvariant, (alpha / 2.0)) - 1.0);
}

// Computes the second component of the hydrodynamic (i.e. thermal) energy density, using the third invariant of the Finger tensor.
double ElasticEquationOfState::computeHydrodynamicInternalEnergySecondComponent(double fingerTensorThirdInvariant, double entropy, double specificHeatCapacity,
                                                                                double initialTemperature, double gamma)
{
    return specificHeatCapacity * initialTemperature * pow(fingerTensorThirdInvariant, (gamma / 2.0)) * (exp(entropy / specificHeatCapacity) - 1.0);
}

// Computes the third component of the internal energy (i.e. the internal energy due to shear deformations), using the three invariants of the Finger tensor.
double ElasticEquationOfState::computeShearDeformationInternalEnergy(double fingerTensorFirstInvariant, double fingerTensorSecondInvariant,
                                                                     double fingerTensorThirdInvariant, double shearWaveSpeedSquared, double beta)
{
    return (shearWaveSpeedSquared / 2.0) * pow(fingerTensorThirdInvariant, (beta / 2.0)) * (((fingerTensorFirstInvariant * fingerTensorFirstInvariant) / 3.0) - fingerTensorSecondInvariant);
}

// Computes the product of the Finger tensor and the derivative of the first component of the hydrodynamic (i.e. thermal) energy density, with respect to the Finger tensor.
vector<vector<double> > ElasticEquationOfState::computeHydrodynamicInternalEnergyFirstComponentDerivative(double fingerTensorThirdInvariant, double bulkSoundSpeedSquared, double alpha)
{
    return MatrixAlgebra::multiplyMatrix((bulkSoundSpeedSquared / (2.0 * alpha)) * (pow(fingerTensorThirdInvariant, alpha) - pow(fingerTensorThirdInvariant, (alpha / 2.0))),
                                         MatrixAlgebra::computeIdentityMatrix(3));
}

// Computes the product of the Finger tensor and the derivative of the second component of the hydrodynamic (i.e. thermal) energy density, with respect to the Finger tensor.
vector<vector<double> > ElasticEquationOfState::computeHydrodynamicInternalEnergySecondComponentDerivative(double fingerTensorThirdInvariant, double entropy, double specificHeatCapacity,
                                                                                                           double initialTemperature, double gamma)
{
    return MatrixAlgebra::multiplyMatrix(specificHeatCapacity * initialTemperature * (gamma / 2.0) * (exp(entropy / specificHeatCapacity) - 1.0) *
                                         pow(fingerTensorThirdInvariant, (gamma / 2.0)), MatrixAlgebra::computeIdentityMatrix(3));
}

// Computes the product of the Finger tensor and the derivative of the third component of the internal energy (i.e. the internal energy due to shear deformations), with respect to the
// Finger tensor.
vector<vector<double> > ElasticEquationOfState::computeShearDeformationInternalEnergyDerivative(vector<vector<double> > fingerTensor, double fingerTensorFirstInvariant,
                                                                                                double fingerTensorSecondInvariant, double fingerTensorThirdInvariant,
                                                                                                double shearWaveSpeedSquared, double beta)
{
    return MatrixAlgebra::multiplyMatrix(
                (shearWaveSpeedSquared / 2.0) * pow(fingerTensorThirdInvariant, (beta / 2.0)), MatrixAlgebra::addMatrices(MatrixAlgebra::subtractMatrices(
                    MatrixAlgebra::multiplyMatrix(
                        (beta / 2.0) * (((fingerTensorFirstInvariant * fingerTensorFirstInvariant) / 3.0) - fingerTensorSecondInvariant), MatrixAlgebra::computeIdentityMatrix(3)),
                    MatrixAlgebra::multiplyMatrix((fingerTensorFirstInvariant / 3.0), fingerTensor)), MatrixAlgebra::multiplyMatrices(fingerTensor, fingerTensor)));
}

// Computes the total energy, using the inverse deformation gradient and the entropy.
double ElasticEquationOfState::computeTotalEnergy(vector<vector<double> > inverseDeformationGradient, double entropy, double xVelocity, double yVelocity, double zVelocity,
                                                  HyperelasticVariables hyperelasticVariables)
{
    vector<vector<double> > fingerTensor = MatrixAlgebra::computeGramianMatrix(inverseDeformationGradient);

    double fingerTensorFirstInvariant = TensorAlgebra::computeFirstInvariant(fingerTensor);
    double fingerTensorSecondInvariant = TensorAlgebra::computeSecondInvariant(fingerTensor);
    double fingerTensorThirdInvariant = TensorAlgebra::computeThirdInvariant(fingerTensor);

    double bulkSoundSpeedSquared = hyperelasticVariables.computeBulkSoundSpeedSquared();
    double alpha = hyperelasticVariables.getAlpha();
    double specificHeatCapacity = hyperelasticVariables.getSpecificHeatCapacity();
    double initialTemperature = hyperelasticVariables.getInitialTemperature();
    double gamma = hyperelasticVariables.getGamma();
    double shearWaveSpeedSquared = hyperelasticVariables.computeShearWaveSpeedSquared();
    double beta = hyperelasticVariables.getBeta();

    double internalEnergy = computeHydrodynamicInternalEnergyFirstComponent(fingerTensorThirdInvariant, bulkSoundSpeedSquared, alpha) +
            computeHydrodynamicInternalEnergySecondComponent(fingerTensorThirdInvariant, entropy, specificHeatCapacity, initialTemperature, gamma) +
            computeShearDeformationInternalEnergy(fingerTensorFirstInvariant, fingerTensorSecondInvariant, fingerTensorThirdInvariant, shearWaveSpeedSquared, beta);
    double velocitySquared = (xVelocity * xVelocity) + (yVelocity * yVelocity) + (zVelocity * zVelocity);

    return internalEnergy + (velocitySquared / 2.0);
}

// Computes the entropy, using the deformation gradient and the total energy.
double ElasticEquationOfState::computeEntropy(double totalEnergy, vector<vector<double> > deformationGradient, double xVelocity, double yVelocity, double zVelocity,
                                              HyperelasticVariables hyperelasticVariables)
{
    double bulkSoundSpeedSquared = hyperelasticVariables.computeBulkSoundSpeedSquared();
    double alpha = hyperelasticVariables.getAlpha();
    double shearWaveSpeedSquared = hyperelasticVariables.computeShearWaveSpeedSquared();
    double beta = hyperelasticVariables.getBeta();
    double specificHeatCapacity = hyperelasticVariables.getSpecificHeatCapacity();
    double initialTemperature = hyperelasticVariables.getInitialTemperature();
    double gamma = hyperelasticVariables.getGamma();

    vector<vector<double> > inverseDeformationGradient = MatrixAlgebra::computeInverseMatrix(deformationGradient);
    vector<vector<double> > fingerTensor = MatrixAlgebra::computeGramianMatrix(inverseDeformationGradient);

    double fingerTensorFirstInvariant = TensorAlgebra::computeFirstInvariant(fingerTensor);
    double fingerTensorSecondInvariant = TensorAlgebra::computeSecondInvariant(fingerTensor);
    double fingerTensorThirdInvariant = TensorAlgebra::computeThirdInvariant(fingerTensor);

    double velocitySquared = (xVelocity * xVelocity) + (yVelocity * yVelocity) + (zVelocity * zVelocity);
    double hydrodynamicInternalEnergySecondComponent = totalEnergy - (computeHydrodynamicInternalEnergyFirstComponent(fingerTensorThirdInvariant, bulkSoundSpeedSquared, alpha) +
                                                                      computeShearDeformationInternalEnergy(fingerTensorFirstInvariant, fingerTensorSecondInvariant, fingerTensorThirdInvariant,
                                                                                                              shearWaveSpeedSquared, beta) + (velocitySquared / 2.0));

    return specificHeatCapacity * log(1.0 + (hydrodynamicInternalEnergySecondComponent / (specificHeatCapacity * initialTemperature * pow(fingerTensorThirdInvariant, (gamma / 2.0)))));
}

// Computes the pressure, using the inverse deformation gradient and the entropy.
double ElasticEquationOfState::computePressure(double density, vector<vector<double> > inverseDeformationGradient, double entropy, HyperelasticVariables hyperelasticVariables)
{
    double bulkSoundSpeedSquared = hyperelasticVariables.computeBulkSoundSpeedSquared();
    double alpha = hyperelasticVariables.getAlpha();
    double specificHeatCapacity = hyperelasticVariables.getSpecificHeatCapacity();
    double initialTemperature = hyperelasticVariables.getInitialTemperature();
    double gamma = hyperelasticVariables.getGamma();
    double shearWaveSpeedSquared = hyperelasticVariables.computeShearWaveSpeedSquared();
    double beta = hyperelasticVariables.getBeta();

    vector<vector<double> > fingerTensor = MatrixAlgebra::computeGramianMatrix(inverseDeformationGradient);

    double fingerTensorFirstInvariant = TensorAlgebra::computeFirstInvariant(fingerTensor);
    double fingerTensorSecondInvariant = TensorAlgebra::computeSecondInvariant(fingerTensor);
    double fingerTensorThirdInvariant = TensorAlgebra::computeThirdInvariant(fingerTensor);

    return 2.0 * density * (((bulkSoundSpeedSquared / (2.0 * alpha)) * (pow(fingerTensorThirdInvariant, alpha) - pow(fingerTensorThirdInvariant, (alpha / 2.0)))) +
                          (specificHeatCapacity * initialTemperature * (gamma / 2.0) * (exp(entropy / specificHeatCapacity) - 1.0) * pow(fingerTensorThirdInvariant, (gamma / 2.0))) +
                          (((shearWaveSpeedSquared / 2.0) * pow(fingerTensorThirdInvariant, (beta / 2.0))) * ((beta / 2.0) * (((fingerTensorFirstInvariant * fingerTensorFirstInvariant) / 3.0) -
                                                                                                          fingerTensorSecondInvariant))));
}

// Compute the total stress tensor, using the deformation gradient and the entropy.
vector<vector<double> > ElasticEquationOfState::computeTotalStressTensor(double density, vector<vector<double> > deformationGradient, double entropy,
                                                                         HyperelasticVariables hyperelasticVariables)
{
    double bulkSoundSpeedSquared = hyperelasticVariables.computeBulkSoundSpeedSquared();
    double alpha = hyperelasticVariables.getAlpha();
    double specificHeatCapacity = hyperelasticVariables.getSpecificHeatCapacity();
    double initialTemperature = hyperelasticVariables.getInitialTemperature();
    double gamma = hyperelasticVariables.getGamma();
    double shearWaveSpeedSquared = hyperelasticVariables.computeShearWaveSpeedSquared();
    double beta = hyperelasticVariables.getBeta();

    vector<vector<double> > inverseDeformationGradient = MatrixAlgebra::computeInverseMatrix(deformationGradient);
    vector<vector<double> > fingerTensor = MatrixAlgebra::computeGramianMatrix(inverseDeformationGradient);

    double fingerTensorFirstInvariant = TensorAlgebra::computeFirstInvariant(fingerTensor);
    double fingerTensorSecondInvariant = TensorAlgebra::computeSecondInvariant(fingerTensor);
    double fingerTensorThirdInvariant = TensorAlgebra::computeThirdInvariant(fingerTensor);

    return MatrixAlgebra::multiplyMatrix(
                -2.0 * density, MatrixAlgebra::addMatrices(
                    MatrixAlgebra::addMatrices(
                        computeHydrodynamicInternalEnergyFirstComponentDerivative(fingerTensorThirdInvariant, bulkSoundSpeedSquared, alpha),
                        computeHydrodynamicInternalEnergySecondComponentDerivative(fingerTensorThirdInvariant, entropy, specificHeatCapacity, initialTemperature, gamma)),
                    computeShearDeformationInternalEnergyDerivative(fingerTensor, fingerTensorFirstInvariant, fingerTensorSecondInvariant, fingerTensorThirdInvariant,
                                                                    shearWaveSpeedSquared, beta)));
}
