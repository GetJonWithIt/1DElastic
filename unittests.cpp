#include "unittests.h"

// This class encapulates the unit tests for the stress tensor computation (i.e. derivatives of the equation of state parameters).
UnitTests::UnitTests()
{
}

// Computes the product of the Finger tensor and the analytic derivative of the first component of the hydrodynamic (i.e. thermal) energy density, with respect to the Finger tensor.
vector<vector<double> > UnitTests::computeAnalyticHydrodynamicInternalEnergyFirstComponentDerivative(ElasticStateVector stateVector, HyperelasticVariables hyperelasticVariables)
{
    vector<vector<double> > deformationGradient = stateVector.getDeformationGradient();
    vector<vector<double> > inverseDeformationGradient = MatrixAlgebra::computeInverseMatrix(deformationGradient);
    vector<vector<double> > fingerTensor = MatrixAlgebra::computeGramianMatrix(inverseDeformationGradient);

    double fingerTensorThirdInvariant = TensorAlgebra::computeThirdInvariant(fingerTensor);
    double bulkSoundSpeedSquared = hyperelasticVariables.computeBulkSoundSpeedSquared();
    double alpha = hyperelasticVariables.getAlpha();

    return ElasticEquationOfState::computeHydrodynamicInternalEnergyFirstComponentDerivative(fingerTensorThirdInvariant, bulkSoundSpeedSquared, alpha);
}

// Computes the product of the Finger tensor and the analytic derivative of the second component of the hydrodynamic (i.e. thermal) energy density, with respect to the Finger tensor.
vector<vector<double> > UnitTests::computeAnalyticHydrodynamicInternalEnergySecondComponentDerivative(ElasticStateVector stateVector, HyperelasticVariables hyperelasticVariables)
{
    vector<vector<double> > deformationGradient = stateVector.getDeformationGradient();
    vector<vector<double> > inverseDeformationGradient = MatrixAlgebra::computeInverseMatrix(deformationGradient);
    vector<vector<double> > fingerTensor = MatrixAlgebra::computeGramianMatrix(inverseDeformationGradient);

    double fingerTensorThirdInvariant = TensorAlgebra::computeThirdInvariant(fingerTensor);
    double entropy = stateVector.getEntropy();
    double specificHeatCapacity = hyperelasticVariables.getSpecificHeatCapacity();
    double initialTemperature = hyperelasticVariables.getInitialTemperature();
    double gamma = hyperelasticVariables.getGamma();

    return ElasticEquationOfState::computeHydrodynamicInternalEnergySecondComponentDerivative(fingerTensorThirdInvariant, entropy, specificHeatCapacity, initialTemperature, gamma);
}

// Computes the product of the Finger tensor and the analytic derivative of the third compont of the internal energy (i.e. the internal energy due to shear deformations), with respect
// to the Finger tensor.
vector<vector<double> > UnitTests::computeAnalyticShearDeformationInternalEnergyDerivative(ElasticStateVector stateVector, HyperelasticVariables hyperelasticVariables)
{
    vector<vector<double> > deformationGradient = stateVector.getDeformationGradient();
    vector<vector<double> > inverseDeformationGradient = MatrixAlgebra::computeInverseMatrix(deformationGradient);
    vector<vector<double> > fingerTensor = MatrixAlgebra::computeGramianMatrix(inverseDeformationGradient);

    double fingerTensorFirstInvariant = TensorAlgebra::computeFirstInvariant(fingerTensor);
    double fingerTensorSecondInvariant = TensorAlgebra::computeSecondInvariant(fingerTensor);
    double fingerTensorThirdInvariant = TensorAlgebra::computeThirdInvariant(fingerTensor);
    double shearWaveSpeedSquared = hyperelasticVariables.computeShearWaveSpeedSquared();
    double beta = hyperelasticVariables.getBeta();

    return ElasticEquationOfState::computeShearDeformationInternalEnergyDerivative(fingerTensor, fingerTensorFirstInvariant, fingerTensorSecondInvariant, fingerTensorThirdInvariant,
                                                                                   shearWaveSpeedSquared, beta);
}

// Computes the product of the Finger tensor and the numerical (second-order, centred) derivative of the first component of the hydrodynamic (i.e. thermal) energy density, with respect
// to the Finger tensor.
vector<vector<double> > UnitTests::computeNumericalHydrodynamicInternalEnergyFirstComponentDerivative(ElasticStateVector stateVector, HyperelasticVariables hyperelasticVariables,
                                                                                                      double epsilon)
{
    vector<vector<double> > hydrodynamicInternalEnergyFirstComponentDerivative(3, vector<double>(3));

    vector<vector<double> > deformationGradient = stateVector.getDeformationGradient();
    vector<vector<double> > inverseDeformationGradient = MatrixAlgebra::computeInverseMatrix(deformationGradient);
    vector<vector<double> > fingerTensor = MatrixAlgebra::computeGramianMatrix(inverseDeformationGradient);

    double bulkSoundSpeedSquared = hyperelasticVariables.computeBulkSoundSpeedSquared();
    double alpha = hyperelasticVariables.getAlpha();

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            vector<vector<double> > fingerTensorForwardStep = fingerTensor;
            vector<vector<double> > fingerTensorBackwardStep = fingerTensor;

            fingerTensorForwardStep[i][j] += epsilon;
            fingerTensorBackwardStep[i][j] -= epsilon;

            double fingerTensorThirdInvariantForwardStep = TensorAlgebra::computeThirdInvariant(fingerTensorForwardStep);
            double fingerTensorThirdInvariantBackwardStep = TensorAlgebra::computeThirdInvariant(fingerTensorBackwardStep);

            hydrodynamicInternalEnergyFirstComponentDerivative[i][j] =
                    (ElasticEquationOfState::computeHydrodynamicInternalEnergyFirstComponent(fingerTensorThirdInvariantForwardStep, bulkSoundSpeedSquared, alpha) -
                     ElasticEquationOfState::computeHydrodynamicInternalEnergyFirstComponent(fingerTensorThirdInvariantBackwardStep, bulkSoundSpeedSquared, alpha)) / (2 * epsilon);
        }
    }

    return MatrixAlgebra::multiplyMatrices(fingerTensor, hydrodynamicInternalEnergyFirstComponentDerivative);
}

// Computes the product of the Finger tensor and the numerical (second-order, centred) derivative of the second component of the hydrodynamic (i.e. thermal) energy density, with
// respect to the Finger tensor.
vector<vector<double> > UnitTests::computeNumericalHydrodynamicInternalEnergySecondComponentDerivative(ElasticStateVector stateVector, HyperelasticVariables hyperelasticVariables,
                                                                                                       double epsilon)
{
    vector<vector<double> > hydrodynamicInternalEnergySecondComponentDerivative(3, vector<double>(3));

    vector<vector<double> > deformationGradient = stateVector.getDeformationGradient();
    vector<vector<double> > inverseDeformationGradient = MatrixAlgebra::computeInverseMatrix(deformationGradient);
    vector<vector<double> > fingerTensor = MatrixAlgebra::computeGramianMatrix(inverseDeformationGradient);

    double entropy = stateVector.getEntropy();
    double specificHeatCapacity = hyperelasticVariables.getSpecificHeatCapacity();
    double initialTemperature = hyperelasticVariables.getInitialTemperature();
    double gamma = hyperelasticVariables.getGamma();

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            vector<vector<double> > fingerTensorForwardStep = fingerTensor;
            vector<vector<double> > fingerTensorBackwardStep = fingerTensor;

            fingerTensorForwardStep[i][j] += epsilon;
            fingerTensorBackwardStep[i][j] -= epsilon;

            double fingerTensorThirdInvariantForwardStep = TensorAlgebra::computeThirdInvariant(fingerTensorForwardStep);
            double fingerTensorThirdInvariantBackwardStep = TensorAlgebra::computeThirdInvariant(fingerTensorBackwardStep);

            hydrodynamicInternalEnergySecondComponentDerivative[i][j] =
                    (ElasticEquationOfState::computeHydrodynamicInternalEnergySecondComponent(fingerTensorThirdInvariantForwardStep, entropy, specificHeatCapacity,
                                                                                              initialTemperature, gamma) -
                     ElasticEquationOfState::computeHydrodynamicInternalEnergySecondComponent(fingerTensorThirdInvariantBackwardStep, entropy, specificHeatCapacity,
                                                                                              initialTemperature, gamma)) / (2 * epsilon);
        }
    }

    return MatrixAlgebra::multiplyMatrices(fingerTensor, hydrodynamicInternalEnergySecondComponentDerivative);
}

// Computes the product of the Finger tensor and the numerical (second-order, centred) derivative of the third component of the internal energy (i.e. the internal energy due to shear
// deformations), with respect to the Finger tensor.
vector<vector<double> > UnitTests::computeNumericalShearDeformationInternalEnergyDerivative(ElasticStateVector stateVector, HyperelasticVariables hyperelasticVariables, double epsilon)
{
    vector<vector<double> > shearDeformationInternalEnergyDerivative(3, vector<double>(3));

    vector<vector<double> > deformationGradient = stateVector.getDeformationGradient();
    vector<vector<double> > inverseDeformationGradient = MatrixAlgebra::computeInverseMatrix(deformationGradient);
    vector<vector<double> > fingerTensor = MatrixAlgebra::computeGramianMatrix(inverseDeformationGradient);

    double shearWaveSpeedSquared = hyperelasticVariables.computeShearWaveSpeedSquared();
    double beta = hyperelasticVariables.getBeta();

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            vector<vector<double> > fingerTensorForwardStep = fingerTensor;
            vector<vector<double> > fingerTensorBackwardStep = fingerTensor;

            fingerTensorForwardStep[i][j] += epsilon;
            fingerTensorBackwardStep[i][j] -= epsilon;

            double fingerTensorFirstInvariantForwardStep = TensorAlgebra::computeFirstInvariant(fingerTensorForwardStep);
            double fingerTensorFirstInvariantBackwardStep = TensorAlgebra::computeFirstInvariant(fingerTensorBackwardStep);
            double fingerTensorSecondInvariantForwardStep = TensorAlgebra::computeSecondInvariant(fingerTensorForwardStep);
            double fingerTensorSecondInvariantBackwardStep = TensorAlgebra::computeSecondInvariant(fingerTensorBackwardStep);
            double fingerTensorThirdInvariantForwardStep = TensorAlgebra::computeThirdInvariant(fingerTensorForwardStep);
            double fingerTensorThirdInvariantBackwardStep = TensorAlgebra::computeThirdInvariant(fingerTensorBackwardStep);

            shearDeformationInternalEnergyDerivative[i][j] =
                    (ElasticEquationOfState::computeShearDeformationInternalEnergy(fingerTensorFirstInvariantForwardStep, fingerTensorSecondInvariantForwardStep,
                                                                                   fingerTensorThirdInvariantForwardStep, shearWaveSpeedSquared, beta) -
                     ElasticEquationOfState::computeShearDeformationInternalEnergy(fingerTensorFirstInvariantBackwardStep, fingerTensorSecondInvariantBackwardStep,
                                                                                   fingerTensorThirdInvariantBackwardStep, shearWaveSpeedSquared, beta)) / (2 * epsilon);
        }
    }

    return MatrixAlgebra::multiplyMatrices(fingerTensor, shearDeformationInternalEnergyDerivative);
}

// Computes the three matrices for the left-hand-side of Zhang's first (contact discontinuity) test, using both analytical and numerical methods.
void UnitTests::computeZhangTest1LeftMatrices()
{
    HyperelasticVariables hyperelasticVariables = HyperelasticVariables(8.93, 4.6, 2.1, 3.9 * pow(10, -4), 300.0, 1.0, 3.0, 2.0);

    vector<vector<double> > leftDeformationGradient(3, vector<double>(3));
    leftDeformationGradient[0][0] = 1.156276139;
    leftDeformationGradient[0][1] = 0.034688284;
    leftDeformationGradient[0][2] = 0.0;
    leftDeformationGradient[1][0] = 0.093190648;
    leftDeformationGradient[1][1] = 1.002195719;
    leftDeformationGradient[1][2] = 0.0;
    leftDeformationGradient[2][0] = 0.0;
    leftDeformationGradient[2][1] = 0.0;
    leftDeformationGradient[2][2] = 1.0;

    vector<vector<double> > analyticHydrodynamicInternalEnergyFirstComponentDerivative =
            computeAnalyticHydrodynamicInternalEnergyFirstComponentDerivative(ElasticStateVector(0.01, 0.0, 0.0, leftDeformationGradient, pow(10, -3)), hyperelasticVariables);
    vector<vector<double> > analyticHydrodynamicInternalEnergySecondComponentDerivative =
            computeAnalyticHydrodynamicInternalEnergySecondComponentDerivative(ElasticStateVector(0.01, 0.0, 0.0, leftDeformationGradient, pow(10, -3)), hyperelasticVariables);
    vector<vector<double> > analyticShearDeformationInternalEnergyDerivative =
            computeAnalyticShearDeformationInternalEnergyDerivative(ElasticStateVector(0.01, 0.0, 0.0, leftDeformationGradient, pow(10, -3)), hyperelasticVariables);

    vector<vector<double> > numericalHydrodynamicInternalEnergyFirstComponentDerivative =
            computeNumericalHydrodynamicInternalEnergyFirstComponentDerivative(ElasticStateVector(0.01, 0.0, 0.0, leftDeformationGradient, pow(10, -3)), hyperelasticVariables, pow(10, -8));
    vector<vector<double> > numericalHydrodynamicInternalEnergySecondComponentDerivative =
            computeNumericalHydrodynamicInternalEnergySecondComponentDerivative(ElasticStateVector(0.01, 0.0, 0.0, leftDeformationGradient, pow(10, -3)), hyperelasticVariables, pow(10, -8));
    vector<vector<double> > numericalShearDeformationInternalEnergyDerivative =
            computeNumericalShearDeformationInternalEnergyDerivative(ElasticStateVector(0.01, 0.0, 0.0, leftDeformationGradient, pow(10, -3)), hyperelasticVariables, pow(10, -8));

    cout << "Analytically-computed matrices:" << endl << endl;

    outputMatrix(analyticHydrodynamicInternalEnergyFirstComponentDerivative);
    outputMatrix(analyticHydrodynamicInternalEnergySecondComponentDerivative);
    outputMatrix(analyticShearDeformationInternalEnergyDerivative);

    cout << "Numerically-computed matrices:" << endl << endl;

    outputMatrix(numericalHydrodynamicInternalEnergyFirstComponentDerivative);
    outputMatrix(numericalHydrodynamicInternalEnergySecondComponentDerivative);
    outputMatrix(numericalShearDeformationInternalEnergyDerivative);
}

void UnitTests::outputMatrix(vector<vector<double> > matrix)
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}
