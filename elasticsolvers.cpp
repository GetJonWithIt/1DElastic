#include "elasticsolvers.h"

// This class encapsulates the methods that are common to all solvers for nonlinear elastic systems.
ElasticSolvers::ElasticSolvers()
{
}

// Returns the computational domain, but with a given number of boundary cells inserted on each side (transmissive boundary conditions).
vector<ElasticStateVector> ElasticSolvers::insertBoundaryCells(vector<ElasticStateVector> & cells, int boundaryCells)
{
    int cellCount = cells.size();
    vector<ElasticStateVector> cellsWithBoundary(cellCount + (2 * boundaryCells));

    if (boundaryCells == 1)
    {
        cellsWithBoundary[0] = cells[0];
        cellsWithBoundary[cellCount + 1] = cells[cellCount - 1];
    }
    else if (boundaryCells == 2)
    {
        cellsWithBoundary[0] = cells[1];
        cellsWithBoundary[1] = cells[0];
        cellsWithBoundary[cellCount + 2] = cells[cellCount - 1];
        cellsWithBoundary[cellCount + 3] = cells[cellCount - 2];
    }

    for (int i = 0; i < cellCount; i++)
    {
        cellsWithBoundary[i + boundaryCells] = cells[i];
    }

    return cellsWithBoundary;
}

// DEPRECATED! Computes the Jacobian of the flux vector, with respect to the conserved variable vector.
vector<vector<double> > ElasticSolvers::computeFluxJacobian(vector<double> conservedVariableVector, HyperelasticVariables hyperelasticVariables, double epsilon)
{
    vector<vector<double> > fluxJacobian(13, vector<double>(13));

    for (int i = 0; i < 13; i++)
    {
        for (int j = 0; j < 13; j++)
        {
            vector<double> conservedVariableVectorForwardStep = conservedVariableVector;
            vector<double> conservedVariableVectorBackwardStep = conservedVariableVector;

            conservedVariableVectorForwardStep[j] += epsilon;
            conservedVariableVectorBackwardStep[j] -= epsilon;

            vector<double> fluxVectorForwardStep = ElasticStateVector::computeFluxVector(conservedVariableVectorForwardStep, hyperelasticVariables);
            vector<double> fluxVectorBackwardStep = ElasticStateVector::computeFluxVector(conservedVariableVectorBackwardStep, hyperelasticVariables);

            fluxJacobian[i][j] = (fluxVectorForwardStep[i] - fluxVectorBackwardStep[i]) / (2.0 * epsilon);
        }
    }

    return fluxJacobian;
}

// DEPRECATED! Computes the largest absolute eigenvalue of the flux Jacobian using power iteration, for the purposes of wavespeed calculation.
double ElasticSolvers::computeFluxAcousticWaveSpeed(vector<double> conservedVariableVector, HyperelasticVariables hyperelasticVariables)
{
    vector<vector<double> > fluxJacobian = computeFluxJacobian(conservedVariableVector, hyperelasticVariables, pow(10, -8));

    vector<double> eigenvectorGuess(13);
    for(int i = 0; i < 13; i++)
    {
        eigenvectorGuess[i] = 1;
    }

    for (int i = 0; i < 20; i++)
    {
        eigenvectorGuess = MatrixAlgebra::multiplyMatrixByVector(fluxJacobian, eigenvectorGuess);
        eigenvectorGuess = VectorAlgebra::multiplyVector((1.0 / VectorAlgebra::computeNorm(eigenvectorGuess)), eigenvectorGuess);
    }

    return VectorAlgebra::computeNorm(MatrixAlgebra::multiplyMatrixByVector(fluxJacobian, eigenvectorGuess));
}

// Computes the centred difference of the stress tensor with respect to the inverse deformation gradient, for the purposes of acoustic tensor calculation.
vector<vector<double> > ElasticSolvers::computeStressTensorCentredDifference(vector<vector<double> > deformationGradient, double entropy, HyperelasticVariables hyperelasticVariables,
                                                                         double epsilon)
{
    vector<vector<double> > stressTensorCentredDifference(3, vector<double>(3));
    vector<vector<double> > inverseDeformationGradient = MatrixAlgebra::computeInverseMatrix(deformationGradient);

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            vector<vector<double> > inverseDeformationGradientForwardStep = inverseDeformationGradient;
            vector<vector<double> > inverseDeformationGradientBackwardStep = inverseDeformationGradient;

            inverseDeformationGradientForwardStep[0][j] += epsilon;
            inverseDeformationGradientBackwardStep[0][j] -= epsilon;

            vector<vector<double> > deformationGradientForwardStep = MatrixAlgebra::computeInverseMatrix(inverseDeformationGradientForwardStep);
            vector<vector<double> > deformationGradientBackwardStep = MatrixAlgebra::computeInverseMatrix(inverseDeformationGradientBackwardStep);

            double densityForwardStep;
            double densityBackwardStep;
            if (MatrixAlgebra::computeDeterminant(deformationGradientForwardStep) <= pow(10, -8))
            {
                densityForwardStep = hyperelasticVariables.getReferenceMassDensity() / pow(10, -8);
            }
            else
            {
                densityForwardStep = hyperelasticVariables.getReferenceMassDensity() / MatrixAlgebra::computeDeterminant(deformationGradientForwardStep);
            }
            if (MatrixAlgebra::computeDeterminant(deformationGradientBackwardStep) <= pow(10, -8))
            {
                densityBackwardStep = hyperelasticVariables.getReferenceMassDensity() / pow(10, -8);
            }
            else
            {
                densityBackwardStep = hyperelasticVariables.getReferenceMassDensity() / MatrixAlgebra::computeDeterminant(deformationGradientBackwardStep);
            }

            vector<vector<double> > stressTensorForwardStep = ElasticEquationOfState::computeTotalStressTensor(densityForwardStep, deformationGradientForwardStep, entropy,
                                                                                                               hyperelasticVariables);
            vector<vector<double> > stressTensorBackwardStep = ElasticEquationOfState::computeTotalStressTensor(densityBackwardStep, deformationGradientBackwardStep, entropy,
                                                                                                                hyperelasticVariables);

            stressTensorCentredDifference[i][j] = stressTensorForwardStep[0][i] - stressTensorBackwardStep[0][i];
        }
    }

    return stressTensorCentredDifference;
}

// Computes the approximated acoustic tensor, using an eighth-order approximation of the derivative of the stress tensor with respect to the deformation gradient.
vector<vector<double> > ElasticSolvers::computeAcousticTensor(vector<vector<double> > deformationGradient, double entropy, HyperelasticVariables hyperelasticVariables, double epsilon)
{
    vector<vector<double> > stressTensorDerivativeFirstStep = MatrixAlgebra::multiplyMatrix(
                (4.0 / (5.0 * epsilon)), computeStressTensorCentredDifference(deformationGradient, entropy, hyperelasticVariables, epsilon));
    vector<vector<double> > stressTensorDerivativeSecondStep = MatrixAlgebra::multiplyMatrix(
                (1.0 / (5.0  * epsilon)), computeStressTensorCentredDifference(deformationGradient, entropy, hyperelasticVariables, 2 * epsilon));
    vector<vector<double> > stressTensorDerivativeThirdStep = MatrixAlgebra::multiplyMatrix(
                (4.0 / (105.0 * epsilon)), computeStressTensorCentredDifference(deformationGradient, entropy, hyperelasticVariables, 3 * epsilon));
    vector<vector<double> > stressTensorDerivativeFourthStep = MatrixAlgebra::multiplyMatrix(
                (1.0 / (280.0 * epsilon)), computeStressTensorCentredDifference(deformationGradient, entropy, hyperelasticVariables, 4 * epsilon));

    vector<vector<double> > stressTensorDerivative = MatrixAlgebra::subtractMatrices(
                MatrixAlgebra::addMatrices(
                    MatrixAlgebra::subtractMatrices(stressTensorDerivativeFirstStep, stressTensorDerivativeSecondStep), stressTensorDerivativeThirdStep), stressTensorDerivativeFourthStep);

    return MatrixAlgebra::multiplyMatrix(-1.0, MatrixAlgebra::multiplyMatrices(stressTensorDerivative, deformationGradient));
}

// Computes the estimated maximum acoustic wave speed, using an explicit eigenvalue decomposition of the acoustic tensor.
double ElasticSolvers::computeAcousticTensorWaveSpeed(vector<vector<double> > deformationGradient, double entropy, HyperelasticVariables hyperelasticVariables)
{
    vector<vector<double> > acousticTensor = computeAcousticTensor(deformationGradient, entropy, hyperelasticVariables, pow(10, -8));

    double A = -MatrixAlgebra::computeTrace(acousticTensor);
    double B = 0.5 * ((MatrixAlgebra::computeTrace(acousticTensor) * MatrixAlgebra::computeTrace(acousticTensor)) -
                      MatrixAlgebra::computeTrace(MatrixAlgebra::multiplyMatrices(acousticTensor, acousticTensor)));
    double C = -MatrixAlgebra::computeDeterminant(acousticTensor);

    double quotient = ((3.0 * B) - (A * A)) / 9.0;
    double remainder = ((9.0 * A * B) - (27.0 * C) - (2.0 * A * A * A)) / 54.0;
    double discriminant = (quotient * quotient * quotient) + (remainder * remainder);

    double theta = acos(remainder / sqrt(-(quotient * quotient * quotient)));

    double pi = atan(1.0) * 4.0;
    double eigenvalue1 = (2.0 * sqrt(-quotient) * cos(theta / 3.0)) - (A / 3.0);
    double eigenvalue2 = (2.0 * sqrt(-quotient) * cos((theta / 3.0) + ((2.0 * pi) / 3.0))) - (A / 3.0);
    double eigenvalue3 =  (2.0 * sqrt(-quotient) * cos((theta / 3.0) + ((4.0 * pi) / 3.0))) - (A / 3.0);

    double density = 0.0;
    if (MatrixAlgebra::computeDeterminant(deformationGradient) < pow(10, -8))
    {
        density = hyperelasticVariables.getReferenceMassDensity() / pow(10, -8);
    }
    else
    {
        density = hyperelasticVariables.getReferenceMassDensity() / MatrixAlgebra::computeDeterminant(deformationGradient);
    }

    return sqrt(max(abs(eigenvalue1), max(abs(eigenvalue2), abs(eigenvalue3))) / density);
}

// Computes the estimated maximum wave speed across the entire computational domain, for the purposes of timestep calculation.
double ElasticSolvers::computeMaximumWaveSpeed(vector<ElasticStateVector> & cells, HyperelasticVariables hyperelasticVariables)
{
    double maximumWaveSpeed = 0;

    for (int i = 0; i < cells.size(); i++)
    {
        vector<vector<double> > deformationGradient = cells[i].getDeformationGradient();
        double entropy = cells[i].getEntropy();
        double waveSpeed = abs(cells[i].getXVelocity()) + computeAcousticTensorWaveSpeed(deformationGradient, entropy, hyperelasticVariables);

        if (waveSpeed > maximumWaveSpeed)
        {
            maximumWaveSpeed = waveSpeed;
        }
    }

    return maximumWaveSpeed;
}

// Computes the maximum stable timestep, in accordance with the CFL condition.
double ElasticSolvers::computeStableTimeStep(vector<ElasticStateVector> & cells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime, int currentIteration,
                                             HyperelasticVariables hyperelasticVariables)
{
    double timeStep = CFLCoefficient * (cellSpacing / computeMaximumWaveSpeed(cells, hyperelasticVariables));

    if (currentIteration <= 5)
    {
        timeStep *= 0.2;
    }

    if (currentTime + timeStep > finalTime)
    {
        timeStep = finalTime - currentTime;
    }

    return timeStep;
}

// Evolves the state given by middleStateVector for half a timestep, using the slope across the neighbouring cells, given by leftStateVector, middleStateVector and rightStateVector,
// respectively.
ElasticStateVector ElasticSolvers::evolveForHalfTimeStep(ElasticStateVector leftStateVector, ElasticStateVector middleStateVector, ElasticStateVector rightStateVector, double cellSpacing,
                                                         double timeStep, double bias, int slopeLimiter, int side, HyperelasticVariables hyperelasticVariables)
{
    vector<double> slopeVector = ElasticSlopeLimiters::computeSlopeVector(leftStateVector, middleStateVector, rightStateVector, bias, slopeLimiter, hyperelasticVariables);
    vector<double> leftExtrapolatedValue = VectorAlgebra::subtractVectors(middleStateVector.computeConservedVariableVector(hyperelasticVariables),
                                                                          VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> rightExtrapolatedValue = VectorAlgebra::addVectors(middleStateVector.computeConservedVariableVector(hyperelasticVariables),
                                                                      VectorAlgebra::multiplyVector(0.5, slopeVector));

    vector<double> fluxEvolution = VectorAlgebra::multiplyVector(0.5 * (timeStep / cellSpacing), VectorAlgebra::subtractVectors(
                                                                     ElasticStateVector::computeFluxVector(leftExtrapolatedValue, hyperelasticVariables),
                                                                     ElasticStateVector::computeFluxVector(rightExtrapolatedValue, hyperelasticVariables)));

    ElasticStateVector evolvedStateVector;
    if (side == 0)
    {
        evolvedStateVector.setConservedVariableVector(VectorAlgebra::addVectors(leftExtrapolatedValue, fluxEvolution), hyperelasticVariables);
    }
    else
    {
        evolvedStateVector.setConservedVariableVector(VectorAlgebra::addVectors(rightExtrapolatedValue, fluxEvolution), hyperelasticVariables);
    }

    return evolvedStateVector;
}
