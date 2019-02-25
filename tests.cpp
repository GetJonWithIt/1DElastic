#include "tests.h"

#include "unittests.h"

// This class encapsulates the four relevant nonlinear elastic impact tests from the paper of Zhang et al.
Tests::Tests()
{
}

// Solves Zhang's first test (the contact discontinuity problem) at the given resolution.
void Tests::solveZhangTest1(int cellCount)
{
    vector<ElasticStateVector> initialCells(cellCount);
    HyperelasticVariables hyperelasticVariables(8.93, 4.6, 2.1, 3.9 * pow(10, -4), 300.0, 1.0, 3.0, 2.0);
    double cellSpacing = 1.0 / cellCount;

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

    vector<vector<double> > rightDeformationGradient(3, vector<double>(3));
    rightDeformationGradient[0][0] = 1.0;
    rightDeformationGradient[0][1] = 0.03;
    rightDeformationGradient[0][2] = 0.0;
    rightDeformationGradient[1][0] = 0.02;
    rightDeformationGradient[1][1] = 1.0;
    rightDeformationGradient[1][2] = 0.0;
    rightDeformationGradient[2][0] = 0.0;
    rightDeformationGradient[2][1] = 0.0;
    rightDeformationGradient[2][2] = 1.0;

    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (cellCount / 2))
        {
            //initialCells[i] = ElasticStateVector(0.01, 0.0, 0.0, leftDeformationGradient, pow(10, -3));
            initialCells[i] = ElasticStateVector(0.0, 0.0, 0.0, leftDeformationGradient, pow(10, -3));
        }
        else
        {
            //initialCells[i] = ElasticStateVector(0.01, 0.0, 0.0, rightDeformationGradient, 0.0);
            initialCells[i] = ElasticStateVector(0.0, 0.0, 0.0, rightDeformationGradient, 0.0);
        }
    }

    outputElasticSolution(ElasticSLICSolver::solve(initialCells, cellSpacing, 0.8, 0.7, 0.0, 0, hyperelasticVariables), hyperelasticVariables);
}

// Solves Zhang's second test (the five wave initial condition problem) at the given resolution.
void Tests::solveZhangTest2(int cellCount)
{
    vector<ElasticStateVector> initialCells(cellCount);
    HyperelasticVariables hyperelasticVariables(8.93, 4.6, 2.1, 3.9 * pow(10, -4), 300.0, 1.0, 3.0, 2.0);
    double cellSpacing = 1.0 / cellCount;

    vector<vector<double> > leftDeformationGradient(3, vector<double>(3));
    leftDeformationGradient[0][0] = 0.95;
    leftDeformationGradient[0][1] = 0.0;
    leftDeformationGradient[0][2] = 0.0;
    leftDeformationGradient[1][0] = 0.05;
    leftDeformationGradient[1][1] = 1.0;
    leftDeformationGradient[1][2] = 0.0;
    leftDeformationGradient[2][0] = 0.0;
    leftDeformationGradient[2][1] = 0.0;
    leftDeformationGradient[2][2] = 1.0;

    vector<vector<double> > rightDeformationGradient(3, vector<double>(3));
    rightDeformationGradient[0][0] = 1.0;
    rightDeformationGradient[0][1] = 0.0;
    rightDeformationGradient[0][2] = 0.0;
    rightDeformationGradient[1][0] = 0.0;
    rightDeformationGradient[1][1] = 1.0;
    rightDeformationGradient[1][2] = 0.0;
    rightDeformationGradient[2][0] = 0.0;
    rightDeformationGradient[2][1] = 0.0;
    rightDeformationGradient[2][2] = 1.0;

    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (cellCount / 2))
        {
            initialCells[i] = ElasticStateVector(0.0, 1.0, 0.0, leftDeformationGradient, pow(10, -3));
        }
        else
        {
            initialCells[i] = ElasticStateVector(0.0, 0.0, 0.0, rightDeformationGradient, 0.0);
        }
    }

    outputElasticSolution(ElasticSLICSolver::solve(initialCells, cellSpacing, 0.8, 0.06, 0.0, 0, hyperelasticVariables), hyperelasticVariables);
}

// Solves Zhang's third test (the seven wave initial condition problem) at the given resolution.
void Tests::solveZhangTest3(int cellCount)
{
    vector<ElasticStateVector> initialCells(cellCount);
    HyperelasticVariables hyperelasticVariables(8.93, 4.6, 2.1, 3.9 * pow(10, -4), 300.0, 1.0, 3.0, 2.0);
    double cellSpacing = 1.0 / cellCount;

    vector<vector<double> > leftDeformationGradient(3, vector<double>(3));
    leftDeformationGradient[0][0] = 0.98;
    leftDeformationGradient[0][1] = 0.0;
    leftDeformationGradient[0][2] = 0.0;
    leftDeformationGradient[1][0] = 0.02;
    leftDeformationGradient[1][1] = 1.0;
    leftDeformationGradient[1][2] = 0.1;
    leftDeformationGradient[2][0] = 0.0;
    leftDeformationGradient[2][1] = 0.0;
    leftDeformationGradient[2][2] = 1.0;

    vector<vector<double> > rightDeformationGradient(3, vector<double>(3));
    rightDeformationGradient[0][0] = 1.0;
    rightDeformationGradient[0][1] = 0.0;
    rightDeformationGradient[0][2] = 0.0;
    rightDeformationGradient[1][0] = 0.0;
    rightDeformationGradient[1][1] = 1.0;
    rightDeformationGradient[1][2] = 0.1;
    rightDeformationGradient[2][0] = 0.0;
    rightDeformationGradient[2][1] = 0.0;
    rightDeformationGradient[2][2] = 1.0;

    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (cellCount / 2))
        {
            initialCells[i] = ElasticStateVector(0.0, 0.5, 1.0, leftDeformationGradient, pow(10, -3));
        }
        else
        {
            initialCells[i] = ElasticStateVector(0.0, 0.0, 0.0, rightDeformationGradient, 0.0);
        }
    }

    outputElasticSolution(ElasticSLICSolver::solve(initialCells, cellSpacing, 0.8, 0.06, 0.0, 0, hyperelasticVariables), hyperelasticVariables);
}

// Solves Zhang's fourth test (the solid-solid "stick" problem) at the given resolution.
void Tests::solveZhangTest4(int cellCount)
{
    vector<ElasticStateVector> initialCells(cellCount);
    HyperelasticVariables hyperelasticVariables(8.93, 4.6, 2.1, 3.9 * pow(10, -4), 300.0, 1.0, 3.0, 2.0);
    double cellSpacing = 1.0 / cellCount;

    vector<vector<double> > leftDeformationGradient(3, vector<double>(3));
    leftDeformationGradient[0][0] = 1.0;
    leftDeformationGradient[0][1] = 0.0;
    leftDeformationGradient[0][2] = 0.0;
    leftDeformationGradient[1][0] = -0.01;
    leftDeformationGradient[1][1] = 0.95;
    leftDeformationGradient[1][2] = 0.02;
    leftDeformationGradient[2][0] = -0.015;
    leftDeformationGradient[2][1] = 0.0;
    leftDeformationGradient[2][2] = 0.9;

    vector<vector<double> > rightDeformationGradient(3, vector<double>(3));
    rightDeformationGradient[0][0] = 1.0;
    rightDeformationGradient[0][1] = 0.0;
    rightDeformationGradient[0][2] = 0.0;
    rightDeformationGradient[1][0] = 0.015;
    rightDeformationGradient[1][1] = 0.95;
    rightDeformationGradient[1][2] = 0.0;
    rightDeformationGradient[2][0] = -0.01;
    rightDeformationGradient[2][1] = 0.0;
    rightDeformationGradient[2][2] = 0.9;

    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (cellCount / 2))
        {
            initialCells[i] = ElasticStateVector(2.0, 0.0, 0.1, leftDeformationGradient, 0.0);
        }
        else
        {
            initialCells[i] = ElasticStateVector(0.0, -0.03, -0.01, rightDeformationGradient, 0.0);
        }
    }

    outputElasticSolution(ElasticSLICSolver::solve(initialCells, cellSpacing, 0.6, 0.06, 0.0, 0, hyperelasticVariables), hyperelasticVariables);
}

// Produces output files of density, velocity, specific internal energy, and the first component of the stress tensor, for a given solution to a nonlinear elastic impact test.
void Tests::outputElasticSolution(vector<ElasticStateVector> elasticSolution, HyperelasticVariables hyperelasticVariables)
{
    ofstream densityFile("density.dat");
    ofstream velocityFile("velocity.dat");
    ofstream energyFile("energy.dat");
    ofstream stressTensorFile("stressTensor.dat");
    int cellCount = elasticSolution.size();
    double cellSpacing = 1.0 / cellCount;

    for (int i = 0; i < cellCount; i++)
    {
        densityFile << (cellSpacing * i) << " " << elasticSolution[i].computeDensity(hyperelasticVariables) << endl;
        velocityFile << (cellSpacing * i) << " " << elasticSolution[i].getXVelocity() << endl;
        energyFile << (cellSpacing * i) << " " << ElasticEquationOfState::computeTotalEnergy(MatrixAlgebra::computeInverseMatrix(elasticSolution[i].getDeformationGradient()),
                                                                                             elasticSolution[i].getEntropy(), 0.0, 0.0, 0.0, hyperelasticVariables) << endl;
        stressTensorFile << (cellSpacing * i) << " " << ElasticEquationOfState::computeTotalStressTensor(elasticSolution[i].computeDensity(hyperelasticVariables),
                                                                                                         elasticSolution[i].getDeformationGradient(), elasticSolution[i].getEntropy(),
                                                                                                         hyperelasticVariables)[0][0] << endl;
    }

    densityFile.close();
    velocityFile.close();
    energyFile.close();
    stressTensorFile.close();
}
