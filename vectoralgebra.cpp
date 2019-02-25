#include "vectoralgebra.h"

// This class encapsulates the routine vector operations required by nonlinear elastic solvers.
VectorAlgebra::VectorAlgebra()
{
}

// Adds two vectors (vector1 and vector2) together.
vector<double> VectorAlgebra::addVectors(vector<double> vector1, vector<double> vector2)
{
    int componentCount = vector1.size();
    vector<double> sumVector(componentCount);

    for (int i = 0; i < componentCount; i++)
    {
        sumVector[i] = vector1[i] + vector2[i];
    }

    return sumVector;
}

// Subtracts one vector (vector2) from another (vector1).
vector<double> VectorAlgebra::subtractVectors(vector<double> vector1, vector<double> vector2)
{
    return addVectors(vector1, multiplyVector(-1.0, vector2));
}

// Multiplies a vector (vector1) by a specified scalar quantity (scalar).
vector<double> VectorAlgebra::multiplyVector(double scalar, vector<double> vector1)
{
    int componentCount = vector1.size();
    vector<double> productVector(componentCount);

    for (int i = 0; i < componentCount; i++)
    {
        productVector[i] = scalar * vector1[i];
    }

    return productVector;
}

// Computes the dot product of two vectors (vector1 and vector2).
double VectorAlgebra::computeDotProduct(vector<double> vector1, vector<double> vector2)
{
    int componentCount = vector1.size();
    double dotProduct = 0;

    for (int i = 0; i < componentCount; i++)
    {
        dotProduct += vector1[i] * vector2[i];
    }

    return dotProduct;
}

// Computes the norm of a given vector (vector1).
double VectorAlgebra::computeNorm(vector<double> vector1)
{
    return sqrt(computeDotProduct(vector1, vector1));
}

// Computes the sum of all components within a given vector (vector1).
double VectorAlgebra::computeComponentSum(vector<double> vector1)
{
    int componentCount = vector1.size();
    double componentSum = 0.0;

    for (int i = 0; i < componentCount; i++)
    {
        componentSum += vector1[i];
    }

    return componentSum;
}
