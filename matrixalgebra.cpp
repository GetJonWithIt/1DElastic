#include "matrixalgebra.h"

// This class encapsulates the routine matrix operations required by nonlinear elastic solvers.
MatrixAlgebra::MatrixAlgebra()
{
}

// Computes the identity matrix of the specified dimension.
vector<vector<double> > MatrixAlgebra::computeIdentityMatrix(int dimension)
{
    vector<vector<double> > identityMatrix(dimension, vector<double>(dimension));

    for (int i = 0; i < dimension; i++)
    {
        for (int j = 0; j < dimension; j++)
        {
            if (i == j)
            {
                identityMatrix[i][j] = 1.0;
            }
            else
            {
                identityMatrix[i][j] = 0.0;
            }
        }
    }

    return identityMatrix;
}

// Adds two matrices (matrix1 and matrix2) together.
vector<vector<double> > MatrixAlgebra::addMatrices(vector<vector<double> > matrix1, vector<vector<double> > matrix2)
{
    int rowCount = matrix1.size();
    int columnCount = matrix1[0].size();
    vector<vector<double> > sumMatrix(rowCount, vector<double>(columnCount));

    for (int i = 0; i < rowCount; i++)
    {
        sumMatrix[i] = VectorAlgebra::addVectors(matrix1[i], matrix2[i]);
    }

    return sumMatrix;
}

// Subtracts one matrix (matrix2) from another (matrix1).
vector<vector<double> > MatrixAlgebra::subtractMatrices(vector<vector<double> > matrix1, vector<vector<double> > matrix2)
{
    return MatrixAlgebra::addMatrices(matrix1, multiplyMatrix(-1.0, matrix2));
}

// Multiplies a matrix (matrix1) by a specified scalar quantity (scalar).
vector<vector<double> > MatrixAlgebra::multiplyMatrix(double scalar, vector<vector<double> > matrix1)
{
    int rowCount = matrix1.size();
    int columnCount = matrix1[0].size();
    vector<vector<double> > productMatrix(rowCount, vector<double>(columnCount));

    for (int i = 0; i < rowCount; i++)
    {
        productMatrix[i] = VectorAlgebra::multiplyVector(scalar, matrix1[i]);
    }

    return productMatrix;
}

// Multiplies two matrices (matrix1 and matrix2) together.
vector<vector<double> > MatrixAlgebra::multiplyMatrices(vector<vector<double> > matrix1, vector<vector<double> > matrix2)
{
    int matrix1RowCount = matrix1.size();
    int matrix1ColumnCount = matrix1[0].size();
    int matrix2ColumnCount = matrix2[0].size();
    vector<vector<double> > productMatrix(matrix1RowCount, vector<double>(matrix2ColumnCount));

    for (int i = 0; i < matrix1RowCount; i++)
    {
        for (int j = 0; j < matrix2ColumnCount; j++)
        {
            for (int k = 0; k < matrix1ColumnCount; k++)
            {
                productMatrix[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }

    return productMatrix;
}

// Multiplies two matrices (matrix1 and matrix2) together, in a componentwise fashion.
vector<vector<double> > MatrixAlgebra::multiplyMatricesElementWise(vector<vector<double> > matrix1, vector<vector<double> > matrix2)
{
    int rowCount = matrix1.size();
    int columnCount = matrix1[0].size();
    vector<vector<double> > productMatrix(rowCount, vector<double>(columnCount));

    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < rowCount; j++)
        {
            productMatrix[i][j] = matrix1[i][j] * matrix2[i][j];
        }
    }

    return productMatrix;
}

// Multiplies a matrix (matrix) and a vector (vector1) together.
vector<double> MatrixAlgebra::multiplyMatrixByVector(vector<vector<double> > matrix, vector<double> vector1)
{
    int componentCount = vector1.size();

    vector<vector<double> > matrix2(1, vector<double>(componentCount));
    matrix2[0] = vector1;
    matrix2 = computeTranspose(matrix2);

    vector<vector<double> > productMatrix = multiplyMatrices(matrix, matrix2);
    return computeTranspose(productMatrix)[0];
}

// Computes the tranpose of a given matrix (matrix).
vector<vector<double> > MatrixAlgebra::computeTranspose(vector<vector<double> > matrix)
{
    int rowCount = matrix.size();
    int columnCount = matrix[0].size();
    vector<vector<double> > transposeMatrix(columnCount, vector<double>(rowCount));

    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            transposeMatrix[j][i] = matrix[i][j];
        }
    }

    return transposeMatrix;
}

// Computes the trace of a given matrix (matrix).
double MatrixAlgebra::computeTrace(vector<vector<double> > matrix)
{
    int rowCount = matrix.size();
    double trace = 0.0;

    for (int i = 0; i < rowCount; i++)
    {
        trace += matrix[i][i];
    }

    return trace;
}

// Computes the determinant of a given 3-by-3 matrix (matrix).
double MatrixAlgebra::computeDeterminant(vector<vector<double> > matrix)
{
    /*
    int rowCount = matrix.size();
    vector<vector<double> > subMatrix(rowCount - 1, vector<double>(rowCount - 1));
    double determinant;

    if (rowCount == 2)
    {
        determinant = ((matrix[0][0] * matrix[1][1]) - (matrix[1][0] * matrix[0][1]));
    }
    else
    {
        for (int k = 0; k < rowCount; k++)
        {
            int subMatrixI = 0;
            for (int i = 1; i < rowCount; i++)
            {
                int subMatrixJ = 0;
                for (int j = 0; j < rowCount; j++)
                {
                    if (j == k)
                    {
                        continue;
                        subMatrix[subMatrixI][subMatrixJ] = matrix[i][j];
                        subMatrixJ += 1;
                    }
                }
                subMatrixI += 1;
            }

            determinant += (pow(-1, k) * matrix[0][k] * computeDeterminant(subMatrix));
        }
    }

    return determinant;
    */

    return abs((matrix[0][0] * ((matrix[1][1] * matrix[2][2]) - (matrix[2][1] * matrix[1][2]))) - (matrix[0][1] * ((matrix[1][0] * matrix[2][2]) - (matrix[1][2] * matrix[2][0]))) +
            (matrix[0][2] * ((matrix[1][0] * matrix[2][1]) - (matrix[1][1] * matrix[2][0]))));
}

// Computes the sum of all components within a given matrix (matrix).
double MatrixAlgebra::computeComponentSum(vector<vector<double> > matrix)
{
    int rowCount = matrix.size();
    double componentSum = 0;

    for (int i = 0; i < rowCount; i++)
    {
        componentSum += VectorAlgebra::computeComponentSum(matrix[i]);
    }

    return componentSum;
}

// Computes the Gramian of a given matrix (matrix).
vector<vector<double> > MatrixAlgebra::computeGramianMatrix(vector<vector<double> > matrix)
{
    return multiplyMatrices(computeTranspose(matrix), matrix);
}

// Computes the inverse of a given matrix (matrix).
vector<vector<double> > MatrixAlgebra::computeInverseMatrix(vector<vector<double> > matrix)
{
    int rowCount = matrix.size();
    vector<vector<double> > inverseMatrix(rowCount, vector<double>(rowCount));
    double determinant = computeDeterminant(matrix);

    if (rowCount == 3)
    {
        inverseMatrix[0][0] = ((matrix[1][1] * matrix[2][2]) - (matrix[2][1] * matrix[1][2])) / determinant;
        inverseMatrix[0][1] = ((matrix[0][2] * matrix[2][1]) - (matrix[0][1] * matrix[2][2])) / determinant;
        inverseMatrix[0][2] = ((matrix[0][1] * matrix[1][2]) - (matrix[0][2] * matrix[1][1])) / determinant;

        inverseMatrix[1][0] = ((matrix[1][2] * matrix[2][0]) - (matrix[1][0] * matrix[2][2])) / determinant;
        inverseMatrix[1][1] = ((matrix[0][0] * matrix[2][2]) - (matrix[0][2] * matrix[2][0])) / determinant;
        inverseMatrix[1][2] = ((matrix[1][0] * matrix[0][2]) - (matrix[0][0] * matrix[1][2])) / determinant;

        inverseMatrix[2][0] = ((matrix[1][0] * matrix[2][1]) - (matrix[2][0] * matrix[1][1])) / determinant;
        inverseMatrix[2][1] = ((matrix[2][0] * matrix[0][1]) - (matrix[0][0] * matrix[2][1])) / determinant;
        inverseMatrix[2][2] = ((matrix[0][0] * matrix[1][1]) - (matrix[1][0] * matrix[0][1])) / determinant;
    }

    return inverseMatrix;
}
