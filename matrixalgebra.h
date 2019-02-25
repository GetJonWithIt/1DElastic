#ifndef MATRIXALGEBRA_H
#define MATRIXALGEBRA_H

#include "vectoralgebra.h"
#include <cmath>
using namespace std;

class MatrixAlgebra
{
public:
    MatrixAlgebra();

    static vector<vector<double> > computeIdentityMatrix(int dimension);

    static vector<vector<double> > addMatrices(vector<vector<double> > matrix1, vector<vector<double> > matrix2);
    static vector<vector<double> > subtractMatrices(vector<vector<double> > matrix1, vector<vector<double> > matrix2);

    static vector<vector<double> > multiplyMatrix(double scalar, vector<vector<double> > matrix1);
    static vector<vector<double> > multiplyMatrices(vector<vector<double> > matrix1, vector<vector<double> > matrix2);
    static vector<vector<double> > multiplyMatricesElementWise(vector<vector<double> > matrix1, vector<vector<double> > matrix2);

    static vector<double> multiplyMatrixByVector(vector<vector<double> > matrix, vector<double> vector1);

    static vector<vector<double> > computeTranspose(vector<vector<double> > matrix);

    static double computeTrace(vector<vector<double> > matrix);
    static double computeDeterminant(vector<vector<double> > matrix);

    static double computeComponentSum(vector<vector<double> > matrix);

    static vector<vector<double> > computeGramianMatrix(vector<vector<double> > matrix);

    static vector<vector<double> > computeInverseMatrix(vector<vector<double> > matrix);
};

#endif // MATRIXALGEBRA_H
