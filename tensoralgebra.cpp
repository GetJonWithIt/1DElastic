#include "tensoralgebra.h"

// This class encapsulates the routine tensor operations required by nonlinear elastic solvers.
TensorAlgebra::TensorAlgebra()
{
}

// Computes the first invariant of a given tensor (tensor).
double TensorAlgebra::computeFirstInvariant(vector<vector<double> > tensor)
{
    return MatrixAlgebra::computeTrace(tensor);
}

// Computes the second invariant of a given tensor (tensor).
double TensorAlgebra::computeSecondInvariant(vector<vector<double> > tensor)
{
    return 0.5 * ((MatrixAlgebra::computeTrace(tensor) * MatrixAlgebra::computeTrace(tensor)) - MatrixAlgebra::computeComponentSum(MatrixAlgebra::multiplyMatricesElementWise(tensor, tensor)));
}

// Computes the third invariant of a given tensor (tensor).
double TensorAlgebra::computeThirdInvariant(vector<vector<double> > tensor)
{
    return MatrixAlgebra::computeDeterminant(tensor);
}
