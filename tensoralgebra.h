#ifndef TENSORALGEBRA_H
#define TENSORALGEBRA_H

#include "matrixalgebra.h"
using namespace std;

class TensorAlgebra
{
public:
    TensorAlgebra();

    static double computeFirstInvariant(vector<vector<double> > tensor);
    static double computeSecondInvariant(vector<vector<double> > tensor);
    static double computeThirdInvariant(vector<vector<double> > tensor);
};

#endif // TENSORALGEBRA_H
