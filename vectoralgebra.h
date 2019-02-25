#ifndef VECTORALGEBRA_H
#define VECTORALGEBRA_H

#include <vector>
#include <cmath>
using namespace std;

class VectorAlgebra
{
public:
    VectorAlgebra();

    static vector<double> addVectors(vector<double> vector1, vector<double> vector2);
    static vector<double> subtractVectors(vector<double> vector1, vector<double> vector2);

    static vector<double> multiplyVector(double scalar, vector<double> vector1);

    static double computeDotProduct(vector<double> vector1, vector<double> vector2);
    static double computeNorm(vector<double> vector1);

    static double computeComponentSum(vector<double> vector1);

};

#endif // VECTORALGEBRA_H
