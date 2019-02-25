#ifndef TESTS_H
#define TESTS_H

#include <fstream>
#include "elasticslicsolver.h"
using namespace std;

class Tests
{
public:
    Tests();

    static void solveZhangTest1(int cellCount);
    static void solveZhangTest2(int cellCount);
    static void solveZhangTest3(int cellCount);
    static void solveZhangTest4(int cellCount);

    static void outputElasticSolution(vector<ElasticStateVector> elasticSolution, HyperelasticVariables hyperelasticVariables);
};

#endif // TESTS_H
