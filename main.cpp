#include <QCoreApplication>

//#include "slicsolver.h"
#include "elasticslicsolver.h"
#include "tests.h"
#include <iostream>
#include <fstream>
#include "unittests.h"
using namespace std;

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    cout << "Please enter which test you would like to run:" << endl;
    cout << "1 - Zhang's First Test (Contact Discontinuity)" << endl;
    cout << "2 - Zhang's Second Test (Five Wave)" << endl;
    cout << "3 - Zhang's Third Test (Seven Wave)" << endl;
    cout << "4 - Zhang's Fourth Test (Solid-Solid Stick Problem)" << endl;
    int test;
    cin >> test;

    cout << "Please enter the desired resolution:" << endl;
    int resolution;
    cin >> resolution;

    cout << "Computing!" << endl;

    if (test == 1)
    {
        Tests::solveZhangTest1(resolution);
    }

    if (test == 2)
    {
        Tests::solveZhangTest2(resolution);
    }

    if (test == 3)
    {
        Tests::solveZhangTest3(resolution);
    }

    if (test == 4)
    {
        Tests::solveZhangTest4(resolution);
    }

    cout << "Computation completed!" << endl;

    return a.exec();
}
