#include <iostream>
#include "mfd.h"
#include "constants.h"

using namespace std;

int main()
{
    cout << "Hello world!" << endl;
    Rectangular_magnet a;
    double* b = new double[3];
    a.Rectangular_MFD(0,0,0,b);
    cout << b[0]<<"\t"<<b[1]<<"\t"<<b[2] << endl;
    delete b;
    return 0;
}
