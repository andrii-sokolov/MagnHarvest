#include <iostream>
#include "mfd.h"
#include "constants.h"
#include "oscillator.h"

using namespace std;

int main()
{
    //mechanical_resonance_debug(-0.003,0.001,0.0001,600,1000,10);
    mechanical_resonance_shift(820,-0.003,0.001,0.0001);
    //Get_Force_File(-2,2,0.1,-0.0002,0.005,0.000001);

    return 0;
}
