#include "constants.h"
#include <math.h>
#include <fstream>
#include <iostream>

double a_parameter(double r)
{
    return(a_constant[0] + a_constant[1]*r + a_constant[2]*r*r + a_constant[3]*r*r*r + a_constant[4]*r*r*r*r);
}

double b_parameter(double r)
{
    return(b_constant[0] + b_constant[1]*r + b_constant[2]*r*r + b_constant[3]*r*r*r + b_constant[4]*r*r*r*r);
}

double Flux_square_interp(double r, double z)
{
    return(a_parameter(r)*a_parameter(r)/exp(2.5*log(b_parameter(r)*b_parameter(r)+z*z)));
}

double Force_EM_interp(double z, double v)
{
    return((v/0.5)*(2.60062e-16*(0.0022 + z)*(0.0022 + z))/exp(3*log((5.55112e-17 + 0.727808*exp(4*log(0.00236194 + z)) + 1.39979e-6*(0.244827 + z)*(0.244827 + z) + 0.0000626621*fabs(0.00118313 + z) + 25.9432*exp(3*log(fabs(0.00237872 + z)))))));
}

void Get_Force_File(double v_min, double v_max, double dv, double z_min, double z_max, double dz)
{
    std::ofstream ofs ("Force_interp.dat", std::ofstream::out);
    for(double z = z_min; z<z_max; z+=dz)
    {
        ofs<<z<<"\t";
        std::cout<<z<<"\t";
        for(double v=v_min; v<v_max; v+=dv)
        {
            std::cout<<Force_EM_interp(z-0.001,v)<<"\t";
            ofs<<Force_EM_interp(z-0.001,v)<<"\t";
        }
        std::cout<<std::endl;
        ofs<<std::endl;
    }
}
