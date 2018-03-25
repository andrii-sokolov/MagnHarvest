#include "oscillator.h"
#include "constants.h"
#include <math.h>
#include <string>
#include <fstream>
#include <iostream>

using namespace std;

void harvester::RK45()
{
        double v = 0.0;
        double z = 0.0;
        double z1,z2;
        double p1,p2,p3,p4,l1,l2,l3,l4;
        double ampl_min = z;
        double ampl_max = z;
        double integr = 0.0;
        double integr_ch = 0.0;
        double integr_ne = -1.0;
        double integr_t = 0.0;
        double power = 0.0;
        int counter = 0;
        for(double t=0.0; fabs(integr_ch - integr_ne)>0.001 ;t+=dt)
        {
            z1 = z;
            p1 = Right_part_dimentionless(z,v,t);
            p2 = Right_part_dimentionless(z + p1*0.5*dt, v, t + 0.5*dt);
            p3 = Right_part_dimentionless(z + p2*0.5*dt, v, t + 0.5*dt);
            p4 = Right_part_dimentionless(z + p3*dt, v, t + dt);
            v = v + dt*(p1 + 2.0*p2 + 2.0*p3 + p4)/6.0;
            l1 = v;
            l2 = v + dt*l1*0.5;
            l3 = v + dt*l2*0.5;
            l4 = v + dt*l3;
            z = z + dt*(l1+2.0*l2+2.0*l3+l4)/6.0;
            z2 = z;
            integr+=fabs(z)*dt;
            power+= (Fm*cos(omega*t*sqrt(m/k)) - Fm*beta/sqrt(k*m)*v - Fm*z)*(z2*Fm/k - z1*Fm/k);
            if(z > ampl_max)
                ampl_max = z;
            if(z < ampl_min)
                ampl_min = z;

            if(integr_t > (2*M_PI/omega)/(sqrt(m/k)) )
            {
                //cout<<counter<<"\t"<<integr_ch<<"\t"<<integr_ne<<ampl_max<<"\t"<<ampl_min<<"\t"<<power<<"\n";
                if(counter%2 == 0)
                    integr_ch = integr;
                else
                    integr_ne = integr;
                counter++;
                integr_t = 0.0;
                integr = 0.0;
                harvester_ampl = (ampl_max - ampl_min)/2;
                ampl_min = 0.0;
                ampl_max = 0.0;
                harvester_power = power/(2*M_PI/omega);
                power = 0.0;
            }
            integr_t+=dt;
        }
}

double harvester::Right_part_dimentionless(double khi, double dkhi, double tau)
{
    return(-(beta/sqrt(k*m))*dkhi - khi - Force_EM_interp(Fm/k*khi-z1,Fm/sqrt(k*m)*dkhi)/(Fm)+cos(omega*sqrt(m/k)*tau));
}

double harvester::GetAmpl()
{
    return(harvester_ampl);
}

double harvester::GetPower()
{
    return(harvester_power);
}

void harvester::SetDt(double t)
{
    dt = t;
}

void harvester::SetShaker(double Force_m, double om)
{
    Fm = Force_m;
    omega = om;
}

void harvester::SetOscillator(double spring_k, double mass, double shift, double bet)
{
    k = spring_k;
    m = mass;
    z1 = shift;
    beta = bet;
}

harvester::harvester()
{
    k = 2665.7;
    m = 9.83e-5;
    z1 = 0.000;
    beta = 8.8e-3;

    Fm = m*1.0*9.8;
    omega = 100.0;

    dt = (2*M_PI/(omega))/sqrt(m/k)/100000;

}

double harvester::GetMass()
{
    return(m);
}

double harvester::GetSpring()
{
    return(k);
}

double harvester::GetBeta()
{
    return(beta);
}

void mechanical_resonance(double z1_min, double z1_max, double dz1, double f_min, double f_max, double df)
{
    for(double z1 = z1_min;z1<z1_max+dz1;z1+=dz1)
    {
        std::ofstream ofs ("Mechanical_resonance_shift"+std::to_string(z1)+".dat", std::ofstream::out);
        for(double o = 2*M_PI*f_min; o < 2*M_PI*f_max; o += 2*M_PI*df)
        {
            harvester harv;
            harv.SetShaker(1.0*9.8*harv.GetMass(),o);
            harv.SetOscillator(harv.GetSpring(),harv.GetMass(),z1,harv.GetBeta());
            harv.SetDt((2*M_PI/(o))/sqrt(harv.GetMass()/harv.GetSpring())/100000);
            harv.RK45();
            cout<<o/(2*M_PI)<<"\t"<<harv.GetAmpl()<<"\t"<<harv.GetPower()<<"\n";
            ofs<<o/(2*M_PI)<<"\t"<<harv.GetAmpl()<<"\t"<<harv.GetPower()<<"\n";
        }
        ofs.close();
    }
}
