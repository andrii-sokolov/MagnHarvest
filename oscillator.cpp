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
        for(double t=0.0; fabs(integr_ch - integr_ne)>0.0001 ;t+=dt)
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

void harvester::RK45_ndl()
{
        dt = (2*M_PI/(omega))/1000000;
        double v = 0.0;
        double z = 0.0;
        double z_1,z_2;
        double p1,p2,p3,p4,l1,l2,l3,l4;
        double ampl_min = z;
        double ampl_max = z;
        double integr = 0.0;
        double integr_ch = 0.0;
        double integr_ne = -1.0;
        double integr_t = 0.0;
        double power = 0.0;
        double power_el = 0.0;
        int counter = 0;
        double EMF;
        for(double t=0.0; fabs(integr_ch-integr_ne)>1e-14 ;t+=dt)
        {
            z_1 = z;
            p1 = Right_part(z,v,t);
            p2 = Right_part(z + p1*0.5*dt, v, t + 0.5*dt);
            p3 = Right_part(z + p2*0.5*dt, v, t + 0.5*dt);
            p4 = Right_part(z + p3*dt, v, t + dt);
            v = v + dt*(p1 + 2.0*p2 + 2.0*p3 + p4)/6.0;
            l1 = v;
            l2 = v + dt*l1*0.5;
            l3 = v + dt*l2*0.5;
            l4 = v + dt*l3;
            z = z + dt*(l1+2.0*l2+2.0*l3+l4)/6.0;
            z_2 = z;
            integr+=fabs(z)*dt;
            power+= Force_EM_interp(z,v)*(z_2 - z_1);
            EMF = 0.0;
            for(double a = coil_a_min; a<=coil_a_max; a+=(coil_a_max-coil_a_min)/coil_n_of_loops)
            {
                EMF+=(Flux_square_interp(a,z_2-z1)-Flux_square_interp(a,z_1-z1))/dt;
            }
            power_el+=EMF*EMF/coil_resistance*dt;
            if(z > ampl_max)
                ampl_max = z;
            if(z < ampl_min)
                ampl_min = z;

            if(integr_t > (2*M_PI/omega) )
            {
                //cout<<counter<<"\t"<<integr_ch<<"\t"<<integr_ne<<ampl_max<<"\t"<<ampl_min<<"\t"<<power<<"\n";
                if(counter%2 == 0)
                    integr_ch = integr;
                else
                    integr_ne = integr;
                counter++;
                cout<<counter<<"\n";
                integr_t = 0.0;
                integr = 0.0;
                harvester_ampl = (ampl_max - ampl_min)/2;
                ampl_min = 0.0;
                ampl_max = 0.0;
                harvester_power = power/(2*M_PI/omega);
                harvester_power_el = power_el/(2*M_PI/omega);
                power = 0.0;
                power_el = 0.0;
            }
            integr_t+=dt;


        }
}

void harvester::RK45_debug()
{
        dt = (2*M_PI/(omega))/1000000;
        double v = 0.0;
        double z = 0.0;
        double z_1,z_2;
        double p1,p2,p3,p4,l1,l2,l3,l4;
        double ampl_min = z;
        double ampl_max = z;
        double integr = 0.0;
        double integr_ch = 0.0;
        double integr_ne = -1.0;
        double integr_t = 0.0;
        double power_sp = 0.0;
        double power_du = 0.0;
        double power_ex = 0.0;
        double power = 0.0;
        int counter = 0;
        for(double t=0.0; fabs(integr_ch-integr_ne)>1e-14 ;t+=dt)
        {
            z_1 = z;
            p1 = Right_part(z,v,t);
            p2 = Right_part(z + p1*0.5*dt, v, t + 0.5*dt);
            p3 = Right_part(z + p2*0.5*dt, v, t + 0.5*dt);
            p4 = Right_part(z + p3*dt, v, t + dt);
            v = v + dt*(p1 + 2.0*p2 + 2.0*p3 + p4)/6.0;
            l1 = v;
            l2 = v + dt*l1*0.5;
            l3 = v + dt*l2*0.5;
            l4 = v + dt*l3;
            z = z + dt*(l1+2.0*l2+2.0*l3+l4)/6.0;
            z_2 = z;
            integr+=fabs(z)*dt;
            power_du+= -beta*v*(z_2 - z_1);
            power_sp+= -k*z*(z_2 - z_1);
            power_ex+= Fm*cos(omega*t)*(z_2 - z_1);
            power+=Force_EM_interp(z-z1,v)*(z_2 - z_1);
            if(z > ampl_max)
                ampl_max = z;
            if(z < ampl_min)
                ampl_min = z;

            if(integr_t > (2*M_PI/omega) )
            {
                //cout<<counter<<"\t"<<integr_ch<<"\t"<<integr_ne<<ampl_max<<"\t"<<ampl_min<<"\t"<<power<<"\n";
                if(counter%2 == 0)
                    integr_ch = integr;
                else
                    integr_ne = integr;
                counter++;
                //cout<<counter<<"\n";
                integr_t = 0.0;
                integr = 0.0;
                harvester_ampl = (ampl_max - ampl_min)/2;
                ampl_min = 0.0;
                ampl_max = 0.0;
                harvester_power_dump = power_du/(2*M_PI/omega);
                harvester_power_spring = power_sp/(2*M_PI/omega);
                harvester_power_ext = power_ex/(2*M_PI/omega);
                harvester_power = power/(2*M_PI/omega);
                power_du = 0.0;
                power_ex = 0.0;
                power_sp = 0.0;
                power = 0.0;
            }
            integr_t+=dt;


        }
}

double harvester::Right_part_dimentionless(double khi, double dkhi, double tau)
{
    return(-(beta/sqrt(k*m))*dkhi - khi - Force_EM_interp(Fm/k*khi-z1,Fm/sqrt(k*m)*dkhi)/(Fm)+cos(omega*sqrt(m/k)*tau));
}

double harvester::Right_part(double x, double v, double t)
{
    return((-beta*v - k*x  +Fm*cos(omega*t)- Force_EM_interp(x-z1,v))/m);
}

double harvester::Right_part_debug(double x, double v, double t)
{
    return((-beta*v - k*x  +Fm*cos(omega*t))/m);
}

double harvester::GetAmpl()
{
    return(harvester_ampl);
}

double harvester::GetPower()
{
    return(harvester_power);
}

double harvester::GetPowerEl()
{
    return(harvester_power_el);
}

double harvester::GetPowerSpring()
{
    return(harvester_power_spring);
}

double harvester::GetPowerExt()
{
    return(harvester_power_ext);
}

double harvester::GetPowerDump()
{
    return(harvester_power_dump);
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
    k1 = 0.0;
    m = 9.83e-5;
    z1 = 0.000;
    beta = 0.01536;

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

harvester::harvester(double z_shift, double om)
{
    k = 2665.7;
    k1 = 0.0;
    m = 9.83e-5;
    z1 = z_shift;
    beta = 0.01536;

    Fm = m*1.0*9.8;
    omega = om;

    dt = (2*M_PI/(omega))/100000;
    harvester::RK45_debug();
}

void mechanical_resonance(double z1_min, double z1_max, double dz1, double f_min, double f_max, double df)
{
    for(double z1 = z1_min;z1<z1_max+dz1;z1+=dz1)
    {
        std::ofstream ofs ("Mechanical_resonance_shift"+std::to_string(z1)+".dat", std::ofstream::out);
        for(double o = 2*M_PI*f_min; o < 2*M_PI*f_max; o += 2*M_PI*df)
        {
            harvester harv;
            harv.SetShaker(10.0*9.8*harv.GetMass(),o);
            harv.SetOscillator(harv.GetSpring(),harv.GetMass(),z1,harv.GetBeta());
            harv.SetDt((2*M_PI/(o))/sqrt(harv.GetMass()/harv.GetSpring())/100000);
            harv.RK45();
            cout<<o/(2*M_PI)<<"\t"<<harv.GetAmpl()<<"\t"<<harv.GetPower()<<"\n";
            ofs<<o/(2*M_PI)<<"\t"<<harv.GetAmpl()<<"\t"<<harv.GetPower()<<"\n";
        }
        ofs.close();
    }
}

void mechanical_resonance_nd(double z1_min, double z1_max, double dz1, double f_min, double f_max, double df)
{
    for(double z1 = z1_min;z1<z1_max+dz1;z1+=dz1)
    {
        std::ofstream ofs ("Mechanical_resonance_shift"+std::to_string(z1)+".dat", std::ofstream::out);
        for(double o = 2*M_PI*f_min; o < 2*M_PI*f_max; o += 2*M_PI*df)
        {
            harvester harv(z1,o);
            cout<<o/(2*M_PI)<<"\t"<<harv.GetAmpl()<<"\t"<<harv.GetPower()<<"\n";
            ofs<<o/(2*M_PI)<<"\t"<<harv.GetAmpl()<<"\t"<<harv.GetPower()<<"\n";
        }
        ofs.close();
    }
}

void mechanical_resonance_debug(double z1_min, double z1_max, double dz1, double f_min, double f_max, double df)
{
    for(double z1 = z1_min;z1<z1_max+dz1;z1+=dz1)
    {
        std::ofstream ofs ("Mechanical_resonance_shift"+std::to_string(z1)+".dat", std::ofstream::out);
        for(double o = 2*M_PI*f_min; o < 2*M_PI*f_max; o += 2*M_PI*df)
        {
            harvester harv(z1,o);
            cout<<o/(2*M_PI)<<"\t"<<harv.GetAmpl()<<"\t"<<harv.GetPowerExt()<<"\t"<<harv.GetPowerSpring()<<"\t"<<harv.GetPowerDump()<<"\t"<<-harv.GetPower()<<"\t"<<harv.GetPowerDump()+harv.GetPowerSpring()+harv.GetPowerExt()-harv.GetPower()<<"\n";
            ofs<<o/(2*M_PI)<<"\t"<<harv.GetAmpl()<<"\t"<<harv.GetPowerExt()<<"\t"<<harv.GetPowerSpring()<<"\t"<<harv.GetPowerDump()<<"\t"<<-harv.GetPower()<<"\t"<<harv.GetPowerDump()+harv.GetPowerSpring()+harv.GetPowerExt()-harv.GetPower()<<"\n";
        }
        ofs.close();
    }
}

void mechanical_resonance_shift(double freq, double z1_min, double z1_max, double dz1)
{
    std::ofstream ofs ("Mechanical_resonance.dat", std::ofstream::out);
    for(double z1 = z1_min;z1<z1_max+dz1;z1+=dz1)
    {
            harvester harv(z1,freq*2*M_PI);
            cout<<z1<<"\t"<<harv.GetAmpl()<<"\t"<<harv.GetPower()<<"\n";
            ofs<<z1<<"\t"<<harv.GetAmpl()<<"\t"<<harv.GetPower()<<"\n";
    }
    ofs.close();
}
