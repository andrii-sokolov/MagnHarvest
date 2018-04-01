#ifndef OSCILLATOR_H_INCLUDED
#define OSCILLATOR_H_INCLUDED

class harvester
{
    private:
        double dt;

        double Fm;
        double omega;

        double m;
        double beta;
        double k;
        double k1;
        double z1;

        double harvester_ampl;
        double harvester_power;
        double harvester_power_el;
        double harvester_power_spring;
        double harvester_power_dump;
        double harvester_power_ext;

        double Right_part_dimentionless(double khi, double dkhi, double tau);
        double Right_part(double x, double v, double t);
        double Right_part_debug(double x, double v, double t);

    public:
        void RK45();
        void RK45_ndl();
        void RK45_debug();
        double GetAmpl();
        double GetPower();
        double GetPowerEl();
        double GetPowerSpring();
        double GetPowerDump();
        double GetPowerExt();
        void SetDt(double t);
        void SetShaker(double Force_m, double om);
        void SetOscillator(double spring_k, double mass, double shift, double beta);
        double GetMass();
        double GetSpring();
        double GetBeta();
        harvester();
        harvester(double z_shift, double om);

};

void mechanical_resonance_nd(double z1_min, double z1_max, double dz1, double f_min, double f_max, double df);
void mechanical_resonance(double z1_min, double z1_max, double dz1, double f_min, double f_max, double df);
void mechanical_resonance_debug(double z1_min, double z1_max, double dz1, double f_min, double f_max, double df);
void mechanical_resonance_shift(double freq, double z1_min, double z1_max, double dz1);
#endif // OSCILLATOR_H_INCLUDED
