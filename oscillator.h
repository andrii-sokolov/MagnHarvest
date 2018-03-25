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
        double z1;

        double harvester_ampl;
        double harvester_power;

        double Right_part_dimentionless(double khi, double dkhi, double tau);


    public:
        void RK45();
        double GetAmpl();
        double GetPower();
        void SetDt(double t);
        void SetShaker(double Force_m, double om);
        void SetOscillator(double spring_k, double mass, double shift, double beta);
        double GetMass();
        double GetSpring();
        double GetBeta();
        harvester();

};

void mechanical_resonance(double z1_min, double z1_max, double dz1, double f_min, double f_max, double df);
#endif // OSCILLATOR_H_INCLUDED
