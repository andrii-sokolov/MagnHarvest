#ifndef CONSTANTS_H_INCLUDED
#define CONSTANTS_H_INCLUDED

/**magnet constant*/
const double mu_0 = 1.25663706e-6;

/** Spring constant*/
const double k = 2665.7;
/** Air dumphing constant*/
const double beta = 8.8e-3;

/** The temporary block
*   time constants for the numerical method*/
//double dt = 0.000001;
//double t_max = 1.5;
//int i_max = 100;

/** Initial shift of a coil corresponding to the center of a magnet*/
const double z0 = 0.0022;
/** Parametric shift*/
//double z1 = 0.000;

/** Geometrical parameters of a magnet*/
//const double magnet_a = 0.0025;
//const double magnet_b = 0.0025;
//const double magnet_c = 0.002;

/** The density of given magnet*/
//const double dens = 7400.0;
/** The volume of magnet*/
//const double Volume = magnet_a*magnet_b*magnet_c;
/** The mass of magnet*/
//const double m = 9.83e-5;

/** Amplitude of the external force*/
//const double Fm = m*1.0*9.8;

/** Meshing of a surface in Flux calcucations*/
const int mesh_surf_x = 10;
const int mesh_surf_y = 10;

/** Meshing of a force in Force calculations*/
const int mesh_force = 5;

/** Maximal size of a FLAT coil*/
const double coil_a_max = 0.0028;
/** Minimak size of a FLAT coil*/
const double coil_a_min = 0.00009;
/** The number of loops in a coil*/
const int coil_n_of_loops = 144;
/** Resistance of a coil*/
const double coil_resistance = 192.0;

/** Interpolation constants of the flux*/
const double a_constant[5] = {1.77062722e-08, -1.92740770e-05, 2.62380282e-03, -1.03433942e-01, 1.40998442e+00};
const double b_constant[5] = {2.16792821e-03, -2.12101602e-01, 2.32880357e+01, -1.09021961e+03, 1.77267708e+04};

double a_parameter(double r);
double b_parameter(double r);
double Flux_square_interp(double r, double z);
double Force_EM_interp(double z, double v);
#endif // CONSTANTS_H_INCLUDED
