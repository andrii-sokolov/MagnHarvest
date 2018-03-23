#ifndef MFD_H_INCLUDED
#define MFD_H_INCLUDED

#include <string>
using namespace std;
//double MFD (double rx, double ry, double rz, double rmx, double rmy, double rmz, double mmx, double mmy, double mmz, double *ans);
//double Rectangular_MFD(double rx, double ry, double rz, double a, double b, double c, double nmx, double nmy, double nmz,int mesh_a, int mesh_b, int mesh_c, double Magn, double* result);

class Rectangular_magnet {
    private:
        double magnet_a;
        double magnet_b;
        double magnet_c;
        double magnet_n_x;
        double magnet_n_y;
        double magnet_n_z;
        int mesh_a;
        int mesh_b;
        int mesh_c;
        double Magn;
        int mesh_surf_a;
        int mesh_surf_b;
        void MFD (double rx, double ry, double rz, double rmx, double rmy, double rmz, double mmx, double mmy, double mmz, double *ans);
    public:
        void Rectangular_MFD(double rx, double ry, double rz, double* result);

        double Rectangular_Flux_flat(double z, double surf_a, double surf_b);
        double Rectangular_Flux(double x, double y, double z, double nx, double ny, double nz, double surf_a, double surf_b); //not implemented yet

        void Set_Geometry(double ma, double mb, double mc);
        void Set_Magnetization(double nx, double ny, double nz, double magn);
        void Set_Mesh(int mesa, int mesb, int mesc);
        void Set_Mesh_surf(int mesa, int mesb);

        void Generate_Flux_flat_file(double z_min, double z_max, double dz, double a_max, double a_min,double da, string filename);

        Rectangular_magnet();
};

#endif // MFD_H_INCLUDED
