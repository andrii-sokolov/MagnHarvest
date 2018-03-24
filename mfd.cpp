#include "mfd.h"
#include "constants.h"
#include <math.h>
#include <string>
#include <fstream>
#include <iostream>

using namespace std;

void Rectangular_magnet::MFD (double rx, double ry, double rz, double rmx, double rmy, double rmz, double mmx, double mmy, double mmz, double *ans)
{   /**
    Function of elementar magnet flux density
    rx, ry, rz - coordinates of observing point
    rmx, rmy, rmz - coordinates of magnet dipole
    mmx, mmy, mmz - coorginates of magnetization vector
    ans - pointer of answers (1-d 3 elements)
    */
    //coordinates of a magnetic moment
    double x = rx-rmx;
    double y = ry-rmy;
    double z = rz-rmz;
    // Distance between magnetic moment and a point
    double r = sqrt(x*x + y*y +z*z);
    // Scalar multiplication
    double mscalar = mmx*x + mmy*y + mmz*z;
    double Bx = 0;
    double By = 0;
    double Bz = 0;

    if(r==0) //the condition to avoid divizion by zero
    {
        Bx = 0.f;
        By = 0.f;
        Bz = 0.f;
    }
    else
    {
        Bx = mu_0/(4*M_PI)*(3*mscalar*x/(r*r*r*r*r) - mmx/(r*r*r));
        By = mu_0/(4*M_PI)*(3*mscalar*y/(r*r*r*r*r) - mmy/(r*r*r));
        Bz = mu_0/(4*M_PI)*(3*mscalar*z/(r*r*r*r*r) - mmz/(r*r*r));
    }
    ans[0] = Bx;
    ans[1] = By;
    ans[2] = Bz;
}


void Rectangular_magnet::Rectangular_MFD(double rx, double ry, double rz, double* result)
{
    /**
    Function of magnet flux density of parallelepiped
    rx, ry, rz - coordinates of observing point
    a,b,c - length, width and height of parallelepiped
    nmx, nmy, nmz - ort of magnetization vectiors in parallelepiped
    mesh_a, mesh_b, mesh_c - numbers of point through a, b and c dimensions of parallelepiped
    Magn - magnetization density of the material of parallelepiped
    result - the pesults pointer (1-d 3 elements)
    */

    //definition of cylinder's mesh
    double da = magnet_a/mesh_a;
    double db = magnet_b/mesh_b;
    double dc = magnet_c/mesh_c;

    //initialization of B-projections
    double Bx = 0.0;
    double By = 0.0;
    double Bz = 0.0;

    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    double mmz = da*db*dc*Magn;

    double *ans = new double[3];
    ans[0]=0.0;
    ans[1]=0.0;
    ans[2]=0.0;

    for(int i=0; i<=mesh_a; i++)
    {
        for(int j=0; j<=mesh_b; j++)
        {
            for(int k=0; k<=mesh_c; k++)
            {
                x = -magnet_a/2+da*i;
                y = -magnet_b/2+db*j;
                z = -magnet_c/2+dc*k;

                MFD(rx,ry,rz,x,y,z,magnet_n_x*mmz,magnet_n_y*mmz,magnet_n_z*mmz,ans);

                Bx += ans[0];
                By += ans[1];
                Bz += ans[2];
            }
        }
    }
    result[0]=Bx;
    result[1]=By;
    result[2]=Bz;
    delete ans;
}

double Rectangular_magnet::Rectangular_Flux_flat(double z, double surf_a, double surf_b)
{
    /**
    Function defines the Flux of magnet flux density through the rectangular surface perpendicular to the z-axis
    z - z coordinate of the surface
    surf_a, surf_b - dimentions of the surface
    mesh_surf_a, mesh_surf_b - number of points through the a-, and b- dimentions of the surface
    a,b,c - length, width and height of parallelepiped
    nmx, nmy, nmz - ort of magnetization vectiors in parallelepiped
    mesh_a, mesh_b, mesh_c - numbers of point through a, b and c dimensions of parallelepiped
    magnez - magnetization density of the material of parallelepiped
    */

    double dx = 2*surf_a/mesh_surf_a;
    double dy = 2*surf_b/mesh_surf_b;
    double Flux = 0.0;
    double dS = 0.0;
    double *ans=new double[3];

    double Bz = 0.0;
    for(double x = 0.0; x<= surf_a/2; x+=dx)
    {
        for(double y = 0.0; y<= surf_b/2; y+=dy)
        {
            dS = dx*dy;
            Rectangular_magnet::Rectangular_MFD(x+dx/2, y+dy/2, z, ans);
            Bz = ans[2];
            Flux += Bz*dS;
        }
    }

    Flux = Flux*4;
    return(Flux);
}

double Rectangular_magnet::Rectangular_Flux(double x, double y, double z, double a1, double a2, double a3, double surf_a, double surf_b)
{
    double da = surf_a/mesh_surf_a;
    double db = surf_b/mesh_surf_b;
    double dS = da*db;
    int i = 0;

    /**generate plane*/
    for(double sx = -surf_a/2.0; sx < (surf_a/2.0); sx+=da)
    {
        for(double sy = -surf_b/2.0; sy < (surf_b/2.0); sy+=db)
        {
            x_surf[i] = da/2 + sx;
            y_surf[i] = db/2 + sy;
            z_surf[i] = 0.0;
            i++;
        }
    }
    double ex,vi,ze;
    /** rotate plane _ x*/
    for(i=0; i<mesh_surf_a*mesh_surf_b; i++)
    {
        ex  = x_surf[i];
        vi  = y_surf[i];
        ze  = z_surf[i];
        x_surf[i] = ex;
        y_surf[i] = vi*cos(a1) + ze*sin(a1);
        z_surf[i] =-vi*sin(a1) + ze*cos(a1);
    }
    /** rotate plane _ y*/
    for(i=0; i<mesh_surf_a*mesh_surf_b; i++)
    {
        ex  = x_surf[i];
        vi  = y_surf[i];
        ze  = z_surf[i];
        x_surf[i] += ex*cos(a2) - ze*sin(a2);
        y_surf[i] += vi;
        z_surf[i] += ex*sin(a2) + ze*cos(a2);
    }
    /** rotate plane _ z*/
    for(i=0; i<mesh_surf_a*mesh_surf_b; i++)
    {
        ex  = x_surf[i];
        vi  = y_surf[i];
        ze  = z_surf[i];
        x_surf[i] += ex*cos(a3) + vi*sin(a3);
        y_surf[i] +=-ex*sin(a3) + vi*cos(a3);
        z_surf[i] += ze;
    }

    /**shift plane*/
    for(i=0; i<mesh_surf_a*mesh_surf_b; i++)
    {
        x_surf[i] += x;
        y_surf[i] += y;
        z_surf[i] += z;
    }

    /**calculation of the normal to the plane*/
    double x1 = x_surf[0];
    double x2 = x_surf[(int)(mesh_surf_a*mesh_surf_b/2)];
    double x3 = x_surf[mesh_surf_a*mesh_surf_b - 1];

    double y1 = y_surf[0];
    double y2 = y_surf[(int)(mesh_surf_a*mesh_surf_b/2)];
    double y3 = y_surf[mesh_surf_a*mesh_surf_b - 1];

    double z1 = z_surf[0];
    double z2 = z_surf[(int)(mesh_surf_a*mesh_surf_b/2)];
    double z3 = z_surf[mesh_surf_a*mesh_surf_b - 1];

    double x_norm = -y2*z1 + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3;
    double y_norm =  x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3;
    double z_norm = -x2*y1 + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3;

    double Flux =0.0;

    double* res = new double(3);

    for(i = 0; i<mesh_surf_a*mesh_surf_b; i++)
    {
        Rectangular_magnet::Rectangular_MFD(x_surf[i],y_surf[i],z_surf[i],res);
        Flux+= (res[0]*x_norm + res[1]*y_norm + res[2]*z_norm)*dS/sqrt(x_norm*x_norm + y_norm*y_norm + z_norm*z_norm);
    }

    return(Flux);
/*    for(i=0; i<mesh_surf_a*mesh_surf_b; i++)
    {
        cout<<x_surf[i]<<"\t"<<y_surf[i]<<"\t"<<z_surf[i]<<endl;
    }
*/
}

void Rectangular_magnet::Generate_Flux_flat_file(double z_min, double z_max, double dz, double a_max, double a_min,double da,   string filename)
{
    /**
    This method generates the file of a magnet flux density
    it calculates the flux through different square surfaces
    from the point z_min to the point z_max with the step dz
    */

    double Flux=0.0;
    std::ofstream ofs (filename, std::ofstream::out);

    ofs<<"##\t";
    cout<<"##\t";
    for(double a = a_min; a<=(a_max+da);a+=da)
        {
            ofs<<a<<"\t";
            cout<<a<<"\t";
        }
    ofs<<endl;
    cout<<endl;
    for(double z = z_min; z<z_max; z+=dz)
    {
        ofs<<z<<"\t";
        cout<<z<<"\t";
        for(double a = coil_a_min*2.0;a<=coil_a_max*2.01;a+=da)
        {
            Flux = Rectangular_magnet::Rectangular_Flux_flat(z, a, a);
            ofs<<Flux<<"\t";
            cout<<Flux<<"\t";
        }
        ofs<<endl;
        ofs<<endl;

    }
    ofs.close();
}

void Rectangular_magnet::Set_Geometry(double ma, double mb, double mc)
{
    magnet_a = ma;
    magnet_b = mb;
    magnet_c = mc;
}

void Rectangular_magnet::Set_Magnetization(double nx, double ny, double nz, double mag)
{
    magnet_n_x = nx;
    magnet_n_y = ny;
    magnet_n_z = nz;
    Magn = mag;
}

void Rectangular_magnet::Set_Mesh(int mesa, int mesb, int mesc)
{
    mesh_a = mesa;
    mesh_b = mesb;
    mesh_c = mesc;
}

void Rectangular_magnet::Set_Mesh_surf(int mesa, int mesb)
{
    mesh_surf_a = mesa;
    mesh_surf_b = mesb;
}

Rectangular_magnet::Rectangular_magnet()
{
    magnet_a = 0.0025;
    magnet_b = 0.0025;
    magnet_c = 0.0020;
    magnet_n_x = 0.0;
    magnet_n_y = 0.0;
    magnet_n_z = 1.0;
    mesh_a = 200;
    mesh_b = 200;
    mesh_c = 200;
    mesh_surf_a = 10;
    mesh_surf_b = 10;
    Magn = 796.0e3;

    x_surf = new double[mesh_surf_a*mesh_surf_b];
    y_surf = new double[mesh_surf_a*mesh_surf_b];
    z_surf = new double[mesh_surf_a*mesh_surf_b];
}

Rectangular_magnet::~Rectangular_magnet()
{
    delete x_surf;
    delete y_surf;
    delete z_surf;
}
