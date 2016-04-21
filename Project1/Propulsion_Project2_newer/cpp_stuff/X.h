#ifndef X_H
#define X_H
#include "DEFS.h"
class X{
    public: 
        X(){
            Vr      =   0.;
            Vv      =   0.;
            r       =   0.;
            v       =   0.;
            m       =   0.;
            Vr_dot  =   0.;
            Vv_dot  =   0.;
            r_dot   =   0.;
            v_dot   =   0.;
            m_dot   =   0.;
            t       =   0.;
//            mu      =   3.9860044E14;      // mu for earth
            mu      =   3.9860000E14;      // mu for earth
        }
        // variables
        double Vr, Vv,r,v,m;
        double Vr_dot, Vv_dot,r_dot,v_dot,m_dot;
        double gamma,Isp,mu;
        double Thrust,dt;
        double x,y,t;
        double a,e; // elliptical orbit
        double DeltaV;  // accumulated DeltaV value for the whole system

        // functions
        void x_dot(double& F_thrust);
        X Trapz (double& F_thrust,double& dt);
        X RK4   (double& F_thrust,double& dt);
        X function(int a,X& Old_f, double& Thrust_f, double& dt_f);
        double  HomannTransfer(int a,int b,X start, X coast, X final, X Old, X New);
        double  ShootingMethod(int a,int b,X start, X coast, X final, X Old, X New);
        void    HomannTransfer_output(int a,int b,X start, X coast, X final, X Old, X New, ofstream& myfile, ofstream& myfile2);
        void xy();
        void print();
        void output_line(stringstream &a);
};
#endif
