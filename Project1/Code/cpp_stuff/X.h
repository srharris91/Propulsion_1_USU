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
        }
        // variables
        double Vr, Vv,r,v,m;
        double Vr_dot, Vv_dot,r_dot,v_dot,m_dot;
        double beta,gamma,C_D,A_ref,Isp,mu,rho,V_inf;
        double x,y,t;

        // functions
        void x_dot(double& F_thrust);
        X Trapz(double& F_thrust,double& dt);
        void x_set(X& a);
        void xy();
        void print();
        void output_line(stringstream &a);
};
#endif
