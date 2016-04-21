#ifndef SOLVER_H
#define SOLVER_H
#include "CylindricalPort.h"
#include "Prandtl_Meyer.h"
#include "RP.h"
#include "Solver.h"
#include <cmath>
#include <iostream>

Solver::Solver(){
    t0=0.;
    t1=0.;
    bates = 0;
}
Solver::Solver(CylindricalPort a, RP b,Prandtl_Meyer c,Prandtl_Meyer d){
    CP = a;
    RP0= b;
    RP1= b;
    Port = c;
    Nozzle = d;

    t0=0.;
    t1=0.;

    bates = 0;
}

void Solver::calc_RP1_dot(){
//    RP1.r_dot = CP.a * pow(RP0.P0,CP.n)/pow(1000.,CP.n);
    Port.A = RP0.V_c(CP,bates)/CP.L0;
    RP1.r_dot = CP.a * pow(RP0.P0,CP.n)
        * ((
            1.+CP.k*(Port.M/CP.M_crit))
        / (1.+CP.k));
    RP1.P0_dot = (RP0.A_burn(CP,bates) * RP1.r_dot / RP0.V_c(CP,bates) )
//    RP1.P0_dot = (2. * RP1.r_dot / RP0.r )
        * (CP.rho_propellant*CP.Rg*CP.T0-RP0.P0)
        - ((CP.A_star/RP0.V_c(CP,bates))
        * (RP0.P0)
        * sqrt(CP.gamma*CP.Rg*CP.T0
                * pow(2./(CP.gamma+1.),(CP.gamma+1.)/(CP.gamma-1.))));
}
void Solver::calc_RP1(){
    RP1.P0 = RP0.P0 + (t1-t0) * RP1.P0_dot;
    RP1.r  = RP0.r  + (t1-t0) * RP1.r_dot;
}
void Solver::writeout(){
//    std::cout<<std::endl;
    if (isnan(RP1.P0)) {
        std::cerr<<"nan on RP1.P0"<<std::endl;
        return;
    }
    std::cout<<t1<<" "<<RP1.P0_dot  <<" "<<RP1.r_dot<<" "<<RP1.P0<<" "<<RP1.r<<" "<<RP1.A_burn(CP,bates)<<" "<< ((CP.D0*CP.D0*M_PI/4.*CP.L0)-RP1.V_c(CP,bates))*CP.rho_propellant<<" "<<((RP1.V_c(CP,bates)-RP0.V_c(CP,bates))*CP.rho_propellant)/(t1-t0)<<" "<<Port.m_dot<<" "<<Nozzle.Thrust(20.,0.)<<" "<<Nozzle.Thrust(20.,0.)/(Port.m_dot*9.81)<<" "<<Port.M<<" "<<" "<<Nozzle.M<<std::endl;
}
void Solver::reset(){
    RP0 = RP1;
    t0  = t1;
}
void Solver::next_step(){
    calc_RP1_dot();
    calc_RP1();
    Port.P0=RP1.P0;
    Nozzle.P0=RP1.P0;
    Port.Newton();
    Port.Pressure();
    Port.calc_m_dot();
    Nozzle.Newton();
    Nozzle.Pressure();
    Nozzle.calc_m_dot();
    writeout();
    reset();
}
    

#endif
