#ifndef PRANDTL
#define PRANDTL
#include <iostream>
#include <cmath>
#include "Prandtl_Meyer.h"
using namespace std;


void    Prandtl_Meyer::set_isentropic(int  a)  { isen=a;   }
// set isen = 0 for isentropic relationships
// set isen = 1 for non-isentropic relationships
// set isen = 2 for Rayleigh Pitot Equation
// set isen = 3 for Oblique Shock Waves
// set isen = 4 for Prandtl-Meyer  Equations
void    Prandtl_Meyer::set_Rg      (double a)  { Rg=a;     }
void    Prandtl_Meyer::set_F(){

    if (isen==0){//isentropic nozzle
        F   =   (1./M)* 
            pow(((2./(gamma+1))*(1.+((gamma-1.)/2.)*M*M)), 
                    ((gamma+1.)/(2.*(gamma-1.))))- 
            (A/A_star); 
    }
    else if (isen==1){//non-isentropic nozzle, ahead of shock wave
        F   =   2.
            /((gamma+1.)*
                    pow(gamma*M*M-((gamma-1.)/2.),
                        1./(gamma-1.)))
            *pow(pow((((gamma+1.)/2.)*M),2)/
                    (1.+((gamma-1.)/2.)*M*M),
                    gamma/(gamma-1.))
            -(P02/P0);
    }
    else if (isen==2){// rayleigh pitot equations
        if (M>1){
            F   =   (
                    pow(((gamma+1.)/2.)*M*M,gamma/(gamma-1.))
                    )
                /(
                        pow(((2.*gamma*M*M)/(gamma + 1.))
                            -((gamma - 1.)/(gamma+1.)),1./(gamma-1.))
                 )
                - (P02/P)
                ;
        }
        else if (M<=1){
            F   =   pow(1. + (M*M*(gamma-1.))/(2.),gamma/(gamma-1.))
                - (P02/P)
                ;
        }
    }
    else if (isen==3){// oblique shock waves
        a = (1.+((gamma-1.)/2.)*M*M)*tan(theta);
        b = (M*M-1.);
        c = (1.+((gamma+1.)/2.)*M*M)*tan(theta);
        xj = tan(beta);
        F   =   2.*a*xj*xj*xj
            - b*xj*xj - 1.;

    }
    else if (isen==4){// prandtl-meyer equations
        F   =   (theta + V(M) - V(M2));
    }
}
void    Prandtl_Meyer::set_dFdM(){
    if (isen==0){
        //if isen=0 then assume isentropic nozzle
        dFdM=   (pow(2., 
                    (1.-3.*gamma)/(2.-2*gamma)))* 
            (M*M-1.)/ 
            (M*M*(2.+M*M*(gamma-1.)))* 
            pow(((1.+((gamma-1.)/2.)*M*M)/(gamma+1.)),
                    ((gamma+1.)/(2*(gamma-1.))));
    }
    //if a=1 then non-isentropic nozzle, ahead of shock wave
    else if (isen==1){
        dFdM=   -(pow(2.,
                    3.-((2.*gamma)/(gamma-1.)))*
                gamma*
                pow(M*M-1.,2)*
                pow(pow((((gamma+1.)/2.)*M),2)/
                    (1.+((gamma-1.)/2.)*M*M),
                    gamma/(gamma-1.))*
                pow(0.5 + gamma*(M*M-0.5),-1./(gamma-1.)))
            /
            ((gamma+1.)*M*(2.+M*M*(gamma-1.))*(1.+gamma*(2*M*M-1.)));
    }
    else if (isen==2){// rayleigh pitot equations
        if (M>1){
            dFdM    =   (
                    gamma * M * (2.*M*M-1.) * pow((M*M*(gamma+1.))/2.,1./(gamma-1.))
                    )
                /(
                        pow(((2.*gamma)/(gamma+1.))*M*M - (gamma-1.)/(gamma+1.),(gamma/(gamma-1.)))
                 );
            //                std::cout<<"M>1"<<std::endl;
        }
        else if (M<=1){
            dFdM    =   (M*(gamma))
                * pow(1.+((gamma-1.)/2.)*M*M,(1./(gamma-1.)));
            //                std::cout<<"M<=1"<<std::endl;
        }
    }
    else if (isen==3){
        dFdM    =   3.*a*xj*xj
            - 2.*b*xj + c;
    }
    else if (isen==4){
        dFdM    =   (1./M2)
            *((sqrt(M2*M2 - 1.))
                    /(1.+((gamma-1.)/2.) *M2*M2))
            ;
    }
}


//calculate the next iteration Mach number  in the newton solver iteration
void Prandtl_Meyer::Newton_iter(){ M-=F/dFdM;    } 
// Newton solver for oblique shock waves
void Prandtl_Meyer::Newton_iter_oblique(){ 
    xnew = F/dFdM;
    beta = atan(xnew);
    //        cout<<iter<<" iteration and beta = "<<beta*180./M_PI<<endl;
}
void Prandtl_Meyer::Newton_iter_PM(){ M2+=F/dFdM;    } 

//calculate the error (and set the new F and dFdM values) on each iteration of Mach number in the newton solver
double Prandtl_Meyer::error(){ 
    set_F();    // set the F = 0 value for the current iteration
    set_dFdM(); // also calculate the dF/dM value for the current iteration
    iter++;     // each time this is run count up the iterations by one increment
    //        std::cout<<iter<<" iteration Ma = "<<M<<std::endl;
    //        return (std::fabs(F) /(A/A_star));    //alternative iterative solver convergence criteria
    //        cout<<"error = "<<F<<" / "<<dFdM<<endl;
    if (isen!=3){
        return (std::abs(F/dFdM));
    }
    else if (isen==3){
        return (std::abs((F/dFdM-xnew)/(F/dFdM)));
    }
    else
        return 999999999.;
}


//calculate the converged value of Mach number using the Newton iteration solver
void Prandtl_Meyer::Newton(){
    if (isen!=3 && isen!=4){
        // for loop to converge on the correct M value
        for (;error()>tol;) Newton_iter();
    }
    else if (isen==3){
        // for loop to converge on the correct beta value (for oblique shock waves)
        for (;error()>tol;) Newton_iter_oblique();
    }
    else if (isen==4){
        // for loop to converge on the correct M2 value for Prandtl_Meyer equations
        for (;error()>tol;) Newton_iter_PM();

    }
    else {
        cout<<"you are in big trouble. Fix your isen value to Newton solve correctly."<<endl;
    }
}

//calculate the temperature for an isentropic nozzle
void Prandtl_Meyer::Temperature(){
    T = T0/
        (1.+(gamma-1.)/2.*M*M);
}
void Prandtl_Meyer::Pressure(){
    P = P0/
        pow((1.+(gamma-1.)/2.*M*M),
                (gamma/(gamma-1.)));
}
// calculate the choked mass flow with the given Gamma, Rg, P0 and T0
void Prandtl_Meyer::calc_m_dot(){
    m_dot = (P0/sqrt(T0))
        * A_star
        * sqrt((gamma/Rg)
                * pow(2./(gamma+1.),
                    (gamma+1.)/(gamma-1.)));
}
// oblique shock wave M to M2
void Prandtl_Meyer::calc_M2_oblique(){
    double Mn2;
    Mn2 = sqrt((1.+((gamma-1.)/2.)*pow(M*sin(beta),2))
            /(gamma*pow(M*sin(beta),2) - (gamma-1.)/2.));
    M2  =   Mn2/(sin(beta-theta));
}
// Prandtl-Meyer equation V(M)
double Prandtl_Meyer::V(double Mach){
    return sqrt((gamma+1.)/(gamma-1.))
        *atan(sqrt((gamma-1.)/(gamma+1.)*(Mach*Mach-1.)))
        - atan(sqrt(Mach*Mach-1.));
}
void Prandtl_Meyer::Stagnation_Pressure(){//isentropic relation to get P0 from P
    P0 = P*
        pow((1.+(gamma-1.)/2.*M*M),
                (gamma/(gamma-1.)));
}
void Prandtl_Meyer::Pressure_Oblique(){
    P2  =   P
        *(1.+(2.*gamma)/(gamma+1.)
                *(pow(M*sin(beta),2)-1.));
}
void Prandtl_Meyer::Temperature_Oblique(){
    T2  =   T*(1.+(2.*gamma)/(gamma+1.)*(pow(M*sin(beta),2)-1.))
        *((2.+(gamma-1.)*(pow(M*sin(beta),2)))
                /((gamma+1.)*pow(M*sin(beta),2)));
}
void Prandtl_Meyer::P0_Oblique(){
    P02   =   P0*2.
        /((gamma+1.)*
                pow(gamma*M*sin(beta)*M*sin(beta)-((gamma-1.)/2.),
                    1./(gamma-1.)))
        *pow(pow((((gamma+1.)/2.)*M*sin(beta)),2)/
                (1.+((gamma-1.)/2.)*M*sin(beta)*M*sin(beta)),
                gamma/(gamma-1.));
}
void Prandtl_Meyer::calc_Oblique(){
    Newton();
    calc_M2_oblique();
    Pressure_Oblique();
    Temperature_Oblique();
    Stagnation_Pressure();
    P0_Oblique();
}
void Prandtl_Meyer::P2_Prandtl_Meyer(){
    P2  =   P
        *pow((1.+((gamma-1.)/2.)*M*M)
                /(1.+((gamma-1.)/2.)*M2*M2)
                ,(gamma)/(gamma-1.));
}
void Prandtl_Meyer::T2_Prandtl_Meyer(){
    T2  =   T
        *(1.+((gamma-1.)/2.)*M*M)
        /(1.+((gamma-1.)/2.)*M2*M2)
        ;
}
void Prandtl_Meyer::mu_Prandtl_Meyer(){
    mu1 =   asin(1./M);
    mu2 =   asin(1./M2)-theta;
}
void Prandtl_Meyer::calc_Prandtl_Meyer(){
    Newton();           // calculate M2 using Prandtl-Meyer equations and a Newton solver iterative method
    P2_Prandtl_Meyer(); // Calculate P2 using Prandtl-Meyer equations
    T2_Prandtl_Meyer(); // Calculate T2 using Prandtl-Meyer equations
    mu_Prandtl_Meyer(); // Calculate mu values before and after using Prandtl-Meyer equations
}
double Prandtl_Meyer::Thrust(double conical_nozzle_theta,double Pamb){ // Thrust as a function of conical nozzle angle (degrees) and pressure at altitude
    
//    std::cout<<conical_nozzle_theta<<" "<<Pamb<<" "<<P0<<" "<<P<<" "<<gamma<<" "<<A_star<<" "<<std::endl;
    return ((1.+cos(conical_nozzle_theta*M_PI/180.))/2.)
        *
        gamma*P0*A_star
        *
        sqrt( 
                (
                 2.
                 /
                 (gamma-1.)
                )
                *
                (
                 pow(2./(gamma+1.),
                     (gamma+1.)/(gamma-1.))
                )
            )
        *
        sqrt(
                (
                 1.
                 -
                 pow(P/P0,
                     (gamma-1.)/gamma)
                )
            )
        +
        A*(P-Pamb)
        ;
}



#endif
