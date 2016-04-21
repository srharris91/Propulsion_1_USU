#ifndef CYLINDRICALPORT_H
#define CYLINDRICALPORT_H
#include <cmath>
#include "CylindricalPort.h"

CylindricalPort::CylindricalPort()
{
            L0              = .035;      // m
            D0              = .066;      // m
            d0              = .03;       // m
            rho_propellant  = 1260.;     // kg/m^3
            A_star          = .0001887;  // m^2   
            A_exit          = 4.*A_star; // m^2
            theta_exit      = 20.*M_PI/180.;   // radians
            gamma           = 1.18;      //
            Mw              = 23.;       // kg/kg-mol
            Rg              = 8.314/(Mw/1000.);
            T0              = 2900.;     // K
            n               = 0.16;
            a               = 0.00132/pow(1000.,n);     // m/(s*Pa^n)
            M_crit          = 0.3;
            k               = 0.2;
            N               = 3;
}





#endif
