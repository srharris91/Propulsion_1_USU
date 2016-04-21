#ifndef RP_H
#define RP_H
#include <cmath>
#include "CylindricalPort.h"
#include "RP.h"

RP::RP(){       // default constructor
    r=0.015;     // m
    r_dot=0.0000001   ;// dm/ds
    P0=0.0001;// Pa
    P0_dot=0.0000001;  // dPa/ds
}
//double RP::A_burn   (CylindricalPort a) { return 2.*M_PI*r*a.L0;    }
double RP::A_burn   (CylindricalPort a,int b) { 
    if (b == 0) return 2.*M_PI*r*a.L0; 
    else {
//        double N=double(b);
        double s=r-a.d0/2.;
        return a.N*M_PI*
            (
             (a.D0*a.D0-(a.d0+2.*s)*(a.d0+2.*s))/2.
             +
             ((a.L0-2.*s)*(a.d0+2.*s))
            );
    }
}
//double RP::V_c      (CylindricalPort a) { return M_PI*r*r*a.L0;     }
double RP::V_c      (CylindricalPort a,int b) { 
    if (b == 0) return M_PI*r*r*a.L0; 
    else{
//        double N=double(b);
        double s=r-a.d0/2.;
        return a.N*M_PI/4.
            * ((a.d0+2.*s)*(a.d0+2.*s) * (a.L0-2.*s)
                    + a.D0*a.D0*2.*s
              );
    }
}


#endif
