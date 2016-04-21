#include "DEFS.h"
#include "Prandtl_Meyer.h"
#include "CylindricalPort.h"
#include "RP.h"
#include "Solver.h"
int main(){

/*
    RP a;
    CylindricalPort b;
    Prandtl_Meyer PM,PM2;
    b.Mw=197.231;// g/mol
    b.Rg = 8.314/(b.Mw/1000.);
    b.A_exit=3.8*3.8*M_PI/4.;
    b.A_star=b.A_exit/7.78;
    b.gamma=1.02;
    b.L0=34.57;
    b.D0=3.66;
    b.d0=1.7;
    b.n=0.172;
    b.a=0.00192/pow(1000.,b.n); // m/(s*Pa^n)
    b.rho_propellant=1760.;
    a.r=b.d0/2.;
    b.T0=24000.;
    b.k=0;
//    cout<<"A_star = "<<b.A_star<<endl;
//    cout<<"A_exit = "<<b.A_exit<<endl;
//    cout<<"launch mass = "<<(b.D0*b.D0*M_PI/4.*b.L0-a.V_c(b))*b.rho_propellant<<endl;
    PM.A = a.V_c(b,0)/b.L0;
    PM.A_star=b.A_star;
    PM.gamma=b.gamma;
    PM.P0 = a.P0;
    PM.T0 = b.T0;
    PM.Rg = b.Rg;
    PM.M = 0.15; // guess subsonic region
    PM.set_isentropic(0);
    PM2=PM;
    PM2.A= b.A_exit;
    PM2.M= 2.;
    Solver S(b,a,PM,PM2);
//    cout<<"launch mass = "<<(b.D0*b.D0*M_PI/4.*b.L0-a.V_c(b,S.bates))*b.rho_propellant<<endl;
//    cout<<"Burn area = "<<(S.RP0.A_burn(S.CP,S.bates))<<endl;



    
    cout.precision(16);
    S.writeout();
//    while (S.t1<125.) {
    while (1){
        S.t1=S.t0+.01;
        S.next_step();
        if (S.RP0.r>=S.CP.D0/2.){ break; }
    }

*/

/*
    // non-erosive vs erosive burning
    RP a;
    CylindricalPort b;
    b.L0=0.35;
    b.D0=0.076;
    b.d0=0.03;
    b.rho_propellant=1314.;
    b.A_star=0.0001887;
    b.A_exit=4.*b.A_star;
    b.gamma=1.2;
    b.Mw=24.26;
    b.Rg = 8.314/(b.Mw/1000.);
    b.T0=2000.;
    b.M_crit=0.11;
    b.k=2.25;
    b.n=0.188;
    b.a=0.00178/pow(1000.,b.n);// m/(s*Pa^n)
    Prandtl_Meyer Port,Nozzle;
    Port.A = a.V_c(b,0)/b.L0;
    Port.A_star=b.A_star;
    Port.gamma=b.gamma;
    Port.P0 = a.P0;
    Port.T0 = b.T0;
    Port.Rg = b.Rg;
    Port.M = 0.15; // guess subsonic region
    Port.set_isentropic(0);
    Nozzle=Port;
    Nozzle.A= b.A_exit;
    Nozzle.M= 2.6;
    Solver S(b,a,Port,Nozzle);

    cout.precision(16);
    S.writeout();
//    S.CP.k=0.;

    while (1){
        S.t1=S.t0+.001;
        S.next_step();
        if (S.RP0.r>=S.CP.D0/2.){ break; }
    }
*/

///*
    // Part 1 cylindrical port
    // initialize parameters
    RP a;
    CylindricalPort b;
    Prandtl_Meyer Port,Nozzle;
    Port.A = a.V_c(b,0)/b.L0;
    Port.A_star=b.A_star;
    Port.gamma=b.gamma;
    Port.P0 = a.P0;
    Port.T0 = b.T0;
    Port.Rg = b.Rg;
    Port.M = 0.16; // guess subsonic region
    Port.set_isentropic(0);
    Nozzle=Port;
    Nozzle.A= b.A_exit;
    Nozzle.M= 2.6;
    Solver S(b,a,Port,Nozzle);
    S.bates=0;  // cylindrical port
    S.Port.Newton();
    // solve and output
    cout.precision(16);
//    S.CP.k=0.; // non-erosive burn
    S.writeout();
    while (1){
        S.t1=S.t0+.0001;
        S.next_step();
        if (S.RP0.r>=S.CP.D0/2.){ break; }
    }
//*/
/*
    // Part 2 bates grain
    // initialize parameters
    RP a;
    CylindricalPort b;
    Prandtl_Meyer Port,Nozzle;
    Port.A = a.V_c(b,0)/b.L0;
    Port.A_star=b.A_star;
    Port.gamma=b.gamma;
    Port.P0 = a.P0;
    Port.T0 = b.T0;
    Port.Rg = b.Rg;
    Port.M = 0.16; // guess subsonic region
    Port.set_isentropic(0);
    Nozzle=Port;
    Nozzle.A= b.A_exit;
    Nozzle.M= 2.6;
    Solver S(b,a,Port,Nozzle);
    S.bates=1;      // solve using bates grain
    S.Port.Newton();
    // solve and output
    cout.precision(16);
    S.CP.k=0.; // non-erosive burn
    S.writeout();
    while (1){
        S.t1=S.t0+.0001;
        S.next_step();
        if (S.RP0.r>=S.CP.D0/2.){ break; }
    }
*/

/*
    // Part 2 bates
    // initialize parameters
    RP a;
    CylindricalPort b;
    b.T0 = 3300.; // higher values of flame temp
    Prandtl_Meyer Port,Nozzle;
    Port.A = a.V_c(b,0)/b.L0;
    Port.A_star=b.A_star;
    Port.gamma=b.gamma;
    Port.P0 = a.P0;
    Port.T0 = b.T0;
    Port.Rg = b.Rg;
    Port.M = 0.16; // guess subsonic region
    Port.set_isentropic(0);
    Nozzle=Port;
    Nozzle.A= b.A_exit;
    Nozzle.M= 2.6;
    Solver S(b,a,Port,Nozzle);
    S.bates=1;
    S.Port.Newton();
    // solve and output
    cout.precision(16);
    S.CP.k=0.; // non-erosive burn
    S.writeout();
    while (1){
        S.t1=S.t0+.0001;
        S.next_step();
        if (S.RP0.r>=S.CP.D0/2.){ break; }
    }
*/




    return 0;
}
