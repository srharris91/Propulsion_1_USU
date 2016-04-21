class CylindricalPort{
    public:
        //Fuel grain geometry
        double L0,D0,d0,rho_propellant;
        // Nozzle Geometry
        double A_star,A_exit,theta_exit;
        // Combustion gas properties
        double gamma,Mw,T0,Rg;
        // Burn Parameters
        double a,n,M_crit,k;
        // how many bates grains
        int N;
       // constructor
       CylindricalPort();
};

