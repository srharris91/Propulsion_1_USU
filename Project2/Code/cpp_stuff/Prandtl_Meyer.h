// this is a newton solver for computing the area-Mach number relation
class Prandtl_Meyer{
public:
    // public values of gamma, A as a function of x, A at the throat, Mach number, the Newton iteration tolerance
    double gamma, A, A_star, M, tol; 
    //public values for F = 0, dF/dM, and the number of iterations for the Newton solver
    double F, dFdM;
    int iter;
    //public values for P0 and T0 (stagnation pressure and temperature for an isentropic nozzle)
    //and local P and T (if an isentropic nozzle)
    double P0,T0,P,T;
    // non-isentropic values
    double P02;
    // x,yp,ym (positive and negative) values within the nozzle
    double x,yp,ym;
    // set a=1 for non-isentropic and a=0 for isentropic nozzle
    double isen;
    // ideal gas Rg property
    double Rg;
    // mass flow m_dot
    double m_dot;
    // oblique shock wave values of x,a,b,c and xnew and beta
    double xj,a,b,c,xnew,beta;
    // oblique shock wave density, pressure and temperature before and after
    double rho1,rho2,P2,T2;
    // M2 and theta and mu1 and mu2 of fan for Prandtl-Meyer equations (theta is reused for oblique shock waves as well)
    double M2,theta,mu1,mu2;
    //constructors
    Prandtl_Meyer(){
        gamma   =   1.18;
        A_star  =   0.0001887;
        A       =   4.*A_star;
        M       =   4.0; //initial assume it is supersonic desired region
        tol     =   0.00000001;
        iter    =   0;
        P0      =   24250000.;
        T0      =   2900.;
    }

    //set public values
    void    set_isentropic(int  a)  ;
        // set isen = 0 for isentropic relationships
        // set isen = 1 for non-isentropic relationships
        // set isen = 2 for Rayleigh Pitot Equation
        // set isen = 3 for Oblique Shock Waves
        // set isen = 4 for Prandtl-Meyer  Equations
    void    set_Rg(double a) ;
    void    set_F();
    void    set_dFdM();
    
    //calculate the next iteration Mach number  in the newton solver iteration
    void Newton_iter();
    // Newton solver for oblique shock waves
    void Newton_iter_oblique();
    void Newton_iter_PM();

    //calculate the error (and set the new F and dFdM values) on each iteration of Mach number in the newton solver
    double error();

    //calculate the converged value of Mach number using the Newton iteration solver
    void Newton();

    //calculate the temperature for an isentropic nozzle
    void Temperature();
    void Pressure();
    // calculate the choked mass flow with the given Gamma, Rg, P0 and T0
    void calc_m_dot();
    // oblique shock wave M to M2
    void calc_M2_oblique();
    // Prandtl-Meyer equation V(M)
    double V(double Mach);
    void Stagnation_Pressure();//isentropic relation to get P0 from P
    void Pressure_Oblique();
    void Temperature_Oblique();
    void P0_Oblique();
    void calc_Oblique();
    void P2_Prandtl_Meyer();
    void T2_Prandtl_Meyer();
    void mu_Prandtl_Meyer();
    void calc_Prandtl_Meyer();
    double Thrust(double conical_nozzle_theta,double Pamb); // Thrust as a function of altitude
};
