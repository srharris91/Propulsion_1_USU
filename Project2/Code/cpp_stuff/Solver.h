class Solver{
    public:
        CylindricalPort CP;
        RP RP0,RP1;
        Prandtl_Meyer Port,Nozzle;
        double t0,t1;
        int bates; // if not bates grain then bates=0, bates grain otherwise
        // constructor
        Solver();
        Solver(CylindricalPort a, RP b,Prandtl_Meyer c, Prandtl_Meyer d);
        
        // solve P_dot and R_dot
        void calc_RP1_dot();
        // Solve P1 and R1 values
        void calc_RP1();
        // write out values
        void writeout();
        // reset values
        void reset();
        // next step
        void next_step();
        


};
