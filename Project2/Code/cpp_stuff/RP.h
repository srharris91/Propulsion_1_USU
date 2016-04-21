class RP{
    public:
        double r,r_dot,P0,P0_dot;
        RP();
//        double A_burn(CylindricalPort a);
        double A_burn(CylindricalPort a,int b);
//        double V_c(CylindricalPort a);
        double V_c(CylindricalPort a,int b);
};
