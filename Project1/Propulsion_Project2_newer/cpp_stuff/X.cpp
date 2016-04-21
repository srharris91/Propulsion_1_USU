#include "X.h"
#include "DEFS.h"
// calculate the x_dot values using F_thrust and theta
void X::x_dot(double& F_thrust){
    Thrust = F_thrust;
    gamma   = atan(Vr/Vv);
    Vr_dot  = (Vv*Vv/r) - (mu / (r*r))
        + (
                (Thrust/m)
          )
        * sin(gamma);
    Vv_dot  = -((Vr*Vv)/r)
        + (
                (Thrust/m)
          )
        * cos(gamma);
    r_dot   = Vr;
    v_dot   = Vv/r;
    m_dot   = -Thrust/(9.806 * Isp);
}
X X::Trapz(double& F_thrust,double& dt){
    X p,c;    // predictor and correcotr X value initialization
    this->x_dot(F_thrust);    // calculate the current x_dot values
//    double Ei  = -this->mu/(2.*this->a); // initial orbit energy

    // predictor
    p       = *this;
    p.r     = this->r  + dt * this->r_dot   ;
    p.v     = this->v  + dt * this->v_dot   ;
    p.Vr    = this->Vr + dt * this->Vr_dot  ;
    p.Vv    = this->Vv + dt * this->Vv_dot  ;
    p.m     = this->m  + dt * this->m_dot   ;
    p.x_dot(F_thrust)   ;
    // corrector
    c       = *this;
    c.r     = this->r  + dt/2. * (this->r_dot  + p.r_dot)  ;
    c.v     = this->v  + dt/2. * (this->v_dot  + p.v_dot)  ;
    c.Vr    = this->Vr + dt/2. * (this->Vr_dot + p.Vr_dot) ;
    c.Vv    = this->Vv + dt/2. * (this->Vv_dot + p.Vv_dot) ;
    c.m     = this->m  + dt/2. * (this->m_dot  + p.m_dot)  ;
    c.xy();
    c.t     = this->t + dt  ;
//    double Ef  = -c.mu/(2.*c.a); // final orbit energy
    c.DeltaV  = this->DeltaV + 9.806 * c.Isp * log(this->m/c.m);
//    if (abs(Ef-Ei)<1) c.DeltaV = this->DeltaV; // if small DeltaV then do not count DeltaV
    return c;
}
X X::RK4(double& F_thrust, double& dt){
    X k1,k2,k3,k4,c;
    this->x_dot(F_thrust);
//    double Ei  = -mu/(2.*a); // initial orbit energy

    // k1
    k1      = *this;

    // k2
    k2      = k1;
    k2.r    = this->r   + (dt/2.) * k1.r_dot  ;
    k2.v    = this->v   + (dt/2.) * k1.v_dot  ;
    k2.Vr   = this->Vr  + (dt/2.) * k1.Vr_dot ;
    k2.Vv   = this->Vv  + (dt/2.) * k1.Vv_dot ;
    k2.m    = this->m   + (dt/2.) * k1.m_dot  ;
    k2.x_dot(F_thrust);

    // k3
    k3      = k2;
    k3.r    = this->r   + (dt/2.) * k2.r_dot  ;
    k3.v    = this->v   + (dt/2.) * k2.v_dot  ;
    k3.Vr   = this->Vr  + (dt/2.) * k2.Vr_dot ;
    k3.Vv   = this->Vv  + (dt/2.) * k2.Vv_dot ;
    k3.m    = this->m   + (dt/2.) * k2.m_dot  ;
    k3.x_dot(F_thrust);

    // k4
    k4      = k3;
    k4.r    = this->r   + dt    * k3.r_dot  ;
    k4.v    = this->v   + dt    * k3.v_dot  ;
    k4.Vr   = this->Vr  + dt    * k3.Vr_dot ;
    k4.Vv   = this->Vv  + dt    * k3.Vv_dot ;
    k4.m    = this->m   + dt    * k3.m_dot  ;
    k4.x_dot(F_thrust);

    // corrected values
    c       = *this;
    c.r     = this->r   + (dt/6.) * (k1.r_dot  + 2.*k2.r_dot  + 2.*k3.r_dot  + k4.r_dot );
    c.v     = this->v   + (dt/6.) * (k1.v_dot  + 2.*k2.v_dot  + 2.*k3.v_dot  + k4.v_dot );
    c.Vr    = this->Vr  + (dt/6.) * (k1.Vr_dot + 2.*k2.Vr_dot + 2.*k3.Vr_dot + k4.Vr_dot);
    c.Vv    = this->Vv  + (dt/6.) * (k1.Vv_dot + 2.*k2.Vv_dot + 2.*k3.Vv_dot + k4.Vv_dot);
    c.m     = this->m   + (dt/6.) * (k1.m_dot  + 2.*k2.m_dot  + 2.*k3.m_dot  + k4.m_dot );
    c.xy();
    c.t     = this->t + dt  ;
//    double Ef  = -c.mu/(2.*c.a); // final orbit energy
    c.DeltaV  = this->DeltaV + 9.806 * c.Isp * log(this->m/c.m);
//    if (abs(Ef-Ei)<1) c.DeltaV =this->DeltaV; // if small DeltaV then do not count DeltaV

    return c;
}
X X::function(int a, X& Old_f, double& Thrust_f, double& dt_f) {
    X stuff;
    if      (a == 0) {  return Old_f.Trapz  (Thrust_f, dt_f); }
    else if (a == 1) {  return Old_f.RK4    (Thrust_f, dt_f); }
    else                return stuff;
}
double X::HomannTransfer(int a,int b,X start, X coast, X final, X Old, X New){
    double m; // output mass value
    Old = start;
    Old.x_dot(start.Thrust);
    // continuous thrust
    for (int i=0; i<20000000; ++i){
        New=function(a,Old,start.Thrust,start.dt);
        Old=New;    // reset for next iteration
        if (New.a*(1. + New.e)>final.r) break; // if apogee is the max desired radius
    }
    // coast
    for (int i=0; i<20000000; ++i){
        New= function(a,Old,coast.Thrust,coast.dt);
        Old=New;    // reset for next iteration
        if (New.r>=final.r) break;
    }
    if (b == 0) {
        // if impulsive final burn
        double FinalDV = sqrt(final.mu/final.r) - sqrt(2.*(final.mu/(final.r) - New.mu/(2.*New.a)));
        Old.Vv  += FinalDV;
        Old.m   -= final.m - final.m/(exp(FinalDV/(9.806*final.Isp)));
    }
    else if (b==1){
        // final thrust kick
        for (int i=0; i<2000000; ++i){
            New= function(a,Old,final.Thrust,final.dt);
            Old=New;    // reset for next iteration
            if (New.a*(1. - New.e)>=final.r) break; // if perigee is the max desired radius
        }
    }
    // coast after
    double temp_r = New.v+2.*M_PI;
    double error=10.;
    for (int i=0; i<20000; ++i){
        New= function(a,Old,coast.Thrust,coast.dt);
        error = sqrt(pow(New.v - temp_r,2));// rms error
        if (error < 0.01) {break;}
        Old=New;    // reset for next iteration
    }
    m = Old.m;

    return m;
}

double  X::ShootingMethod(int a,int b,X start, X coast, X final, X Old, X New){
    double mass1,mass2,mass_next;
    X start2;
    mass1   = start.HomannTransfer(a,b,start,coast,final,Old,New);
    start2  = start;
    start2.m += 10.; // add 10 kg for start of this shooting method
    mass2 = start2.HomannTransfer(a,b,start2,coast,final,Old,New);
    mass_next = start.m  + ( (start2.m - start.m) / (mass2 - mass1)) * (final.m - start.m);
    cout<<"mass guess = "<<start2.m<<" to get final mass = "<<mass2<<endl;
    cout.precision(14);
    double error = 10.;
    for (; error>4.1;){
        start = start2;
        start2.m = mass_next;
        mass1 = mass2;
        mass2 = start2.HomannTransfer(a,b,start2,coast,final,Old,New);
        mass_next = start.m  + ((start2.m - start.m)/(mass2 - mass1)) * (final.m - mass1);
        error = abs(final.m - mass2);
        cout<<"mass guess = "<<start2.m<<" to get final mass = "<<mass2<<endl;
    }
    return start2.m;
}
void X::HomannTransfer_output(int a,int b,X start, X coast, X final, X Old, X New, ofstream &myfile, ofstream &myfile2){
    stringstream line;
    double error,temp_r;
    Old = start;
    Old.x_dot(start.Thrust);
    Old.output_line(line);
    myfile<<line.str()<<endl;
    // continuous thrust
    for (int i=0; i<200000000; ++i){
        New=function(a,Old,start.Thrust,start.dt);
        New.output_line(line);
        myfile<<line.str()<<endl;
        Old=New;    // reset for next iteration
        if (New.a*(1. + New.e)>final.r) break; // if apogee is the max desired radius
    }
    // coast
    for (int i=0; i<20000000; ++i){
        New= function(a,Old,coast.Thrust,coast.dt);
        New.output_line(line);
        myfile2<<line.str()<<endl;
        Old=New;    // reset for next iteration
        if (New.r>=final.r) break;
    }
    if (b == 0) {
        // if impulsive final burn
        double FinalDV = sqrt(final.mu/final.r) - sqrt(2.*(final.mu/(final.r) - New.mu/(2.*New.a)));
        Old.Vv  += FinalDV;
        Old.m   -= final.m - final.m/(exp(FinalDV/(9.806*final.Isp)));
        cout<<"final time = "<<New.t<<endl;
        cout<<"final mass = "<<Old.m<<endl;
        cout<<"Final DV = "<<FinalDV<<endl;
        cout<<"a and e = "<<New.a<<" "<<New.e<<endl;
        Old.output_line(line);
        myfile<<line.str()<<endl;
        // end if impulsive final burn
    }
    else if (b == 1) {
        // final non-impulsive thrust kick
        for (int i=0; i<20000000; ++i){
            New= function(a,Old,final.Thrust,final.dt);
            New.output_line(line);
            myfile<<line.str()<<endl;
            Old=New;    // reset for next iteration
            if (New.a*(1. - New.e)>=final.r) break; // if perigee is the max desired radius
        }
    }
    // coast after
    temp_r = New.v+2.*M_PI;
    for (int i=0; i<20000; ++i){
        New= function(a,Old,coast.Thrust,coast.dt);
        error = sqrt(pow(New.v - temp_r,2));// rms error
        cout<<'\r'<<"Coasting iteration = "<<i<<" now at "<<New.r<<" "<<New.v;
        if (error < 0.01) {break;}
        New.output_line(line);
        myfile2<<line.str()<<endl;
        Old=New;    // reset for next iteration
    }
    cout<<endl;
    cout<<"Delta V = "<<Old.DeltaV<<endl;
    cout<<"Now V vs start.vv = "<<Old.Vr<<" "<<Old.Vv<<" "<<start.Vv<<endl;
    cout<<"supposed to be =    "<<sqrt(Old.mu/Old.r)<<endl;

}
void X::xy(){
    x = r*cos(v);
    y = r*sin(v);
    a = mu/((2.*mu/r) - (Vr*Vr + Vv*Vv));
    e = (r/mu) * sqrt( pow(Vv*Vv - (mu/r),2) + pow(Vr*Vv,2));

}
void X::output_line(stringstream &a){
    a.str(string());    // clear contents
    a<<fixed<<setprecision(12);
    a<<"  "<<t      ;
    a<<"  "<<x      ;
    a<<"  "<<y      ;
    a<<"  "<<this->a;
    a<<"  "<<e      ;
    a<<"  "<<m      ;
    a<<"  "<<DeltaV ;
}
// print out values
void X::print(){
    cout<<""<<endl;
    cout<<"x      = "<<x      <<endl;
    cout<<"y      = "<<y      <<endl;
    cout<<"a      = "<<a      <<endl;
    cout<<"e      = "<<e      <<endl;
    cout<<"Vr     = "<<Vr     <<endl;
    cout<<"Vv     = "<<Vv     <<endl;
    cout<<"r      = "<<r      <<endl;
    cout<<"v      = "<<v      <<endl;
    cout<<"m      = "<<m      <<endl;
    cout<<"Vr_dot = "<<Vr_dot <<endl;
    cout<<"Vv_dot = "<<Vv_dot <<endl;
    cout<<"r_dot  = "<<r_dot  <<endl;
    cout<<"v_dot  = "<<v_dot  <<endl;
    cout<<"m_dot  = "<<m_dot  <<endl;
    cout<<"gamma  = "<<gamma  <<endl;
    cout<<"Isp    = "<<Isp    <<endl;
    cout<<"mu     = "<<mu     <<endl;
    cout<<"t      = "<<t      <<endl;
    cout<<""<<endl;
}
