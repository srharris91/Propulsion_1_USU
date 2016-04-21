#include "X.h"
#include "DEFS.h"
// calculate the x_dot values using F_thrust and theta
void X::x_dot(double& F_thrust){
    beta    = m/(C_D*A_ref);
    gamma   = atan(Vr/Vv);
    Vr_dot  = (Vv*Vv/r) - (mu / (r*r))
        + (
                (F_thrust/m)
                - (rho*V_inf*V_inf)/(2.*beta)
          )
        * sin(gamma);
    Vv_dot  = -((Vr*Vv)/r)
        + (
                (F_thrust/m)
                - (rho*V_inf*V_inf)/(2.*beta)
          )
        * cos(gamma);
    r_dot   = Vr;
    v_dot   = Vv/r;
    m_dot   = -F_thrust/(9.81 * Isp);
}
X X::Trapz(double& F_thrust,double& dt){
    X p,c;    // predictor and correcotr X value initialization
    this->x_dot(F_thrust);    // calculate the current x_dot values

    // predictor
    p       = *this;
    p.r     = this->r  + dt * this->r_dot   ;
    p.v     = this->v  + dt * this->v_dot   ;
    p.Vr    = this->Vr + dt * this->Vr_dot  ;
    p.Vv    = this->Vv + dt * this->Vv_dot  ;
    p.m     = this->m  + dt * this->m_dot   ;
    p.x_dot(F_thrust)   ;
    // corrector
    c       = p;
    c.r     = this->r  + dt/2. * (this->r_dot  + p.r_dot)  ;
    c.v     = this->v  + dt/2. * (this->v_dot  + p.v_dot)  ;
    c.Vr    = this->Vr + dt/2. * (this->Vr_dot + p.Vr_dot) ;
    c.Vv    = this->Vv + dt/2. * (this->Vv_dot + p.Vv_dot) ;
    c.m     = this->m  + dt/2. * (this->m_dot  + p.m_dot)  ;
    c.xy();
    c.t     = this->t + dt  ;
    return c;
}
void X::x_set(X& a){
    Vr      = a.Vr      ;
    Vv      = a.Vv      ;
    r       = a.r       ;
    v       = a.v       ;
    m       = a.m       ;
    Vr_dot  = a.Vr_dot  ;
    Vv_dot  = a.Vv_dot  ;
    r_dot   = a.r_dot   ;
    v_dot   = a.v_dot   ;
    m_dot   = a.m_dot   ;
    beta    = a.beta    ;
    gamma   = a.gamma   ;
    C_D     = a.C_D     ;
    A_ref   = a.A_ref   ;
    Isp     = a.Isp     ;
    mu      = a.mu      ;
    t       = a.t       ;
}
void X::xy(){
    x = r*cos(v);
    y = r*sin(v);

}
// print out values
void X::print(){
    cout<<""<<endl;
    cout<<"x      = "<<x      <<endl;
    cout<<"y      = "<<y      <<endl;
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
    cout<<"beta   = "<<beta   <<endl;
    cout<<"gamma  = "<<gamma  <<endl;
    cout<<"C_D    = "<<C_D    <<endl;
    cout<<"A_ref  = "<<A_ref  <<endl;
    cout<<"Isp    = "<<Isp    <<endl;
    cout<<"mu     = "<<mu     <<endl;
    cout<<"t      = "<<t      <<endl;
    cout<<""<<endl;
}
void X::output_line(stringstream &a){
    a.str(string());    // clear contents
    a<<fixed<<setprecision(14);
    a<<"  "<<t      ;
    a<<"  "<<x      ;
    a<<"  "<<y      ;
    a<<"  "<<Vr     ;
    a<<"  "<<Vv     ;
    a<<"  "<<r      ;
    a<<"  "<<v      ;
    a<<"  "<<m      ;
    a<<"  "<<Vr_dot ;
    a<<"  "<<Vv_dot ;
    a<<"  "<<r_dot  ;
    a<<"  "<<v_dot  ;
    a<<"  "<<m_dot  ;
    a<<"  "<<beta   ;
    a<<"  "<<gamma  ;
    a<<"  "<<C_D    ;
    a<<"  "<<A_ref  ;
    a<<"  "<<Isp    ;
    a<<"  "<<mu     ;
//    return a;
}
