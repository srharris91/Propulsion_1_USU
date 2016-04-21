#include "DEFS.h"
#include "X.h"
int main(){
    X start,coast,final,Old,New; // initialize data classes
    stringstream line;

    // read input file
    ifstream input("input",ios::in);
    string LINE,NAME;
    double MASS,RADIUS,THRUST,ISP,DT;
    while ( getline(input,LINE)){
        if (LINE[0] != '#'){
            istringstream ss(LINE) ;
            ss>> NAME >> MASS >> RADIUS >> THRUST >> ISP >> DT;
            if (NAME == "start") {
                start.m     = MASS;
                start.r     = RADIUS;
                start.Thrust= THRUST;
                start.Isp   = ISP;
                start.dt    = DT;
            }
            else if (NAME == "coast"){
                coast.m     = MASS;
                coast.r     = RADIUS;
                coast.Thrust= THRUST;
                coast.Isp   = ISP;
                coast.dt    = DT;
            }
            else if (NAME == "final"){
                final.m     = MASS;
                final.r     = RADIUS;
                final.Thrust= THRUST;
                final.Isp   = ISP;
                final.dt    = DT;
            }
        }
    }
    start.Vv    = sqrt(start.mu/start.r); // initial circular orbit velocity
    start.xy();             // calc xy coordinates from r and nu

/* example
    ofstream myfile,myfile2; // output to file
    myfile.open ("continuous.txt", ios::out | ios::trunc);
    myfile2.open ("coast.txt", ios::out | ios::trunc);
    start.HomannTransfer_output(0,0,start,coast,final,Old,New,myfile,myfile2);
    myfile.close();
    myfile2.close();
*/
/*
    // part a
    cout<<"start.V = "<<start.Vr<<" "<<start.Vv<<endl;
    start.m = start.ShootingMethod(0,0,start,coast,final,Old,New);
    // now make output files
    ofstream myfile,myfile2; // output to file
    myfile.open ("continuous.txt", ios::out | ios::trunc);
    myfile2.open ("coast.txt", ios::out | ios::trunc);
    // output to files
    start.HomannTransfer_output(0,0,start,coast,final,Old,New,myfile,myfile2);
    // close files
    myfile.close();
    myfile2.close();
*/
/*
    // part b
    // impulsive Dv
    double DV1,DV2;
    DV1 = sqrt(start.mu/start.r)*(sqrt(2.*(1.-(1./((final.r/start.r)+1.))))-1.);
    DV2 = sqrt(start.mu/start.r)*(sqrt(start.r/final.r) - sqrt(2.*(start.r/final.r - 1./(final.r/start.r + 1.))));
    start.m   -= final.m - final.m/(exp((DV1+DV2)/(9.806*final.Isp)));
    cout<<"DV1 = "<<DV1<<endl;
    cout<<"DV2 = "<<DV2<<endl;
    cout<<"End mass and propellent = "<<start.m<<" "<< final.m - final.m/(exp((DV1+DV2)/(9.806*final.Isp)))<<endl;
*/
/*
    // part c
    start.m = start.ShootingMethod(1,1,start,coast,final,Old,New);
    // now make output files
    ofstream myfile,myfile2; // output to file
    myfile.open ("continuous.txt", ios::out | ios::trunc);
    myfile2.open ("coast.txt", ios::out | ios::trunc);
    // output to files
    start.HomannTransfer_output(1,1,start,coast,final,Old,New,myfile,myfile2);
    // close files
    myfile.close();
    myfile2.close();
*/
// /*
    // part d
    // custom HomannTransfer equation from X.cpp in order to get customized orbit transfer
    ofstream myfile,myfile2; // output to file
    myfile.open ("continuous.txt", ios::out | ios::trunc);
    myfile2.open ("coast.txt", ios::out | ios::trunc);
    double error,temp_r;
    Old = start;
    Old.x_dot(start.Thrust);
    Old.output_line(line);
    myfile<<line.str()<<endl;
    // continuous thrust
    for (int i=0; i<200000000; ++i){
        New=New.function(1,Old,start.Thrust,start.dt);
        New.output_line(line);
        myfile<<line.str()<<endl;
        Old=New;    // reset for next iteration
        if (New.a*(1. + New.e)>final.r) break; // if apogee is the max desired radius
    }
    // coast
    for (int i=0; i<20000000; ++i){
        New= New.function(1,Old,coast.Thrust,coast.dt);
        New.output_line(line);
        myfile2<<line.str()<<endl;
        Old=New;    // reset for next iteration
        if (New.r>=final.r) break;
        if (New.t>=5400.) break;
    }
    // final non-impulsive thrust kick
    for (int i=0; i<20000000; ++i){
        New= New.function(1,Old,final.Thrust,final.dt);
        New.output_line(line);
        myfile<<line.str()<<endl;
        if (Old.e - New.e < 0.) break;  // break if increasing e
        Old=New;    // reset for next iteration
//        if (New.a*(1. - New.e)>=final.r) break; // if perigee is the max desired radius
    }
    // coast after
    temp_r = New.v+2.*M_PI;
    for (int i=0; i<20000; ++i){
        New= New.function(1,Old,coast.Thrust,coast.dt);
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
// */
    return 0;
}
