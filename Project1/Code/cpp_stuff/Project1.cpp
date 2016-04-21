#include "DEFS.h"
#include "X.h"

// assume no drag
int main(){

    vector<X> Data;
    X start;
    stringstream line;
    double dt = 0.01;
    // part a parameters
    double Thrust = 10.;    // Newtons
    start.Isp   = 2000.;    // seconds
    start.C_D   = 0.0000000000001;  // neglect atmospheric drag
    start.A_ref = 0.0000000000001;  // still neglecting atmospheric drag...
    start.rho   = 0.0000000000001;  // still neglecting atmospheric drag
    start.V_inf = 0.0000000000001;  // still neglecting atmospheric drag
    start.mu    = 3.9860044E9;      // mu for earth
    start.r     = 8530000.; // meters
    start.v     = 0.;       // start nu 
    start.xy();             // calc xy coordinates from r and nu
    start.Vr    = 0.;       // circular orbit
    start.Vv    = sqrt(start.mu * start.r)/start.r; // 
    start.m     = 10000.;    // kg
    start.t     = 0.;
    Data.push_back(start);
//    Data.push_back(start);
//    Data[1] = start.Trapz(Thrust,dt); // first step
//    Data[0].output_line(line);
//    cout<<line.str()<<endl;
//    Data[1].output_line(line);
//    cout<<endl;
//    cout<<line.str()<<endl;
//    cout<<endl;
//    cout<<Data.size()<<endl;
    ofstream myfile; // output to file
    myfile.open ("output_file.txt", ios::out | ios::trunc);

    for (int i=0; i<20000000; ++i){
        Data.push_back(Data[i].Trapz(Thrust,dt));
        Data[i+1].output_line(line);
        myfile<<line.str()<<endl;
    }


        myfile.close();
    




    return 0;
}
