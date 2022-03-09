#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <conio.h>
#include <fstream>
#include <ctime>  

using namespace std;
#define PI 3.1415926535
enum WHEEL{
    FL,
    FR,
    RL,
    RR,
    F,
    R
};
enum WHEEL_SPEED{
    XFL,
    YFL,
    XFR,
    YFR,
    XRL,
    YRL,
    XRR,
    YRR
};
enum SS{
    Vx_,
    Vy_,
    w_,
    Vx,
    Vy,
    w,
    X,
    Y,
    fi,
    gamfl,
    gamfr,
    gamrl,
    gamrr,
    MFL,
    MFR,
    MRL,
    MRR,
    UGFL,
    UGFR,
    UGRL,
    UGRR
};
double sign(double a){
    if (a > 0) return 1;
    else if (a < 0) return 0;
    else return 0;
}
class CarDynModel{
    double L=0.94, C=0.63, H=0.2;
    double lf=L/2, lr=L/2, dr=C/2, dl=C/2;
    double r=0.25/2, m=70, J=228, g=9.8;
    double mu=0.8, P=m*g;
    double K=1, TauG=0.1;
    double Dax=0, Day=0;
    double long_angle=0, lateral_angle=0;
    double dt,freq;
    double T=-1;
    public:
    CarDynModel(double _dt, double _freq){
        dt=_dt;
        freq=_freq;
    }
    vector<vector<double>> car(vector<double> u, vector<double> x0, vector<double> xf){
        vector<vector<double>> trec;
        double P_arr[6], Pa;
        double v_speed[8], v_acc[4];
        vector<double> M, ug;
        double Fx[4], Fy[4], Fy_[4];
        double alfa[4], asl[4], tau[4], lam[4], Cp[4];
        double Mz[4];
        int k = u.size()/2;
        trec.push_back(x0);
        
        T=-1;

        for(int i = 0; i < k/4; i++){
            M.clear(); ug.clear();
            copy(u.begin()+i*4, u.begin()+(i+1)*4, back_inserter(M));
            copy(u.begin()+k+i*4, u.begin()+k+(i+1)*4, back_inserter(ug));
            
            for(double t=0; t<1/freq; t+=dt){ 
                x0=trec[trec.size()-1];   
                P_arr[F]=((P*lr)*cos(long_angle)-P/g*x0[Vx_]*H-Dax*H-P*H*sin(long_angle))/L;
                P_arr[R]=((P*lf)*cos(long_angle)+P/g*x0[Vx_]*H+Dax*H+P*H*sin(long_angle))/L;

                P_arr[FL]=((P_arr[F]*dr)*cos(lateral_angle)-P_arr[F]/g*x0[Vy_]*H-Day*H-P_arr[F]*H*sin(lateral_angle))/C;
                P_arr[FR]=((P_arr[F]*dl)*cos(lateral_angle)+P_arr[F]/g*x0[Vy_]*H+Day*H+P_arr[F]*H*sin(lateral_angle))/C;
                P_arr[RL]=((P_arr[R]*dr)*cos(lateral_angle)-P_arr[R]/g*x0[Vy_]*H-Day*H-P_arr[R]*H*sin(lateral_angle))/C;
                P_arr[RR]=((P_arr[R]*dl)*cos(lateral_angle)+P_arr[R]/g*x0[Vy_]*H+Day*H+P_arr[R]*H*sin(lateral_angle))/C;
                Pa=P_arr[FL]+P_arr[FR]+P_arr[RL]+P_arr[RR];
                
                v_speed[XFL]=(x0[Vx]-dl*x0[w])*cos(x0[gamfl])+(x0[Vy]+lf*x0[w])*sin(x0[gamfl]);
                v_speed[YFL]=-(x0[Vx]-dl*x0[w])*sin(x0[gamfl])+(x0[Vy]+lf*x0[w])*cos(x0[gamfl]);
                v_speed[XFR]=(x0[Vx]+dr*x0[w])*cos(x0[gamfr])+(x0[Vy]+lf*x0[w])*sin(x0[gamfr]);
                v_speed[YFR]=-(x0[Vx]+dr*x0[w])*sin(x0[gamfr])+(x0[Vy]+lf*x0[w])*cos(x0[gamfr]);
                v_speed[XRL]=(x0[Vx]-dl*x0[w])*cos(x0[gamrl])+(x0[Vy]-lr*x0[w])*sin(x0[gamrl]);
                v_speed[YRL]=-(x0[Vx]-dl*x0[w])*sin(x0[gamrl])+(x0[Vy]-lr*x0[w])*cos(x0[gamrl]);
                v_speed[XRR]=(x0[Vx]+dr*x0[w])*cos(x0[gamrr])+(x0[Vy]-lr*x0[w])*sin(x0[gamrr]);
                v_speed[YRR]=-(x0[Vx]+dr*x0[w])*sin(x0[gamrr])+(x0[Vy]-lr*x0[w])*cos(x0[gamrr]);

                v_acc[FL]=-(x0[Vx_]-dl*x0[w_])*sin(x0[gamfl])+(x0[Vy_]+lf*x0[w_])*cos(x0[gamfl]);
                v_acc[FR]=-(x0[Vx_]+dr*x0[w_])*sin(x0[gamfr])+(x0[Vy_]+lf*x0[w_])*cos(x0[gamfr]);
                v_acc[RL]=-(x0[Vx_]-dl*x0[w_])*sin(x0[gamrl])+(x0[Vy_]-lr*x0[w_])*cos(x0[gamrl]);
                v_acc[RR]=-(x0[Vx_]+dr*x0[w_])*sin(x0[gamrr])+(x0[Vy_]-lr*x0[w_])*cos(x0[gamrr]);

                for (uint8_t j=0; j<4; j++){
                    Fx[j]=M[j]/r;
                    alfa[j]=atan2(v_speed[j*2+1], -v_speed[j*2]);
                    Fy_[j]=v_acc[j]*m*P_arr[j]/Pa;
                    if (v_speed[j*2+1] == 0)
                        Cp[j]=PI/2;
                    else
                        Cp[j]=abs(atan(Fy_[j]/(-v_speed[j*2+1])*v_speed[j*2]));
                    tau[j]=Cp[j]/3/mu/P_arr[j];
                    asl[j]=atan(tau[j]);
                    lam[j]=1-tau[j]*abs(tan(alfa[j]));
                    if (abs(alfa[j])<=asl[j]){
                        Fy[j]=mu*P_arr[j]*(1-pow(lam[j],3))*sign(alfa[j]);
                        Mz[j]=-mu*P_arr[j]*pow(lam[j],3)*sqrt(P_arr[j]/3/100000)*(1-lam[j])*sign(alfa[j]);
                    }
                    else{
                        Fy[j]=-mu*P_arr[j]*sign(alfa[j]);
                        Mz[j]=0;
                    }
                }
                x0[Vx_]=(Fx[FR]*cos(x0[gamfr])+Fx[RR]*cos(x0[gamrr])+Fx[FL]*cos(x0[gamfl])+Fx[RL]*cos(x0[gamrl]) 
                    +Fy[RR]*cos(x0[gamrr]+PI/2)+Fy[FR]*cos(x0[gamfr]+PI/2)+Fy[RL]*cos(x0[gamrl]+PI/2)+Fy[FL]*cos(x0[gamfl]+PI/2) 
                    -P*sin(long_angle))/m-x0[w]*x0[Vy];
                x0[Vy_]=(Fx[RR]*cos(x0[gamrr]-PI/2)+Fx[FR]*cos(x0[gamfr]-PI/2)+Fx[RL]*cos(x0[gamrl]-PI/2)+Fx[FL]*cos(x0[gamfl]-PI/2) 
                    +Fy[RR]*cos(x0[gamrr])+Fy[FR]*cos(x0[gamfr])+Fy[RL]*cos(x0[gamrl])+Fy[FL]*cos(x0[gamfl])-P*sin(lateral_angle))/m+x0[w]*x0[Vx];
                x0[w_]=(Fx[RL]*(-lr*sin(x0[gamrl])-dl*cos(x0[gamrl]))+Fy[RL]*(-lr*cos(x0[gamrl])+dl*sin(x0[gamrl])) 
                    +Fx[RR]*(-lr*sin(x0[gamrr])+dr*cos(x0[gamrr]))+Fy[RR]*(-lr*cos(x0[gamrr])-dr*sin(x0[gamrr])) 
                    +Fx[FL]*(lf*sin(x0[gamfl])-dl*cos(x0[gamfl]))+Fy[FL]*(lf*cos(x0[gamfl])+dl*sin(x0[gamfl])) 
                    +Fx[FR]*(lf*sin(x0[gamfr])+dr*cos(x0[gamfr]))+Fy[FR]*(lf*cos(x0[gamfr])-dr*sin(x0[gamfr])) 
                    +(+Mz[FL]+Mz[FR]+Mz[RL]+Mz[RR]))/J;
                x0[Vx]=x0[Vx]+dt*x0[Vx_];
                x0[Vy]=x0[Vy]+dt*x0[Vy_];
                x0[w]=x0[w]+dt*x0[w_];
                x0[X]=x0[X]+dt*(cos(x0[fi])*x0[Vx]-sin(x0[fi])*x0[Vy]);
                x0[Y]=x0[Y]+dt*(sin(x0[fi])*x0[Vx]+cos(x0[fi])*x0[Vy]);
                x0[fi]=x0[fi]+dt*x0[w];
                x0[gamfl]=x0[gamfl]+dt*(-(1/TauG)*x0[gamfl]+(K/TauG)*ug[FL]);
                x0[gamfr]=x0[gamfr]+dt*(-(1/TauG)*x0[gamfr]+(K/TauG)*ug[FR]);
                x0[gamrl]=x0[gamrl]+dt*(-(1/TauG)*x0[gamrl]+(K/TauG)*ug[RL]);
                x0[gamrr]=x0[gamrr]+dt*(-(1/TauG)*x0[gamrr]+(K/TauG)*ug[RR]);
                x0[MFL]=M[FL];
                x0[MFR]=M[FR];
                x0[MRL]=M[RL];
                x0[MRR]=M[RR];
                x0[UGFL]=ug[FL];
                x0[UGFR]=ug[FR];
                x0[UGRL]=ug[RL];
                x0[UGRR]=ug[RR];
                trec.push_back(x0);
                // if (powf(xf[0]-x0[X], 2)+powf(xf[1]-x0[Y], 2)+powf(xf[2]-x0[fi], 2) < 0.5) {
                //     T=t+i/freq;
                //     break;
                // }
            }
            
            // cout << "M: " << M[0] << "  " << M[1] << "  " << M[2] << "  " << M[3] << "  " << endl; 
            // cout << "ug: " << ug[0] << "  " << ug[1] << "  " << ug[2] << "  " << ug[3] << "  " << endl; 
            // cout << "Fx: " << Fx[0] << "  " << Fx[1] << "  " << Fx[2] << "  " << Fx[3] << "  " << endl; 
            // cout << "Fy: " << Fy[0] << "  " << Fy[1] << "  " << Fy[2] << "  " << Fy[3] << "  " << endl; 
            // if (T!=-1) break;
        }
        // if (T==-1) T=k/freq;
        
        return trec;
    }
};

class MPCOpt{
    double a=0.03,b1=0.9,b2=0.999,e=1e-8;
    int n=5,k=5,l=1;
    int iter_num=20000,min_quality=99999;
    vector<double> x0, xs, xf={0,0,PI/2};
    vector<vector<double>> obst;
    CarDynModel c=CarDynModel(0.05, 2);
    public:
    vector<vector<double>> mpc_traect(){
        vector<double> tt,tt0,u;
        vector<vector<double>> loc_traect, traect;
        for (uint8_t i=Vx_; i<=UGRR; i++) xs.push_back(0);
        xs[X] = 5;
        xs[Y] = 5;
        vector<double> obst_i={0,0};
        for (float i=1; i<=5; i+=0.5){
            obst_i[0] = i;
            obst_i[1] = 4;
            obst.push_back(obst_i);
        }
        for (float i=1; i<=4; i+=0.5){
            obst_i[0] = 1;
            obst_i[1] = i;
            obst.push_back(obst_i);
        }
        
        for (float i=-1; i>=-5; i-=0.5){
            obst_i[0] = i;
            obst_i[1] = 4;
            obst.push_back(obst_i);
        }
        for (float i=1; i<=4; i+=0.5){
            obst_i[0] = -1;
            obst_i[1] = i;
            obst.push_back(obst_i);
        }

        for (int i=0; i<n*4; i++) tt0.push_back(-10);
        for (int i=0; i<n*4; i++) tt0.push_back(0.3);
        copy(xs.begin(),xs.end(),back_inserter(x0));
        tt.resize(tt0.size());
        for (int i=0; i<l; i++){
            tt=mpc_opt(tt0);
            u.clear();
            for (int j=0; j<k*4; j++) u.push_back(tt[j]);
            for (int j=n*4; j<n*4+k*4; j++) u.push_back(tt[j]);
            // for(double j: u) cout << j << " ";
            // cout << endl;
            loc_traect=c.car(u,x0,xf);
            copy(loc_traect.begin(), loc_traect.end(), back_inserter(traect));
            x0=loc_traect[loc_traect.size()-1];            
        }
        return traect;
    }

    vector<double> mpc_opt(vector<double> tt){
        vector<double> tt0;        
        vector<double> G;
        vector<double> m,v,v_;
        for (double i: tt){
           tt0.push_back(i-e); 
           m.push_back(0);
           v.push_back(0);
           v_.push_back(0);
        } 

        unsigned int start_time; // начальное время
        unsigned int end_time; // конечное время
        unsigned int search_time;

        min_quality=9999;
        //while (min_quality>100){
        for (int i=0; i<iter_num; i++){
            start_time = clock();

            G=grad_bike(tt, tt0);   
            for (int j=0; j<tt.size(); j++){
                m[j]=b1*m[j]+(1-b1)*G[j];
                v[j]=b2*v[j]+(1-b2)*pow(G[j],2);
                if (v_[j]<v[j]) v_[j]=v[j];
                tt[j]=tt0[j]-a*m[j]/(sqrt(v_[j])+e);
                tt0[j]=tt[j]-e;
                // cout << tt0[j] << " ";
            }
            // cout << "tt: ";
            // for (int j=0; j<tt.size(); j++){cout << tt0[j] << " ";}
            // cout << endl;
            // cout << "G: ";
            // for (int j=0; j<tt.size(); j++){cout << tt0[j] << " ";}
            // cout << endl;

            end_time = clock(); // конечное время
            search_time = end_time - start_time; // искомое время
            cout << "[" << i << "]" << " time  = " << search_time/1000.0 << endl;

        }
        return tt;
    }

    vector<double> grad_bike(vector<double> u, vector<double> u_prev){
        vector<double> J_arr, G;
        vector<double> ui;
        double g;
        min_quality=0;
        for (int i=0; i<=u.size(); i++){
            ui.clear();
            copy(u.begin(), u.begin()+i, back_inserter(ui));
            copy(u_prev.begin()+i, u_prev.end(), back_inserter(ui));
            J_arr.push_back(J_finder(ui));
        }    
        for(double i: J_arr) {
            min_quality+=i;    
        }
        for (int i=0; i<u.size(); i++){
            g=(J_arr[i+1]-J_arr[i])/(e);
            if (isnan(g)) g=0;
            G.push_back(g);
        }
        return G;
    }

    double J_finder(vector<double> u){
        vector<vector<double>> loc_traect=c.car(u,x0,xf);
        vector<double> x;
        x=loc_traect[loc_traect.size()-1];
        double J=0;
        // time constrain
        // J+=bike.T*10;

        // coordinate constrain
        double cond;
        for(int i=0; i<loc_traect.size(); i++){    
            for (int j=0; j<obst.size(); j++){
                cond = sqrt((obst[j][0]-loc_traect[i][X])*(obst[j][0]-loc_traect[i][X]) + (obst[j][1]-loc_traect[i][Y])*(obst[j][1]-loc_traect[i][Y]));    
                if (cond<0.5) 
                    J+=0.5/cond*20-20;    
            }   
        }
        J+=10*pow(xf[0]-x[X], 2)+10*pow(xf[1]-x[Y], 2);//+3*10*pow(xf[2]-x[fi], 2);

        // control constrain
        // int k=u.size()/2;
        // vector<double> uv, ug;
        // copy(u.begin(), u.begin()+k, back_inserter(uv));
        // copy(u.begin()+k, u.end(), back_inserter(ug));
        // for (int i = 0; i < k; i++){
        //     if (abs(uv[i])>5) J+=10*(abs(uv[i])-5);
        //     if (abs(ug[i])>1.5) J+=10*10*(abs(ug[i])-1.5);
        // }
        return J;

    }
};

int main(){
    MPCOpt mpc;
    unsigned int start_time =  clock(); // начальное время
    cout << "### START COMPUTING ###" << endl;
    vector<vector<double>> traect = mpc.mpc_traect();

    unsigned int end_time = clock(); // конечное время
    unsigned int search_time = end_time - start_time; // искомое время
    cout << " computation time  = " << search_time/1000.0 << endl;
    cout << "### STOP COMPUTING ###" << endl;
    
    // CarDynModel c=CarDynModel(0.001, 1);
    // vector<double> x0, u, xf={5,5,0};
    // u={7,7,7,7,1.5,1.5,1.5,1.5};
    // for (uint8_t i=Vx_; i<=UGRR; i++) x0.push_back(0);
    // vector<vector<double>> traect = c.car(u,x0,xf);
    ofstream fout("dyn1.txt");
    fout << "X  Y  U1  U2  U3  U4  UG1  UG2  UG3  UG4" << endl;
    for (vector<double> i : traect) 
        fout << i[X] << "  " << i[Y] << "  " << i[MFL] << "  " << i[MFR] << "  " << i[MRL] << "  " << i[MRR] <<  
            "  " << i[UGFL] << "  " << i[UGFR] << "  " << i[UGRL] << "  " << i[UGRR] << endl;
    fout.close();
    return 0;
}