#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <conio.h>
#include <fstream>

using namespace std;
#define PI 3.1415926535
class BikeModel{
    double dt,L,freq;
    double T=-1;
    public:
    BikeModel(double _dt, double _L, double _freq){
        dt=_dt;
        L=_L;
        freq=_freq;
    }
    vector<vector<double>> bike(vector<double> u, vector<double> x0, vector<double> xf){
        vector<vector<double>> X;
        double x_, y_, fi, fi_;
        int k = u.size()/2;
        vector<double> uv, ug;
        copy(u.begin(), u.begin()+k, back_inserter(uv));
        copy(u.begin()+k, u.end(), back_inserter(ug));
        X.push_back(x0);
        T=-1;

        for(int i = 0; i < k; i++){
            double v = uv[i];
            double gam = ug[i];
            for(double t=0; t<1/freq; t+=dt){
                fi=x0[2];
                x_=v*cos(fi);
                y_=v*sin(fi);
                fi_=v/L*tan(gam);
                x0[0]=x0[0]+dt*x_;
                x0[1]=x0[1]+dt*y_;
                x0[2]=x0[2]+dt*fi_;
                X.push_back(x0);
                if (powf(xf[0]-x0[0], 2)+powf(xf[1]-x0[1], 2)+powf(xf[2]-x0[2], 2) < 0.5) {
                    T=t+i/freq;
                    break;
                }
            }
            if (T!=-1) break;
        }
        if (T==-1) T=k/freq;
        
        return X;
    }
};
class MPCOpt{
    double a=0.01,b1=0.9,b2=0.999,e=1e-8;
    int n=10,k=5,l=1;
    int iter_num=1,min_quality=99999;;
    vector<double> x0, xs={5,5,0}, xf={0,0,-PI/2};
    BikeModel model=BikeModel(0.05, 1, 5);
    public:
    vector<vector<double>> mpc_traect(){
        vector<double> tt,tt0,u;
        vector<vector<double>> X, traect;
        for (int i=0; i<n; i++) tt0.push_back(0.1);
        for (int i=0; i<n; i++) tt0.push_back(0);
        copy(xs.begin(),xs.end(),back_inserter(x0));
        tt.resize(tt0.size());
        double dx=0;
        for (int i=0; i<xf.size(); i++) dx+=abs(xf[i]-x0[i]);
        //while(dx>1){
        for (int i=0; i<l; i++){
            tt=mpc_opt(tt0);
            u.erase(u.begin(), u.end());
            for (int j=0; j<k; j++) u.push_back(tt[j]);
            for (int j=n; j<n+k; j++) u.push_back(tt[j]);
            for(double j: u) cout << j << " ";
            cout << endl;
            X=model.bike(u,x0,xf);
            copy(X.begin(), X.end(), back_inserter(traect));
            x0[0]=X[X.size()-1][0];x0[1]=X[X.size()-1][1];x0[2]=X[X.size()-1][2];            
            dx=0;
            for (int j=0; j<xf.size(); j++) dx+=abs(xf[j]-x0[j]);
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
        min_quality=9999;
        //while (min_quality>100){
        for (int i=0; i<iter_num; i++){
            G=grad_bike(tt, tt0);   
            for (int j=0; j<tt.size(); j++){
                m[j]=b1*m[j]+(1-b1)*G[j];
                v[j]=b2*v[j]+(1-b2)*pow(G[j],2);
                if (v_[j]<v[j]) v_[j]=v[j];
                tt[j]=tt[j]-a*m[j]/(sqrt(v_[j])+e);
                tt0[j]=(tt[j]-e);
            }
        }
        return tt;
    }

    vector<double> grad_bike(vector<double> u, vector<double> u_prev){
        vector<double> J_arr, G;
        vector<double> ui;
        double g;
        min_quality=0;
        for (int i=0; i<=u.size(); i++){
            vector<double> ui;
            copy(u.begin(), u.begin()+i, back_inserter(ui));
            copy(u_prev.begin()+i, u_prev.end(), back_inserter(ui));
            J_arr.push_back(J_finder(ui));
        }    
        for(double i: J_arr) {
            min_quality+=i;            
            cout << i << "  ";
        }
        for (int i=0; i<u.size(); i++){
            g=(J_arr[i+1]-J_arr[i])/(u[i]-u_prev[i]);
            if (isnan(g)) g=0;
            G.push_back(g);
        }
        cout << endl;
        return G;
    }

    double J_finder(vector<double> u){
        vector<vector<double>> X=model.bike(u,x0,xf);
        vector<double> x;
        x=X[X.size()-1];
        double J=0;
        // time constrain
        // J+=bike.T*10;
        
        // coordinate constrain
        for(int i=0; i<X.size(); i++){            
            if (X[i][0]>0.5 && X[i][1]<4.7){
                J=J+10*pow((X[i][0]-0.5),2)+10*pow((4.7-X[i][1]),2);
            }
            if (X[i][0]<-0.5 && X[i][1]<4.7){
                J=J+10*pow((0.5+X[i][0]),2)+10*pow((4.7-X[i][1]),2);
            }
        }
        J+=10*pow(xf[0]-x[0], 2)+10*pow(xf[1]-x[1], 2)+25*10*pow(xf[2]-x[2], 2);

        // control constrain
        int k=u.size()/2;
        vector<double> uv, ug;
        copy(u.begin(), u.begin()+k, back_inserter(uv));
        copy(u.begin()+k, u.end(), back_inserter(ug));
        for (int i = 0; i < k; i++){
            //if (abs(uv[i])>5) J+=10*(abs(uv[i])-5);
            if (abs(ug[i])>0.3) J+=10*10*(abs(ug[i])-0.3);
        }
        return J;

    }
};

int main(){
    MPCOpt mpc;
    vector<vector<double>> X = mpc.mpc_traect();
    ofstream fout("cppstudio.txt");
    fout << "X" << "    " << "Y" << "   " << "fi" << endl;
    for (int i=0; i<X.size(); i++){
        fout << X[i][0] << "    " << X[i][1] << "   " << X[i][2] << endl; 
    }
    fout.close();
    return 0;
}