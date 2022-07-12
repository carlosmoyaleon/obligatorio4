#include <iostream>
#include<cmath>
#include<fstream>
#include <algorithm>
#include <complex>



using namespace std;




#define G 6.67*pow(10,-11)
#define Mt 5.9736*pow(10,24)
#define M_l 0.07349*pow(10,24)
#define dtl 3.844*pow(10,8)
#define w 2.6617*pow(10,-6)
#define Rt 6.378*pow(10,6)
#define Rl 1.7374*pow(10,6)
#define pasos 250000
#define m 10000
#define angulo 0.26


double CalculoR(double r_cohete, double phi_c, double pphi, double t);

double CalculoPHI(double r_cohete, double phi_c, double t);


int main( ){

    //Declaramos las variables

    double r_luna, phi_luna, r_cohete, phi_c, p_r, pphi, h=0.8, k[4][4],v;
    int i,j;
    ofstream datos;

    //Inicializamos los parametro


    v=11194.444;
    r_luna=dtl;
    phi_luna=0;
    r_cohete=Rt;
    phi_c=0.;
     phi_luna=0;
    

    //Reescalamos los valores

    r_luna=r_luna/(1*dtl);
    r_cohete=Rt/(1*dtl);
    v=v/(1*dtl);
    p_r=v*cos(angulo);
    pphi=v*r_cohete*sin(angulo);  


    datos.open("planets_data.dat");

    for ( i = 0; i < pasos; i++)
    {
    

        if((i%500)==0){

        datos << r_cohete*cos(phi_c) << ", " << r_cohete*sin(phi_c) << endl << r_luna*cos(w*i*h*1.) << ", " << r_luna*sin(w*i*h*1.) << endl << endl;


        }    


        //Evaluamos los k 

        k[0][0]=h*p_r;
        k[0][1]=h*pphi/(r_cohete*r_cohete);
        k[0][2]=h*CalculoR(r_cohete, phi_c, pphi, i*h);
        k[0][3]=h*CalculoPHI(r_cohete, phi_c, i*h);

        k[1][0]=h*(p_r+k[0][2]/2);
        k[1][1]=h*(pphi+k[0][3]/2)/(((r_cohete+k[0][0]/2)*(r_cohete+k[0][0]/2)));
        k[1][2]=h*CalculoR(r_cohete+ k[0][2]/2, phi_c + k[0][1]/2, pphi+ k[0][3]/2 , i*h);
        k[1][3]=h*CalculoPHI(r_cohete+ k[0][0]/2, phi_c+k[0][1]/2, i*h);

        k[2][0]=h*(p_r+k[1][2]/2);
        k[2][1]=h*(pphi+k[1][3]/2)/(((r_cohete+k[1][0]/2)*(r_cohete+k[1][0]/2)));
        k[2][2]=h*CalculoR(r_cohete+ k[1][2]/2, phi_c + k[1][1]/2, pphi+ k[1][3]/2 , i*h);
        k[2][3]=h*CalculoPHI(r_cohete+ k[1][0]/2, phi_c+k[1][1]/2, i*h);

        k[3][0]=h*(p_r+k[2][2]); 
        k[3][1]=h*(pphi+k[2][3])/(((r_cohete+k[2][0])*(r_cohete+k[2][0])));
        k[3][2]=h*CalculoR(r_cohete+ k[2][2], phi_c + k[2][1], pphi+ k[2][3] , i*h);
        k[3][3]=h*CalculoPHI(r_cohete+ k[2][0], phi_c+k[2][1], i*h);

        //Ahora calcularemos los valores de phi y r

        pphi= pphi + (k[0][3]+2*k[1][3]+2*k[2][3]+k[3][3])/6;
        r_cohete=r_cohete + (k[0][0]+2*k[1][0]+2*k[2][0]+k[3][0])/6;
        phi_c= phi_c + (k[0][1]+2*k[1][1]+2*k[2][1]+k[3][1])/6;
        p_r= p_r + (k[0][2]+2*k[1][2]+2*k[2][2]+k[3][2])/6;


    }
    

    datos.close();


    return 0;
}



double CalculoPHI(double r_cohete, double phi_c, double t){


    double ppunto_phi,DELTA, mu,rprima;

    DELTA=G*Mt/(1.*dtl*dtl*dtl);
    mu=M_l/(1.*Mt);
    rprima=sqrt(1+r_cohete*r_cohete-2*r_cohete*cos(phi_c-1.*w*t));

    ppunto_phi=-DELTA*mu*r_cohete*sin(phi_c-1.*w*t)/(rprima*rprima*rprima);



    return ppunto_phi;
}


double CalculoR(double r_cohete, double phi_c, double pphi, double t){

    double ppunto_r,DELTA, mu,rprima;

    DELTA=G*Mt/(1.*dtl*dtl*dtl);
    mu=M_l/(1.*Mt);
    rprima=sqrt(1+r_cohete*r_cohete-2*r_cohete*cos(phi_c-1.*w*t));

    ppunto_r=pphi*pphi/(r_cohete*r_cohete*r_cohete)-DELTA*(1/(r_cohete*r_cohete)+mu/(rprima*rprima*rprima)*(r_cohete-cos(phi_c-1.*w*t)));



    return ppunto_r;
}