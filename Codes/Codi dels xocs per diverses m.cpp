#include<iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <list>
#include <vector>
#include <fstream>

using namespace std;

double Fx(double x, double y, double z, double Vy, double Vz);
double fx(double Vx);
double Fy(double x, double y, double z, double Vx, double Vz);
double fy(double Vy);
double Fz(double x, double y, double z, double Vx, double Vy);
double fz(double Vz);

int main(){
	
	// Massa (kg)
	double m = 1.67*pow(10,-27);
	
	// Càrrega (c)
	double q = 1.6*pow(10,-19);
	
	// Radi de la nau (m) //
	double r = 5;
	
	// Longitud de la nau (m) //
	double l = 50;
	
	// Moment dipolar magnètic de l'iman (A*m^2) //
	double M = pow(10,6);
	
	// Permeabilitat magnètica del buit
	double mu = 4*M_PI*pow(10,-7);
	
	// Velocitat partícules carregades (m/s) //
	double Vp = 550000;
	
	// Velocitat de la nau (m/s) //
	double Vn = 11000;
	
	// Temps característic (s) //
	double t_0 = pow(10,4)/561000;
	
	int N_particules = 100;
	double dt = pow(10,-7);
	int n = 100000;
	int n_imprimir=400;
	int h;
	int j;
	int s;
	ofstream fitxer("I vs xocs.txt");
	cout  << "m(m^2A)" <<"________" << "Xocs(%)" << endl;
	for(s=1; s<=1500; s++){
		
		M = M + 10000;
		
	    // Longitud característica (m)
	   double l_0 = pow(((mu*q*M*t_0)/(4*M_PI*m)),1.0/3.0);
	   double Impactes = 0;
	   double min = 100;
	   for(h=1; h<=N_particules; h++){
	
	      // Posició inicial
	      double x_0=-90, y_0=0, z_0=h*r/N_particules;
	
	      // Normalització
	      x_0 = x_0/l_0;
	      y_0 = y_0/l_0;
	      z_0 = z_0/l_0;

	      double x = x_0;
          double y = y_0;
          double z = z_0;
       
          // Velocitat inicial
	      double Vx_0=Vp+Vn, Vy_0=0, Vz_0=0;

	      // Normalització
	      Vx_0 = Vx_0*t_0/l_0;
	      Vy_0 = Vy_0*t_0/l_0;
	      Vz_0 = Vz_0*t_0/l_0;
	
	      double Vx = Vx_0;
	      double Vy = Vy_0;
	      double Vz = Vz_0;
	      
          double Senyal = 1;
          for(j=0; j<=n; j++){
		
		   double K1X = fx(Vx);
		   double K1Y = fy(Vy);
		   double K1Z = fz(Vz);
		
	  	   double L1X = Fx(x,y,z,Vy,Vz);
		   double L1Y = Fy(x,y,z,Vx,Vz);
		   double L1Z = Fz(x,y,z,Vx,Vy);
		
		   double K2X = fx(Vx + dt*L1X/2);
		   double K2Y = fy(Vy + dt*L1Y/2);
		   double K2Z = fz(Vz + dt*L1Z/2);
		
		   double L2X = Fx(x + dt*K1X/2, y + dt*K1Y/2, z + dt*K1Z/2, Vy + dt*L1Y/2, Vz + dt*L1Z/2);
		   double L2Y = Fy(x + dt*K1X/2, y + dt*K1Y/2, z + dt*K1Z/2, Vx + dt*L1X/2, Vz + dt*L1Z/2);
		   double L2Z = Fz(x + dt*K1X/2, y + dt*K1Y/2, z + dt*K1Z/2, Vx + dt*L1X/2, Vy + dt*L1Y/2);
		
		   double K3X = fx(Vx + dt*L2X/2);
		   double K3Y = fy(Vy + dt*L2Y/2);
		   double K3Z = fz(Vz + dt*L2Z/2);
		
		   double L3X = Fx(x + dt*K2X/2, y + dt*K2Y/2, z + dt*K2Z/2, Vy + dt*L2Y/2, Vz + dt*L2Z/2);
	       double L3Y = Fy(x + dt*K2X/2, y + dt*K2Y/2, z + dt*K2Z/2, Vx + dt*L2X/2, Vz + dt*L2Z/2);
	 	   double L3Z = Fz(x + dt*K2X/2, y + dt*K2Y/2, z + dt*K2Z/2, Vx + dt*L2X/2, Vy + dt*L2Y/2);
		
		   double K4X = fx(Vx + dt*L3X);
		   double K4Y = fy(Vy + dt*L3Y);
		   double K4Z = fz(Vz + dt*L3Z);
		
		   double L4X = Fx(x+ dt*K3X, y+ dt*K3Y, z+ dt*K3Z, Vy + dt*L3Y, Vz + dt*L3Z);
		   double L4Y = Fy(x+ dt*K3X, y+ dt*K3Y, z+ dt*K3Z, Vx + dt*L3X, Vz + dt*L3Z);
		   double L4Z = Fz(x+ dt*K3X, y+ dt*K3Y, z+ dt*K3Z, Vx + dt*L3X, Vy + dt*L3Y);
		
		   x = x + dt*(K1X + 2*K2X + 2*K3X + K4X)/6;
		   y = y + dt*(K1Y + 2*K2Y + 2*K3Y + K4Y)/6;
		   z = z + dt*(K1Z + 2*K2Z + 2*K3Z + K4Z)/6;
		
		   Vx = Vx + dt*(L1X + 2*L2X + 2*L3X + L4X)/6;
		   Vy = Vy + dt*(L1Y + 2*L2Y + 2*L3Y + L4Y)/6;
		   Vz = Vz + dt*(L1Z + 2*L2Z + 2*L3Z + L4Z)/6;

           if(sqrt(pow(x,2)+pow(y,2)+pow(z,2)) < min){
           	min = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
           }
           
           if(sqrt(pow(z,2)+pow(y,2)) < (r/l_0) && x < l/(2*l_0) && x > -l/(2*l_0) && Senyal==1){
			Impactes = Impactes + 1;
			Senyal = 0;
		   }
		   
		  }
	
	   }
	   cout<< M;
	   cout<< "________";
	   cout<< Impactes << endl;
	   fitxer << M <<"\t" << Impactes << "\n";
	}      
	fitxer.close();
}


double Fx(double x, double y, double z, double Vy, double Vz){
	double R = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
	double Valor = Vy*((3*pow(z,2)/(pow(R,5))) - 1/(pow(R,3))) - Vz*3*z*y/(pow(R,5));
	return Valor;	
	}
	
double fx(double Vx){
	return Vx;
}

double Fy(double x, double y, double z, double Vx, double Vz){
	double R = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
	double Valor = (Vz*3*z*x/(pow(R,5))) - Vx*(3*pow(R,2)/(pow(R,5)) - 1/(pow(R,3)));
	return Valor;
}

double fy(double Vy){
	return Vy;
}

double Fz(double x, double y, double z, double Vx, double Vy){
	double R = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
	double Valor = (Vx*3*z*y/(pow(R,5))) - Vy*3*z*x/(pow(R,5));
	return Valor;
}

double fz(double Vz){
	return Vz;
}








	

