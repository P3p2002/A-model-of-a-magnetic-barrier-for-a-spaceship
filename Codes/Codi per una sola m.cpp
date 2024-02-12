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
	double M = 7*pow(10,6);
	
	// Permeabilitat magnètica del buit
	double mu = 4*M_PI*pow(10,-7);
	
	// Velocitat partícules carregades (m/s) //
	double Vp = 550000;
	
	// Velocitat de la nau (m/s) //
	double Vn = 11000;
	
	// Temps característic (s) //
	double t_0 = pow(10,4)/539000;
	
	int N_particules = 10;
	double dt = pow(10,-7);
	int n = 100000;
	int n_imprimir=800;
	int h;
	int j;
	int s;
	int g=0;
	double x[N_particules], y[N_particules], z[N_particules];
	double Vx[N_particules], Vy[N_particules], Vz[N_particules];
	double xocs[N_particules];
	ofstream fitxer("Posicions.txt");
	for(s=1; s<=1; s++){
		
		M = M + 10000;
		
	    // Longitud característica (m)
	   double l_0 = pow(((mu*q*M*t_0)/(4*M_PI*m)),1.0/3.0);
	   cout << "Longitud de la nau" << l/l_0 <<"Radi de la nau" <<r/l_0 <<  endl;
	   double Impactes = 0;
	   double min = 100;
	   
	   for(h=1; h<=N_particules; h++){
	
	      // Posició inicial
	      double x_0=-90 ,z_0=-r+2*h*r/N_particules;
	      double y_0=0;
	      xocs[h]=0;
		  /*for(int p=0;p<5;p++){
		  	double y_0=-r+2*h*r/5;
		  	// Normalització
	        x_0 = x_0/l_0;
	        y_0 = y_0/l_0;
	        z_0 = z_0/l_0;

	      x[h] = x_0;
          y[h] = y_0;
          z[h] = z_0;
		  }*/
       	 // Normalització
	     x_0 = x_0/l_0;
	     y_0 = y_0/l_0;
	     z_0 = z_0/l_0;

	     x[h] = x_0;
         y[h] = y_0;
         z[h] = z_0;
          // Velocitat inicial
	      double Vx_0=Vp-Vn, Vy_0=0, Vz_0=0;
	
	      // Normalització
	      Vx_0 = Vx_0*t_0/l_0;
	      Vy_0 = Vy_0*t_0/l_0;
	      Vz_0 = Vz_0*t_0/l_0;
	
	      Vx[h] = Vx_0;
	      Vy[h] = Vy_0;
	      Vz[h] = Vz_0;
}
	   
	   for(j=0; j<=n; j++){
	   for(h=1; h<=N_particules; h++){
		
		   double K1X = fx(Vx[h]);
		   double K1Y = fy(Vy[h]);
		   double K1Z = fz(Vz[h]);
		
	  	   double L1X = Fx(x[h],y[h],z[h],Vy[h],Vz[h]);
		   double L1Y = Fy(x[h],y[h],z[h],Vx[h],Vz[h]);
		   double L1Z = Fz(x[h],y[h],z[h],Vx[h],Vy[h]);
		
		   double K2X = fx(Vx[h] + dt*L1X/2);
		   double K2Y = fy(Vy[h] + dt*L1Y/2);
		   double K2Z = fz(Vz[h] + dt*L1Z/2);
		
		   double L2X = Fx(x[h] + dt*K1X/2, y[h] + dt*K1Y/2, z[h] + dt*K1Z/2, Vy[h] + dt*L1Y/2, Vz[h] + dt*L1Z/2);
		   double L2Y = Fy(x[h] + dt*K1X/2, y[h] + dt*K1Y/2, z[h] + dt*K1Z/2, Vx[h] + dt*L1X/2, Vz[h] + dt*L1Z/2);
		   double L2Z = Fz(x[h] + dt*K1X/2, y[h] + dt*K1Y/2, z[h] + dt*K1Z/2, Vx[h] + dt*L1X/2, Vy[h] + dt*L1Y/2);
		
		   double K3X = fx(Vx[h] + dt*L2X/2);
		   double K3Y = fy(Vy[h] + dt*L2Y/2);
		   double K3Z = fz(Vz[h] + dt*L2Z/2);
		
		   double L3X = Fx(x[h] + dt*K2X/2, y[h] + dt*K2Y/2, z[h] + dt*K2Z/2, Vy[h] + dt*L2Y/2, Vz[h] + dt*L2Z/2);
	       double L3Y = Fy(x[h] + dt*K2X/2, y[h] + dt*K2Y/2, z[h] + dt*K2Z/2, Vx[h] + dt*L2X/2, Vz[h] + dt*L2Z/2);
	 	   double L3Z = Fz(x[h] + dt*K2X/2, y[h] + dt*K2Y/2, z[h] + dt*K2Z/2, Vx[h] + dt*L2X/2, Vy[h] + dt*L2Y/2);
		
		   double K4X = fx(Vx[h] + dt*L3X);
		   double K4Y = fy(Vy[h] + dt*L3Y);
		   double K4Z = fz(Vz[h] + dt*L3Z);
		
		   double L4X = Fx(x[h]+ dt*K3X, y[h]+ dt*K3Y, z[h]+ dt*K3Z, Vy[h] + dt*L3Y, Vz[h] + dt*L3Z);
		   double L4Y = Fy(x[h]+ dt*K3X, y[h]+ dt*K3Y, z[h]+ dt*K3Z, Vx[h] + dt*L3X, Vz[h] + dt*L3Z);
		   double L4Z = Fz(x[h]+ dt*K3X, y[h]+ dt*K3Y, z[h]+ dt*K3Z, Vx[h] + dt*L3X, Vy[h] + dt*L3Y);
		
		   x[h] = x[h] + dt*(K1X + 2*K2X + 2*K3X + K4X)/6;
		   y[h] = y[h] + dt*(K1Y + 2*K2Y + 2*K3Y + K4Y)/6;
		   z[h] = z[h] + dt*(K1Z + 2*K2Z + 2*K3Z + K4Z)/6;
		
		   Vx[h] = Vx[h] + dt*(L1X + 2*L2X + 2*L3X + L4X)/6;
		   Vy[h] = Vy[h] + dt*(L1Y + 2*L2Y + 2*L3Y + L4Y)/6;
		   Vz[h] = Vz[h] + dt*(L1Z + 2*L2Z + 2*L3Z + L4Z)/6;
		   //cout << j%n_imprimir << "\n";
		    if(j%n_imprimir==0){
		   		fitxer<< x[h] << "\t" << y[h] << "\t"<< z[h] << "\t";
			}

           if(sqrt(pow(x[h],2)+pow(y[h],2)+pow(z[h],2)) < min){
           	min = sqrt(pow(x[h],2)+pow(y[h],2)+pow(z[h],2));
           }
           
           if(sqrt(pow(z[h],2)+pow(y[h],2)) < (r/l_0) && x[h] < l/(2*l_0) && x[h] > -l/(2*l_0) && xocs[h]==0){
			Impactes = Impactes + 1;
			xocs[h]=h;
		   }
		  }
		  if(j%n_imprimir==0){
		  	fitxer << "\n";
		   }
	
		  if(j%n_imprimir==0){
			char buffer[150];
			sprintf(buffer,"data%d",g);
			ofstream fitxer2(buffer);
			g+=1;
			for(int h=0;h<N_particules;h++){
				fitxer2 << x[h] << "\t" << y[h] << "\t"<< z[h] << "\t";
			}
		  	fitxer2.close();
		  }
	   }
	   cout<< M;
	   cout<< "________";
	   cout<< min*l_0;
	   cout<< "________";
	   cout<< Impactes << endl;
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












	

