#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdarg.h>

#include "magnetic_field_on_particle_axis_v2.h"


#define ME 511000.0     // in eV
#define CLIGHT 3E8      // vacuum speed of light in m per second
#define BSOURCE (3.6*0.7)  // Bfield at source in Tesla
#define ZSTART -39.0 // middle of WGTS
#define ZEND -15.0
#define PRINTFLAG 0
double E_source, theta_source; // values at source

// %%%%%%%%%%%%%%%%%% INTERPOLATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double interpolate_downwards(double x, double *x_array, double *y_array, 
                             int N_array){
	/* x_array is decreasing */
	int i;
	double y_interpol;
	
	y_interpol = -1.;
	if (x>=x_array[0]) {
		y_interpol = y_array[0];
	} else {
		if (x<=x_array[N_array-1]) {
			y_interpol = y_array[N_array-1];
		} else {
			i = 0;
			while ((x<x_array[i])&&(i<N_array-1)) {
				i++;
			}
			y_interpol =  (  (x-x_array[i])  *y_array[i-1] 
						   + (x_array[i-1]-x)*y_array[i]  )
			/ (x_array[i-1]-x_array[i]);
		}
	}
	return y_interpol;
}  

double interpolate_upwards(double x, double *x_array, double *y_array, 
						   int N_array){
	/* x_array is increasing */
	int i;
	double y_interpol;
	
	y_interpol = -1.;
	if (x<=x_array[0]) {
		y_interpol = y_array[0];
	} else {
		if (x>=x_array[N_array-1]) {
			y_interpol = y_array[N_array-1];
		} else {
			i = 0;
			while ((x>x_array[i])&&(i<N_array-1)) {
				i++;
			}
			y_interpol =  (  (x_array[i]-x)  *y_array[i-1] 
						   + (x-x_array[i-1])*y_array[i]  )
			/ (x_array[i]-x_array[i-1]);
		}
	}
	return y_interpol;
}  
// %%%%%%%%%%%%%%%% END OF INTERPOLATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double riemann_integral(double a, double b, 
						double const (*func)(double), double eps){
	double x,deltax;
	deltax = (b-a)*eps;
	double sum=0.0;
	x=a;
	double zwi=func(x)*deltax;
    if (zwi>=0){
		sum+=zwi;
	} else {
        if (PRINTFLAG>0) printf("in Riemann-Start: %f %f\n",deltax,func(x));
		return -1.0;
	}
	while (x<b){
	  x+=deltax;
   	  zwi=func(x)*deltax;
	  if (zwi>=0){
        if (PRINTFLAG>0) printf("ok in Riemann: %f %f %f\n",x,deltax,func(x));
		sum+=zwi;
	  } else {
        if (PRINTFLAG>0) printf("in Riemann: %f %f %f\n",x,deltax,func(x));
		return -10.0;
	  }
	}
	return sum;
}



double bfield(double z){
    return interpolate_downwards(z,s_vec,B_vec,N_VEC);
}

double p_par_square(double z){
	double zwi;
	zwi = (E_source*E_source + 2.0*E_source*ME)*(1-sin(theta_source)*sin(theta_source)*bfield(z)/BSOURCE);
	// zwi += deltaU(z)*deltaU(z) - 2*deltaU(z)*(E_source+ME);
	return zwi;
}

double const tof_integrand(double z){
	double p_par2 = p_par_square(z);
    double b = bfield(z);
	double one_over_v, deltaE_synchrotron;
	if (p_par2>0){
		one_over_v = (E_source + ME)/sqrt(p_par2)/CLIGHT;
        deltaE_synchrotron = 0.39*b*b*b/BSOURCE * sin(theta_source)*sin(theta_source)*E_source;
        if (PRINTFLAG>0) printf("%f %f %f %f %f\n",z,p_par2, E_source, deltaE_synchrotron,one_over_v*deltaE_synchrotron);
        return one_over_v*deltaE_synchrotron;
        
	} else {
        if (PRINTFLAG>0) printf("%f %f %f\n",z,p_par2, E_source);
		return -1.0;
	}
	
}	




//______________________________________________________________________________
int main(){
    double DeltaE;
//    E_source = 30474.0;
//    E_source = 17824.0;
    E_source = 18600.0;
    for (int i=0;i<204;i++){
        theta_source = i*0.25/180.0*M_PI;    // 50.77/180.0*M_PI;
        DeltaE = riemann_integral(ZSTART,ZEND,tof_integrand,1.0e-4);
        printf("%f %f\n",theta_source/M_PI*180.0,DeltaE);
    }
	return 0;
}
