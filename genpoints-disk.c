// Kuzmin disk

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#define G 1.0 // gravitational constant

// uniformly distributed random variable
double
u (void) {
	return (double)rand()/RAND_MAX;
}
// normally distributed random variable
double
n (void) {
	const double e=2.71828182845904523536;
	double U,V,x;
	do {
		U = u();
		V = u();
		x = sqrt(8/e)*(V-0.5)/U;
	} while (x*x>-4*log(U));
	return x;
}

int
main (void) {
	const int N=10000; // number of points
	const double a=10; // length scale
	const double M=N; // total mass
	const double sigma_x=1.0; // standard deviation for x
	const double rsigma_v=0.1; // relative standard deviation for v
	const double sigma_phi=5.0; // standard deviation for phi in degrees
	int i;
	double R,r,phi,v;
	double x,y,z,vx,vy,vz;

	srand(getpid());
	for (i=0; i<N; ++i) {
		R = a*sqrt(1/pow(1-u(),2)-1);
		phi = 2*M_PI*u();
		x = sigma_x*n();
		y = R*cos(phi);
		z = R*sin(phi);
		r = sqrt(R*R+a*a);
		v = sqrt(G*M/r)*R/r;
		v += rsigma_v*n()*v;
		phi += sigma_phi*n()*2*M_PI/360;
		vx = 0;
		vy = -v*sin(phi);
		vz = v*cos(phi);
		printf("%g %g %g \t %g %g %g \t %g\n",x,y,z,vx,vy,vz,M/N);
	}
	return 0;
}
