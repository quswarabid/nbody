// mass distributions with spherical symmetry

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
	const int N=5000;
	int i;
	double r,v,phi,theta;
	double x,y,z,vx,vy,vz;

	srand(getpid());

	for (i=0; i<N; ++i) {
		// mass distribution
#if 1
		// uniform sphere
		const double M=N; // total mass
		const double rmax=500; // radius of the sphere
		double rho; // volume density
		rho = 3*M/(4*M_PI*rmax*rmax*rmax);
		r = rmax*cbrt(u());
		v = r*sqrt(4*M_PI*rho*G/3);
#endif
#if 0
		// singular isothermal sphere
		const double M=N; // total mass
		const double r0=500; // radius of the sphere
		double rho0; // volume density scale
		rho0 = M/(4*M_PI*r0*r0*r0);
		r = r0*u();
		v = r0*sqrt(4*M_PI*G*rho0);
#endif
#if 0
		// Hernquist density profile
		const double M=N; // total mass
		const double a=10; // length scale
		double z; // random variable
		z = sqrt(u());
		r = a*z/(1-z);
		v = sqrt(G*M*r)/(r+a);
#endif

		phi = 2*M_PI*u();
		theta = acos(1-2*u());
		x = r*sin(theta)*cos(phi);
		y = r*sin(theta)*sin(phi);
		z = r*cos(theta);

		// velocity distribution
#if 1
		// random tangential velocity
		double alpha; // direction angle
		double dphi[3]; // basis vector of a tangential plane
		double dtheta[3]; // basis vector of a tangential plane
		double a,b; // coefficients for the basis vectors
		dphi[0] = -sin(phi);
		dphi[1] = cos(phi);
		dphi[2] = 0;
		dtheta[0] = cos(theta)*cos(phi);
		dtheta[1] = cos(theta)*sin(phi);
		dtheta[2] = -sin(theta);
#if 1
		alpha = 2*M_PI*u(); // random angle
#else
		alpha = 0; // in parallel with the xy plane
#endif
		a = v*cos(alpha);
		b = v*sin(alpha);
		vx = a*dphi[0]+b*dtheta[0];	
		vy = a*dphi[1]+b*dtheta[1];	
		vz = a*dphi[2]+b*dtheta[2];	
#endif
#if 0
		// weird case: random normal directions
		phi = 2*M_PI*u();
#if 1
		theta = acos(1-2*u()); // random
#else
		theta = M_PI/2; // in parallel with the xy plane
#endif
		vx = v*sin(theta)*cos(phi);
		vy = v*sin(theta)*sin(phi);
		vz = v*cos(theta);
#endif
#if 0
		// quite weird case: normal directions
		vx = v*sin(theta)*cos(phi);
		vy = v*sin(theta)*sin(phi);
		vz = v*cos(theta);
#endif
		printf("%g %g %g \t %g %g %g \t %g\n",x,y,z,vx,vy,vz,M/N);
	}
	return 0;
}
