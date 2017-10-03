/*
 * Straightforward and naive N-body code.
 * Initial version.
 * 	- Euler's integrator with adaptive step;
 *	- Singularities are prevented by dropping forces to zero
 *	  in the vicinity of points;
 *	- Show speed by color.
 * Aug-Sep 2017
 * Alexander Mukhin
 */

#define M 1000 // phase cube |x,y,z|<=M
#define Px 1000 // projection plane x=Px
#define Ex 4000 // eye point (Ex,0,0)
#define WH 700 // screen width=height
#define	G 1 // gravitational constant
#define FPS 100 // frames per second
#define EPS 0.001 // maximum absolute increment of position and speed
#define DRMIN 10 // distance below which the force is set to zero
#define ORDS 2.5 // orders of magnitude for color output
#define DMAX 5 // maximum displacement of points between consecutive frames

#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <X11/Xlib.h>
#include <X11/keysym.h>
#include <time.h>
#include <float.h>

// Point
struct p {
	double x,y,z; // position
	double px,py,pz; // previous position
	double vx,vy,vz; // velocity
	double ax,ay,az; // acceleration
	double m; // mass
	struct p *next; // link
};

// Length of a vector
double
len (double x, double y, double z) {
	return sqrt(x*x+y*y+z*z);
}

double T,U,P_x,P_y,P_z;	// kinetic energy, potential energy, and momentum
int C; // number of collisions (approaches below DRMIN)

// Central projection
void
proj (double x, double y, double z, int *u, int *v) {
	const double A=M*(Ex-Px)/(double)(Ex-M); // projection screen size: |yp,zp|<=A
	double t; // parameter
	double yp,zp; // coordinates of the projection point

	// find projection coordinates in the phase space
	t = (Ex-Px)/(Ex-x);
	yp = y*t;
	zp = z*t;
	// convert them to the pixel space
	*u = (int)floor((A+yp)*WH/(2*A));
	*v = (int)floor((A-zp)*WH/(2*A));
}

// Trace evolution from one frame to another
double
step (struct p *ppl) {
	struct p *p,*q;	// pointers to population
	double dx,dy,dz,dr; // distance vector and its length
	double v,vmax,a,amax; // current and maximum velocity and acceleration
	double t,dt; // time elapsed and time step
	double d,dmax; // current and maximum displacement

	// Save current positions
	for (p=ppl; p; p=p->next) {
		p->px = p->x;
		p->py = p->y;
		p->pz = p->z;
	}
	// Keep advancing until we accumulate enough displacement
	t = 0;
	do {
		// One step of Euler's method
		T = U = P_x = P_y = P_z = 0;
		vmax = amax = 0;
		for (p=ppl; p; p=p->next)
			p->ax = p->ay = p->az = 0;
		// Find accelerations
		for (p=ppl; p; p=p->next) {
			// Loop for the pairs
			for (q=p->next; q; q=q->next) {
				dx = q->x-p->x;
				dy = q->y-p->y;
				dz = q->z-p->z;
				dr = len(dx, dy, dz);
				if (dr < DRMIN) {
					// Singularity
					++C;
					U += -G*p->m*q->m/DRMIN;
					continue;
				}
				// Acceleration of p caused by q
				p->ax += (dx/dr)*G*q->m/(dr*dr);
				p->ay += (dy/dr)*G*q->m/(dr*dr);
				p->az += (dz/dr)*G*q->m/(dr*dr);
				// Acceleration of q caused by p
				q->ax += (-dx/dr)*G*p->m/(dr*dr);
				q->ay += (-dy/dr)*G*p->m/(dr*dr);
				q->az += (-dz/dr)*G*p->m/(dr*dr);
				// Update potential energy
				U += -G*p->m*q->m/dr;
			}
			// Collect statistics
			v = len(p->vx, p->vy, p->vz);
			if (v > vmax)
				vmax = v;
			a = len(p->ax, p->ay, p->az);
			if (a > amax)
				amax = a;
			// Update kinetic energy and momentum
			T += p->m*v*v/2;
			P_x += p->m*p->vx;
			P_y += p->m*p->vy;
			P_z += p->m*p->vz;
		}
		// Find appropriate step
		// Assume amax or vmax is not zero
		dt = DBL_MAX;
		if (amax != 0)
			dt = EPS/amax;
		if (vmax!=0 && EPS/vmax<dt)
			dt = EPS/vmax;
		// Modify positions and velocities
		for (p=ppl; p; p=p->next) {
			p->x += p->vx*dt;
			p->y += p->vy*dt;
			p->z += p->vz*dt;
			p->vx += p->ax*dt;
			p->vy += p->ay*dt;
			p->vz += p->az*dt;
		}
		// Count time elapsed
		t += dt;
		// Find maximum displacement
		dmax = 0;
		for (p=ppl; p; p=p->next) {
			dx = p->x-p->px;
			dy = p->y-p->py;
			dz = p->z-p->pz;
			d = len(dx, dy, dz);
			if (d > dmax)
				dmax = d;
		}
	} while (dmax < DMAX);
	// Return time elapsed
	return t;
}

// Color space conversion
int
hsv2rgb (double h, double s, double v) {
	double c,x,m;
	double r=0,g=0,b=0;
	// Using definitions from Wiki
	h *= 6;
	c = s*v;
	x = c*(1-fabs(h-2*floor(h/2)-1));
	if (0<=h && h<1)
		r=c, g=x, b=0;
	else if (1<=h && h<2)
		r=x, g=c, b=0;
	else if (2<=h && h<3)
		r=0, g=c, b=x;
	else if (3<=h && h<4)
		r=0, g=x, b=c;
	else if (4<=h && h<5)
		r=x, g=0, b=c;
	else if (5<=h && h<=6)
		r=c, g=0, b=x;
	m = v-c;
	r += m, g += m, b += m;
        return ((int)floor(r*255)<<16)+((int)floor(g*255)<<8)+(int)floor(b*255);
}

// Magnitude to color
int
mag2clr (double v) {
	double m;
	m = log10(1+v);
	if (m > ORDS)
		m = ORDS;
	return hsv2rgb(0.66*m/ORDS, 1.0, 1.0);
}

int
main (void) {
	struct p *ppl; // population
	struct p *p,*q;	// pointers to population
	double t; // time elapsed
	int u,v,pu,pv; // coordinates in pixel space
	double s; // speed
	struct timespec ts; // for delay between frames
	// X stuff
	Display *disp;
	Window win;
	XEvent evt;
        int scr;
        GC gc;
	Pixmap pm;
	KeySym ks;

	// Initialize the population
	ppl = q = NULL;
	while (!feof(stdin)) {
		p = (struct p *)malloc(sizeof(struct p));
		scanf("%lg %lg %lg %lg %lg %lg %lg\n",
		       &p->x, &p->y, &p->z, &p->vx, &p->vy, &p->vz, &p->m);
		p->next = NULL;
		if (q)
			q->next = p;
		else
			ppl = p; // first point
		q = p;
	}

	// X init
	disp = XOpenDisplay(NULL);
	if (disp == NULL) {
		fprintf(stderr, "Cannot open display\n");
		return 1;
	}
	scr = DefaultScreen(disp);
	win = XCreateSimpleWindow(disp, RootWindow(disp,scr),
				  10, 10, WH, WH,
				  1, WhitePixel(disp,scr),
				  BlackPixel(disp,scr));
	XSelectInput(disp, win, ExposureMask|KeyPressMask);
	XMapWindow(disp, win);
	gc = DefaultGC(disp, scr);
	pm = XCreatePixmap(disp, win, WH, WH, DefaultDepth(disp,scr));
	XFillRectangle(disp, pm, gc, 0, 0, WH, WH);
	XSetForeground(disp, gc, WhitePixel(disp,scr));

	// Event loop
	t = 0;
	ts.tv_sec = 0;
	ts.tv_nsec = 1000000000/FPS;
	for (;;) {
		XNextEvent(disp, &evt);
		if (evt.type == Expose) {
			// Drop all pending Exposure events
			while (XCheckMaskEvent(disp, ExposureMask, &evt));
			// Advance the system from one frame to another
			t += step(ppl);
			// Draw lines between consecutive positions
			for (p=ppl; p; p=p->next) {
				proj(p->x, p->y, p->z, &u, &v);
				proj(p->px, p->py, p->pz, &pu, &pv);
				s = len(p->vx, p->vy, p->vz);
				XSetForeground(disp, gc, mag2clr(s));
				XDrawLine(disp, pm, gc, pu, pv, u, v);
				XSetForeground(disp, gc, 0xFFFFFF);
				XDrawPoint(disp, pm, gc, u, v);
			}
			// Show pixmap
			XCopyArea(disp, pm, win, gc, 0, 0, WH, WH, 0, 0);
			// Wait
//			usleep(1000000/FPS);
			nanosleep(&ts, NULL);
			// Send itself an exposure event
			// to display the next frame
			evt.type = Expose;
			XSendEvent(disp, win, False, 0, &evt);
		}
		if (evt.type == KeyPress) {
			ks = XLookupKeysym(&evt.xkey, 0);
			if (ks == XK_q)
				// Quit
				break;
			if (ks == XK_c) {
				// Check conservation laws
				printf("t=%g, C=%d\n", t, C);
				printf("T+U=%g, ", T+U);
				printf("P=%g\n", len(P_x,P_y,P_z));
				printf("\n");
			}
			if (ks == XK_d) {
				// Dump state
				printf("t=%g, C=%d\n", t, C);
				for (p=ppl; p; p=p->next) {
					printf("%g,%g,%g ",p->x,p->y,p->z);
					printf("%g,%g,%g ",p->vx,p->vy,p->vz);
					printf("%g,%g,%g\n",p->ax,p->ay,p->az);
				}
				printf("\n");
			}
		}
	}

	// Clean up
	p = ppl;
	do {
		q = p->next;
		free(p);
		p = q;
	} while (q);
	XFreePixmap(disp, pm);
	XCloseDisplay(disp);

	return 0;
}
