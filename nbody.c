/*
	Simple collisionless N-body code:
	- fixed-step leapfrog integrator;
	- Plummer softening;
	- Barnes-Hut tree.

	Copyright (c) 2018 Alexander Mukhin
	MIT License
 */

#define DT 0.1 // time step
#define M 1500 // phase cube |x,y,z|<=M
#define Px 2000 // projection plane distance from the centre
#define Ox 5000 // projection origin distance from the centre
#define WH 700 // screen width=height
#define	G 1 // gravitational constant
#define FPS 25 // frames per second
#define ORDS 2.5 // orders of magnitude for color output
#define DMAX 5 // maximum displacement of points between consecutive frames
#define EPS 5.0 // gravity softening radius
#define MAXLVL 20 // maximal height of the octree
#define NAVG 10000 // number of steps to average T and V
// thresholds below which the CM approximation is used
#define ATHR 0.5 // opening angle
#define MTHR 0.1 // fraction of the total mass

#define _POSIX_C_SOURCE 200809L
#define M_PI 3.14159265358979323846

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <X11/Xlib.h>
#include <X11/keysym.h>
#include <time.h>
#include <sys/time.h>

// main data structures
// point
struct p {
	double x,y,z; // position
	double px,py,pz; // previous position
	double vx,vy,vz; // velocity
	double ax,ay,az; // acceleration
	double m; // mass
	struct p *next; // link to the next point in a leaf node
	int ins; // ins() result (1=success,0=fail)
};
// octree node
struct o {
	int lvl; // level
	double xmin,xmax,ymin,ymax,zmin,zmax; // bounds
	struct p *ps; // list of points
	struct o *os[8]; // octants
	double m; // total mass
	double x_cm,y_cm,z_cm; // centre of mass
	double rmax; // maximal distance between octree points and their CM
};

// basic vector algebra
struct vec {
	double x,y,z;
};
// vector addition
struct vec
add (struct vec a, struct vec b) {
	struct vec c;
	c.x = a.x+b.x;
	c.y = a.y+b.y;
	c.z = a.z+b.z;
	return c;
}
// vector subtraction
struct vec
sub (struct vec a, struct vec b) {
	struct vec c;
	c.x = a.x-b.x;
	c.y = a.y-b.y;
	c.z = a.z-b.z;
	return c;
}
// multiplication by a scalar
struct vec
mul (struct vec a, double b) {
	struct vec c;
	c.x = a.x*b;
	c.y = a.y*b;
	c.z = a.z*b;
	return c;
}
// inner product
double
prod (struct vec a, struct vec b) {
	return a.x*b.x+a.y*b.y+a.z*b.z;
}
// length of a vector
double
len (double x, double y, double z) {
	return sqrt(x*x+y*y+z*z);
}

// global variables
// projection parameters
int Phi=0,Theta=90; // angles in degrees
struct vec O={Ox,0,0}; // projection origin
struct vec P={Px,0,0}; // projection plane centre
struct vec Dphi={0,1,0}; // basis vector of the projection plane
struct vec Dtheta={0,0,1}; // basis vector of the projection plane
struct vec Dr={1,0,0}; // normal vector of the projection plane
// physical parameters
double T,U; // kinetic and potential energy
double P_x,P_y,P_z; // linear momentum
double L_x,L_y,L_z; // angular momentum
double V; // Clausius' virial
double Ts[NAVG],Vs[NAVG]; // samples of T and V
double Tsum,Vsum; // sum of T and V over NAVG steps
// statistics
struct {
	int steps; // number of steps
	int reins,prune,lift,purge; // event counters
	int cur_rq,max_rq; // current and max re-insert queue length
	int acc_cm,acc_direct; // tree algorithm performance
	int cur_t,max_t; // accel() timings
	int cur_leaves,max_leaves; // current and max number of leaves
	int cur_ints,max_ints; // current and max number of internal nodes
	int cur_h,max_h; // current and max tree height
	int cur_maxps,max_maxps; // current and maximal ps list length
	int folen,req_malloc,req_fo; // memory statistics
} Stat;
// list of free octree nodes
struct o *FO;

// octree-related functions
// check if point *p is in the node *o
int
in (struct o *o, struct p *p) {
	return (o->xmin<=p->x && p->x<o->xmax) && \
	       (o->ymin<=p->y && p->y<o->ymax) && \
	       (o->zmin<=p->z && p->z<o->zmax);
}
// insert a point into an octree
void
ins (struct o *o, struct p *p) {
	double xmin[8],xmax[8]; // octant bounds
	double ymin[8],ymax[8]; // octant bounds
	double zmin[8],zmax[8]; // octant bounds
	int no,no2; // octant numbers
	struct p *q,*qn; // pointers
	// initialize octant boundaries
	// octant #0
	xmin[0] = (o->xmin+o->xmax)/2;
	xmax[0] = o->xmax;
	ymin[0] = (o->ymin+o->ymax)/2;
	ymax[0] = o->ymax;
	zmin[0] = (o->zmin+o->zmax)/2;
	zmax[0] = o->zmax;
	// octant #1
	xmin[1] = o->xmin;
	xmax[1] = xmin[0];
	ymin[1] = ymin[0];
	ymax[1] = ymax[0];
	zmin[1] = zmin[0];
	zmax[1] = zmax[0];
	// octant #2
	xmin[2] = xmin[1];
	xmax[2] = xmax[1];
	ymin[2] = o->ymin;
	ymax[2] = ymin[1];
	zmin[2] = zmin[1];
	zmax[2] = zmax[1];
	// octant #3
	xmin[3] = xmax[2];
	xmax[3] = o->xmax;
	ymin[3] = ymin[2];
	ymax[3] = ymax[2];
	zmin[3] = zmin[2];
	zmax[3] = zmax[2];
	// octant #4
	xmin[4] = xmin[0];
	xmax[4] = xmax[0];
	ymin[4] = ymin[0];
	ymax[4] = ymax[0];
	zmin[4] = o->zmin;
	zmax[4] = zmin[0];
	// octant #5
	xmin[5] = xmin[1];
	xmax[5] = xmax[1];
	ymin[5] = ymin[1];
	ymax[5] = ymax[1];
	zmin[5] = zmin[4];
	zmax[5] = zmax[4];
	// octant #6
	xmin[6] = xmin[2];
	xmax[6] = xmax[2];
	ymin[6] = ymin[2];
	ymax[6] = ymax[2];
	zmin[6] = zmin[5];
	zmax[6] = zmax[5];
	// octant #7
	xmin[7] = xmin[3];
	xmax[7] = xmax[3];
	ymin[7] = ymin[3];
	ymax[7] = ymax[3];
	zmin[7] = zmin[6];
	zmax[7] = zmax[6];
	// check bounds
	if (!in(o,p)) {
		p->ins = 0; // ins() failed
		return;
	}
	if (o->lvl==MAXLVL) {
		// append to the list in this node
		p->next = o->ps;
		o->ps = p;
		p->ins = 1; // ins() succeeded
	} else {
		if (o->ps) {
			// leaf node
			// it becomes internal
			q = o->ps;
			o->ps = NULL;
			while (q) {
				qn = q->next;
				ins(o,q);
				q = qn;
			}
			ins(o,p);
		} else {
			// internal node
			// find octant
			for (no=0; no<8; ++no)
				if ((xmin[no]<=p->x && p->x<xmax[no]) && \
				    (ymin[no]<=p->y && p->y<ymax[no]) && \
				    (zmin[no]<=p->z && p->z<zmax[no]))
					break;
			if (o->os[no]) {
				// call recursively for the octant found
				ins(o->os[no],p);
			} else {
				// allocate leaf node
				// first, try to get it from the FO list
				if (FO) {
					o->os[no] = FO;
					FO = FO->os[0];
					--Stat.folen;
					++Stat.req_fo;
				} else {
					// as a last resort, call malloc
					o->os[no] = malloc(sizeof(struct o));
					++Stat.req_malloc;
				}
				// initialize octree node
				o->os[no]->lvl = o->lvl+1;
				o->os[no]->xmin = xmin[no];
				o->os[no]->xmax = xmax[no];
				o->os[no]->ymin = ymin[no];
				o->os[no]->ymax = ymax[no];
				o->os[no]->zmin = zmin[no];
				o->os[no]->zmax = zmax[no];
				for (no2=0; no2<8; ++no2)
					o->os[no]->os[no2] = NULL;
				o->os[no]->ps = p;
				p->next = NULL;
				p->ins = 1; // ins() succeeded
			}
		}
	}
}
// drift points
struct p *
drift (struct o *o, double dt, struct p *r) {
	struct p *p,*pprev,*q; // point pointers
	int no; // octant number
	if (o->ps) {
		// leaf node
		// drift points
		// delete and place them on the re-insert queue
		// if they left their octant
		p = o->ps;
		pprev = NULL;
		while (p) {
			// drift the point
			p->x += p->vx*dt;
			p->y += p->vy*dt;
			p->z += p->vz*dt;
			// see if it moved out of its octant
			if (!in(o,p)) {
				// point p is not in its place
				q = p; // save pointer to it
				// remove the point from the list
				if (pprev)
					// mid-list
					p = pprev->next = p->next;
				else
					// head
					o->ps = p = p->next;
				// add it to the re-insert queue
				q->next = r;
				r = q;
				// update statistics
				++Stat.reins;
			} else {
				// point stayed in its octant
				// proceed to the next point
				pprev = p;
				p = p->next;
			}
		}
	} else {
		// internal node
		for (no=0; no<8; ++no)
			// call for non-empty suboctants
			// and concatenate re-insert queues from them
			if (o->os[no])
				r = drift(o->os[no],dt,r);
	}
	return r;
}
// prune the octree: remove empty leaf nodes
int
prune (struct o *o) {
	int no; // octant number
	if (o->ps)
		// skip leaf nodes
		return 0;
	// call for non-empty subnodes
	for (no=0; no<8; ++no)
		if (o->os[no])
			if (prune(o->os[no]))
				o->os[no] = NULL;
	// should we prune itself?
	if (o->lvl==0)
		// never prune the root node
		return 0;
	for (no=0; no<8; ++no)
		if (o->os[no])
			break;
	if (no==8) {
		// prune if no subnodes left
		// add o to the free list
		o->os[0] = FO; // link via os[0]
		FO = o;
		++Stat.folen;
		// update statistics
		++Stat.prune;
		return 1;
	} else
		return 0;
}
// lift up leaf nodes
void
lift (struct o *o) {
	int no,no2; // octant number
	int c; // counter
	if (o->ps)
		// skip leaf nodes
		return;
	// call for non-empty subnodes
	for (no=0; no<8; ++no)
		if (o->os[no])
			lift(o->os[no]);
	// should we turn itself into a leaf?
	c = 0;
	for (no=0; no<8; ++no)
		if (o->os[no]) {
			++c;
			no2 = no;
		}
	if (c==1) {
		// only one descendant left
		// if it's a leaf, lift it up
		if (o->os[no2]->ps) {
			// grab its points, so we become a leaf
			o->ps = o->os[no2]->ps;
			// free subnode
			// add o->os[no2] to the free list
			o->os[no2]->os[0] = FO; // link via os[0]
			FO = o->os[no2];
			++Stat.folen;
			// zeroize subnode pointer
			o->os[no2] = NULL;
			// update statistics
			++Stat.lift;
		}
	}
}

// purge runaway points
int
purge (struct p **ppl, int np) {
	struct p *p; // point pointer
	int n,m; // point counters
	m = 0;
	for (n=0; n<np; ++n) {
		p = ppl[n];
		if (p->ins) {
			// point is within bounds
			ppl[m] = p;
			++m;
		} else {
			// point is out of bounds, skip it
			// update purge statistics
			++Stat.purge;
		}
	}
	return m;
}

// compute centres of mass
void
mkcms (struct o *o) {
	struct p *p; // pointer
	int no; // octant number
	double m,smx,smy,smz; // accumulators
	double r,rmax; // current and maximal distances between points and CM
	double dx,dy,dz; // distance vector components
	int c; // counter
	m = smx = smy = smz = 0;
	if (o->ps) {
		// leaf node
		// compute total mass and the centre of mass
		c = 0;
		for (p=o->ps; p; p=p->next) {
			m += p->m;
			smx += p->m*p->x;
			smy += p->m*p->y;
			smz += p->m*p->z;
			++c;
		}
		o->m = m;
		o->x_cm = smx/m;
		o->y_cm = smy/m;
		o->z_cm = smz/m;
		// compute maximal distance between points and their CM
		rmax = 0;
		for (p=o->ps->next; p; p=p->next) {
			dx = p->x-o->x_cm;
			dy = p->y-o->y_cm;
			dz = p->z-o->z_cm;
			r = len(dx,dy,dz);
			if (r>rmax)
				rmax = r;
		}
		o->rmax = rmax;
		// update statistics
		++Stat.cur_leaves;
		if (c>Stat.cur_maxps)
			Stat.cur_maxps = c;
		if (o->lvl>Stat.cur_h)
			Stat.cur_h = o->lvl;
	} else {
		// internal node
		// call for subnodes
		for (no=0; no<8; ++no)
			if (o->os[no])
				mkcms(o->os[no]);
		// sum results of subnodes
		for (no=0; no<8; ++no)
			if (o->os[no]) {
				m += o->os[no]->m;
				smx += o->os[no]->m*o->os[no]->x_cm;
				smy += o->os[no]->m*o->os[no]->y_cm;
				smz += o->os[no]->m*o->os[no]->z_cm;
			}
		// compute total mass and the centre of mass
		o->m = m;
		o->x_cm = smx/m;
		o->y_cm = smy/m;
		o->z_cm = smz/m;
		// compute maximal distance between points and their CM
		rmax = 0;
		for (no=0; no<8; ++no)
			if (o->os[no]) {
				dx = o->os[no]->x_cm-o->x_cm;
				dy = o->os[no]->y_cm-o->y_cm;
				dz = o->os[no]->z_cm-o->z_cm;
				r = len(dx,dy,dz)+o->os[no]->rmax;
				if (r>rmax)
					rmax = r;
			}
		o->rmax = rmax;
		// update statistics
		++Stat.cur_ints;
	}
}

// compute accelerations infulenced on p by the ensemble o
// and return the potential energy of the point p
// in the gravitational field of the ensemble o
double
acc (struct p *p, struct o *o, double totalmass) {
	double dx,dy,dz,dr,dr1; // distance vectors
	struct p *q; // pointer
	int no; // octant number
	double U; // potential energy
	U = 0;
	if (o->ps) {
		// leaf node
		// sum for all points in the list
		for (q=o->ps; q; q=q->next) {
			if (q==p)
				continue;
			dx = q->x-p->x;
			dy = q->y-p->y;
			dz = q->z-p->z;
			dr = len(dx,dy,dz);
			dr1 = sqrt(dr*dr+EPS*EPS); // Plummer softening
			// acceleration of p caused by q
			p->ax += (dx/dr1)*G*q->m/(dr1*dr1);
			p->ay += (dy/dr1)*G*q->m/(dr1*dr1);
			p->az += (dz/dr1)*G*q->m/(dr1*dr1);
			// accumulate potential energy
			U += -G*p->m*q->m/dr1;
		}
		// update statistics
		++Stat.acc_direct;
	} else {
		// internal node
		// is CM approximation applicable?
		dx = o->x_cm-p->x;
		dy = o->y_cm-p->y;
		dz = o->z_cm-p->z;
		dr = len(dx,dy,dz);
		if (o->rmax<ATHR*dr && o->m<MTHR*totalmass) {
			// apply CM approximation
			dr1 = sqrt(dr*dr+EPS*EPS); // apply Plummer softening
			// acceleration of p caused by points of o
			p->ax += (dx/dr1)*G*o->m/(dr1*dr1);
			p->ay += (dy/dr1)*G*o->m/(dr1*dr1);
			p->az += (dz/dr1)*G*o->m/(dr1*dr1);
			// compute potential energy
			U = -G*p->m*o->m/dr1;
			// update statistics
			++Stat.acc_cm;
		} else {
			// call for all non-empty subnodes
			for (no=0; no<8; ++no)
				if (o->os[no])
					U += acc(p,o->os[no],totalmass);
		}
	}
	return U;
}

// update velocities for every point of the population
// and return the sum of potential energies of interaction
// between each pair of points
double
kick (struct p **ppl, int np, struct o *o, double dt) {
	struct p *p; // point pointer
	int n; // point counter
	struct timeval t1,t2; // time values
	int d; // time difference in us
	double U; // potential energy
	U = 0;
	// reset acc statistics
	Stat.acc_cm = Stat.acc_direct = 0;
	// record timestamp
	gettimeofday(&t1,NULL);
	// loop for points
#pragma omp parallel for private(p)
	for (n=0; n<np; ++n) {
		p = ppl[n];
		p->ax = p->ay = p->az = 0;
		U += acc(p,o,o->m);
		p->vx += p->ax*dt;
		p->vy += p->ay*dt;
		p->vz += p->az*dt;
	}
	// record timestamp
	gettimeofday(&t2,NULL);
	// update statistics
	d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
	Stat.cur_t = d;
	if (d>Stat.max_t)
		Stat.max_t = d;
	return U/2;
}

// trace evolution from one frame to another
int
step (struct p **ppl, int np, struct o *o, double dt) {
	struct p *p,*q; // point pointers
	int n; // point counter
	struct p *r; // re-insert queue
	double v; // speed
	double dx,dy,dz; // displacement components
	double d,dmax; // current and maximum displacement
	// save the current positions
	for (n=0; n<np; ++n) {
		p = ppl[n];
		p->px = p->x;
		p->py = p->y;
		p->pz = p->z;
	}
	// keep advancing until we accumulate enough displacement
	do {
		// update step count
		++Stat.steps;
		// find kinetic energy, linear and angular momenta,
		// and virial
		T = P_x = P_y = P_z = L_x = L_y = L_z = V = 0;
		for (n=0; n<np; ++n) {
			p = ppl[n];
			v = len(p->vx,p->vy,p->vz);
			T += p->m*v*v/2;
			P_x += p->m*p->vx;
			P_y += p->m*p->vy;
			P_z += p->m*p->vz;
			L_x += p->m*(p->y*p->vz-p->vy*p->z);
			L_y += p->m*(-p->x*p->vz+p->vx*p->z);
			L_z += p->m*(p->x*p->vy-p->vx*p->y);
			V += p->m*(p->x*p->ax+p->y*p->ay+p->z*p->az);
		}
		// compute the sum of T and V over the last NAVG steps
		Tsum += T;
		Tsum -= Ts[Stat.steps%NAVG];
		Ts[Stat.steps%NAVG] = T;
		Vsum += V;
		Vsum -= Vs[Stat.steps%NAVG];
		Vs[Stat.steps%NAVG] = V;
		// one step of the leapfrog method
		// drift
		r = drift(o,dt,NULL);
		// keep the tree tidy
		prune(o);
		lift(o);
		// insert points placed on the re-insert queue
		Stat.cur_rq = 0;
		p = r;
		while (p) {
			q = p->next; // save the next pointer
			ins(o,p); // NB: ins() changes p->next
			p = q; // proceed to the next point
			// update statistics
			++Stat.cur_rq;
		}
		// track max_rq
		if (Stat.cur_rq > Stat.max_rq)
			Stat.max_rq = Stat.cur_rq;
		// purge runaway points
		np = purge(ppl,np);
		// exit if none or one point left
		if (np<=1)
			break;
		// reset tree statistics
		Stat.cur_leaves = Stat.cur_ints = 0;
		Stat.cur_h = Stat.cur_maxps = 0;
		// recompute the centres of mass
		// as a side effect, this updates tree statistics
		mkcms(o);
		// track maximum values
		if (Stat.cur_leaves>Stat.max_leaves)
			Stat.max_leaves = Stat.cur_leaves;
		if (Stat.cur_ints>Stat.max_ints)
			Stat.max_ints = Stat.cur_ints;
		if (Stat.cur_h>Stat.max_h)
			Stat.max_h = Stat.cur_h;
		if (Stat.cur_maxps>Stat.max_maxps)
			Stat.max_maxps = Stat.cur_maxps;
		// kick
		// it computes the potential energy as well
		U = kick(ppl,np,o,dt);
		// find the maximum displacement
		dmax = 0;
		for (n=0; n<np; ++n) {
			p = ppl[n];
			dx = p->x-p->px;
			dy = p->y-p->py;
			dz = p->z-p->pz;
			d = len(dx,dy,dz);
			if (d>dmax)
				dmax = d;
		}
	} while (dmax<DMAX);
	// return the (possibly updated) number of points
	return np;
}

// central projection
void
proj (double x, double y, double z, short *u, short *v) {
	struct vec A; // point being projected
	double t; // parameter
	struct vec p; // intersection with the projection plane
	struct vec r; // vector in the projection plane
	double up,vp; // coordinates of the projection point
	A.x=x,A.y=y,A.z=z;
	t = prod(sub(P,O),Dr)/prod(sub(A,O),Dr);
	p = add(O,mul(sub(A,O),t));
	r = sub(p,P);
	up = prod(r,Dphi);
	vp = prod(r,Dtheta);
	// convert to the pixel space
	*u = (short)floor(WH*(up+M)/(2*M));
	*v = (short)floor(WH*(M-vp)/(2*M));
}

// recompute global projection parameters after angle change
void
reproj (void) {
	double phi,theta;
	phi=((double)Phi/360)*2*M_PI;
	theta=((double)Theta/360)*2*M_PI;
	P.x = Px*sin(theta)*cos(phi);
	P.y = Px*sin(theta)*sin(phi);
	P.z = Px*cos(theta);
	O.x = Ox*sin(theta)*cos(phi);
	O.y = Ox*sin(theta)*sin(phi);
	O.z = Ox*cos(theta);
	Dphi.x = -sin(phi);
	Dphi.y = cos(phi);
	Dphi.z = 0;
	Dtheta.x = -cos(theta)*cos(phi);
	Dtheta.y = -cos(theta)*sin(phi);
	Dtheta.z = sin(theta);
	Dr.x = sin(theta)*cos(phi);
	Dr.y = sin(theta)*sin(phi);
	Dr.z = cos(theta);
}

// color space conversion
int
hsv2rgb (double h, double s, double v) {
	double c,x,m;
	double r=0,g=0,b=0;
	// using definitions from Wiki
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
        return ((int)floor(r*255)<<16) + \
	       ((int)floor(g*255)<<8) + \
	        (int)floor(b*255);
}

// magnitude to color
int
mag2clr (double v) {
	double m;
	m = log10(1+v);
	if (m>ORDS)
		m = ORDS;
	return hsv2rgb(0.66*m/ORDS,1.0,1.0);
}

// draw octree edges
void
dredges (Display *disp, Pixmap pm, GC gc, struct o *o, int mode) {
	XSegment e[12]; // edges
	int no; // octant number
	int i; // edge index

	// draw internal nodes if mode==0
	// draw leaf nodes if mode==1
	if ((o->ps==NULL && mode==0) || (o->ps!=NULL && mode==1)) {
		// find projections of vertices
		proj(o->xmax,o->ymin,o->zmax,&e[0].x1,&e[0].y1);
		proj(o->xmax,o->ymax,o->zmax,&e[1].x1,&e[1].y1);
		proj(o->xmin,o->ymax,o->zmax,&e[2].x1,&e[2].y1);
		proj(o->xmin,o->ymin,o->zmax,&e[3].x1,&e[3].y1);
		proj(o->xmax,o->ymin,o->zmin,&e[4].x1,&e[4].y1);
		proj(o->xmax,o->ymax,o->zmin,&e[5].x1,&e[5].y1);
		proj(o->xmin,o->ymax,o->zmin,&e[6].x1,&e[6].y1);
		proj(o->xmin,o->ymin,o->zmin,&e[7].x1,&e[7].y1);
		// connect segments
		for (i=0; i<4; ++i) {
			e[i].x2 = e[(i+1)%4].x1;
			e[i].y2 = e[(i+1)%4].y1;
			e[i+4].x2 = e[(i+1)%4+4].x1;
			e[i+4].y2 = e[(i+1)%4+4].y1;
			e[i+8].x1 = e[i].x1;
			e[i+8].y1 = e[i].y1;
			e[i+8].x2 = e[i+4].x1;
			e[i+8].y2 = e[i+4].y1;
		}
		// draw
		XDrawSegments(disp,pm,gc,e,12);
	}
	// call for subnodes
	for (no=0; no<8; ++no)
		if (o->os[no])
			dredges(disp,pm,gc,o->os[no],mode);
}

int
main (int argc, char **argv) {
	struct p **ppl; // population
	struct p *p0,*p,*pprev; // pointers to points
	int np,n; // total number of points and point counter
	struct o *o; // octree
	short u,v,pu,pv; // coordinates in pixel space
	double s; // speed
	struct timespec ts; // for delay between frames
	enum {O_TRACES=0x1,O_EDGES} flags=0; // flags

	// X stuff
	Display *disp;
	Window win;
	XEvent evt;
        int scr;
        GC gc;
	Pixmap pm;
	KeySym ks;

	// initialize the population into a temporary linked list
	p0 = pprev = NULL;
	np = 0;
	while (!feof(stdin)) {
		p = (struct p *)malloc(sizeof(struct p));
		scanf("%lg %lg %lg %lg %lg %lg %lg\n",
		       &p->x,&p->y,&p->z,&p->vx,&p->vy,&p->vz,&p->m);
		p->next = NULL;
		if (pprev)
			pprev->next = p;
		else
			p0 = p; // first point
		pprev = p;
		++np;
	}

	// create linear array of pointers to the points
	ppl = (struct p **)malloc(np*sizeof(struct p *));
	// fill it
	p = p0;
	n = 0;
	while (p) {
		ppl[n] = p;
		p = p->next;
		ppl[n]->next = NULL;
		++n;
	}
	
	// allocate the root node of the octree
	o = malloc(sizeof(struct o));
	o->lvl = 0;
	o->ps = NULL;
	o->xmin = -M;
	o->xmax = M;
	o->ymin = -M;
	o->ymax = M;
	o->zmin = -M;
	o->zmax = M;
	// insert points to the octree
	for (n=0; n<np; ++n)
		ins(o,ppl[n]);
	// purge points not fitting the phase space
	np = purge(ppl,np);
	// compute the centres of mass
	mkcms(o);
	// do an initial half-step kick
	kick(ppl,np,o,DT/2);

	// initialize the free octree nodes list
	FO = NULL;

	// X init
	disp = XOpenDisplay(NULL);
	if (disp==NULL) {
		fprintf(stderr,"Cannot open display\n");
		return 1;
	}
	scr = DefaultScreen(disp);
	win = XCreateSimpleWindow(disp,RootWindow(disp,scr),
				  10,10,WH,WH,
				  1,WhitePixel(disp,scr),
				  BlackPixel(disp,scr));
	XSelectInput(disp,win,ExposureMask|KeyPressMask);
	XMapWindow(disp,win);
	gc = DefaultGC(disp,scr);
	pm = XCreatePixmap(disp,win,WH,WH,DefaultDepth(disp,scr));
	XFillRectangle(disp,pm,gc,0,0,WH,WH);
	XSetForeground(disp,gc,WhitePixel(disp,scr));

	// event loop
	ts.tv_sec = 0;
	ts.tv_nsec = 1000000000L/FPS;
	for (;;) {
		XNextEvent(disp,&evt);
		if (evt.type==Expose) {
			// drop all pending Exposure events
			while (XCheckMaskEvent(disp,ExposureMask,&evt));
			// advance the system from one frame to another
			if (np>1)
				np = step(ppl,np,o,DT);
			// draw image
			if (!(flags&O_TRACES)) {
				// Clear screen
				XSetForeground(disp,gc,0);
				XFillRectangle(disp,pm,gc,0,0,WH,WH);
			}
			for (n=0; n<np; ++n) {
				p = ppl[n];
				proj(p->x,p->y,p->z,&u,&v);
				proj(p->px,p->py,p->pz,&pu,&pv);
				s = len(p->vx,p->vy,p->vz);
				if (flags&O_EDGES) {
					// draw edges of octrees
					// internal nodes
					XSetForeground(disp,gc,0x808080);
					dredges(disp,pm,gc,o,0);
					// leaf nodes
					XSetForeground(disp,gc,0xFFFFFF);
					dredges(disp,pm,gc,o,1);
				}
				XSetForeground(disp,gc,mag2clr(s));
				if (flags&O_TRACES)
					// draw trace, if asked to
					XDrawLine(disp,pm,gc,pu,pv,u,v);
				XDrawPoint(disp,pm,gc,u,v);
			}
			// show pixmap
			XCopyArea(disp,pm,win,gc,0,0,WH,WH,0,0);
			// wait
			nanosleep(&ts,NULL);
			// send itself an exposure event
			// to display the next frame
			evt.type = Expose;
			XSendEvent(disp,win,False,0,&evt);
		}
		if (evt.type==KeyPress) {
			ks = XLookupKeysym(&evt.xkey,0);
			if (ks==XK_q)
				// quit
				break;
			if (ks==XK_t) {
				// toggle traces
				flags ^= O_TRACES;
			}
			if (ks==XK_e) {
				// toggle edges
				flags ^= O_EDGES;
			}
			if (ks==XK_c) {
				// check conservation laws
				printf("points = %d\n",np);
				printf("steps = %d\n",Stat.steps);
				printf("T+U = %g\n",T+U);
				printf("P = %g\n",len(P_x,P_y,P_z));
				printf("L = %g\n",len(L_x,L_y,L_z));
				printf("<T>+<V>/2 = %g\n", \
				(Tsum+Vsum/2)/(Stat.steps<NAVG?Stat.steps:NAVG));
				printf("\n");
			}
			if (ks==XK_s) {
				// show statistics
				printf("points = %d\n",np);
				printf("steps = %d\n",Stat.steps);
				printf("CM/total = %g%%\n",Stat.acc_cm*100.0/ \
					(Stat.acc_cm+Stat.acc_direct));
				printf("kick time cur/max = %d/%d us\n",
					Stat.cur_t,Stat.max_t);
				printf("leaf nodes cur/max = %d/%d\n", \
					Stat.cur_leaves,Stat.max_leaves);
				printf("internal nodes cur/max = %d/%d\n", \
					Stat.cur_ints,Stat.max_ints);
				printf("height cur/max = %d/%d\n", \
					Stat.cur_h,Stat.max_h);
				printf("maxps cur/max = %d/%d\n", \
					Stat.cur_maxps,Stat.max_maxps);
				printf("reinserts = %d\n",Stat.reins);
				printf("reinsert queue cur/max = %d/%d\n", \
					Stat.cur_rq,Stat.max_rq);
				printf("prunes = %d\n",Stat.prune);
				printf("lifts = %d\n",Stat.lift);
				printf("purges = %d\n",Stat.purge);
				printf("folen = %d\n",Stat.folen);
				printf("reqs malloc/fo = %d/%d\n", \
					Stat.req_malloc,Stat.req_fo);
				printf("\n");
			}
			if (ks==XK_d) {
				// dump state to stderr
				for (n=0; n<np; ++n) {
					p = ppl[n];
					fprintf(stderr,"%g %g %g %g %g %g %g\n",
					p->x,p->y,p->z,p->vx,p->vy,p->vz,p->m);
				}
			}
			if (ks==XK_Left) {
				// 15 degrees left
				Phi = (360+Phi-15)%360;
				printf("Phi=%d, Theta=%d\n",Phi,Theta);
				printf("\n");
				reproj();
				flags &= ~O_TRACES;
			}
			if (ks==XK_Right) {
				// 15 degrees right
				Phi = (360+Phi+15)%360;
				printf("Phi=%d, Theta=%d\n",Phi,Theta);
				printf("\n");
				reproj();
				flags &= ~O_TRACES;
			}
			if (ks==XK_Down) {
				// 15 degrees down
				if (Theta>0) {
					Theta -= 15;
					printf("Phi=%d, Theta=%d\n",Phi,Theta);
					printf("\n");
					reproj();
					flags &= ~O_TRACES;
				}
			}
			if (ks==XK_Up) {
				// 15 degrees up
				if (Theta<180) {
					Theta += 15;
					printf("Phi=%d, Theta=%d\n",Phi,Theta);
					printf("\n");
					reproj();
					flags &= ~O_TRACES;
				}
			}
		}
	}

	// clean up
	XFreePixmap(disp,pm);
	XCloseDisplay(disp);
	return 0;
}
