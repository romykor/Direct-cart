//
// Geodesics using cartesian representation
//
// Direct problem
//
// 25-07-2018  -  v. 8  - quad-precision
// attempt to parallelize
// less arrays - No intermediate azimuths

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <cstring>
#include <quadmath.h>
#include <omp.h>

using namespace std;

#define MATR_DIM	6
#define STR_SIZE	50

#define tID		omp_get_thread_num()
typedef __float128 quadfloat;

	quadfloat const zer=0.00000000000000000000000000000;
	quadfloat const one=1.00000000000000000000000000000;
	quadfloat const a=6378137.0;
	quadfloat const finv=298.257223563;
	quadfloat const pi=4.*atan(1.);
	quadfloat const ras=648000./pi;	// Conversion from rad to arcsec
	quadfloat me2,a2,b2;
	long int sub,ns;			// number of intervals
//	quadfloat Res[][MATR_DIM];


quadfloat angqua( quadfloat x,  quadfloat y);

void direct( quadfloat **res, quadfloat x0, quadfloat y0,
		quadfloat z0, quadfloat a0, quadfloat s0, quadfloat &xf,
		quadfloat &yf, quadfloat &zf, quadfloat &azf,
		quadfloat &sf, quadfloat &C0, quadfloat &C, quadfloat &dC,
		quadfloat &maxdC,  quadfloat &CC, quadfloat &maxCC);

void RK4( quadfloat &s, quadfloat s0, quadfloat V0[], quadfloat **res);

int main()
{
	long int i,*id,ng;			// number of geodesics
	int cps = CLOCKS_PER_SEC;
	long int dd = -1, dc = -1;
	clock_t t1, t2;

	int threads = omp_get_max_threads();

	cout << "Threads enabled: " << threads << endl;


	quadfloat f,e2;
	quadfloat delc,maxdC,maxCC,gmaxdC,gmaxCC,C,dC,C0,CC;
	quadfloat thread_gmaxdC[threads], thread_gmaxCC[threads] , thread_dd[threads], thread_dc[threads];
    quadfloat *big_unified;
	quadfloat *x0,*y0,*z0,*xe,*ye,*ze,*a0,*ae,*si,*cp,*cq,*dKar,*dis;
	quadfloat *dx,*dy,*dz,*da,*dr,*x,*y,*z,*az,*s;
    quadfloat **Res;

	cout<<"\n Enter number of intervals : ";
	cin>>sub;
	ns = sub+1;

	cout<<"\n Enter number of sample geodesics : ";
	cin>>ng;

    big_unified = (quadfloat *)calloc((ns) * MATR_DIM, sizeof(quadfloat));
	Res = (quadfloat **)calloc(ns, sizeof(quadfloat *));

	id    = (long int *)calloc(ng,sizeof(long int));
    x0   = (quadfloat *)calloc(ng, sizeof(quadfloat));
	y0   = (quadfloat *)calloc(ng, sizeof(quadfloat));
	z0   = (quadfloat *)calloc(ng, sizeof(quadfloat));
	xe   = (quadfloat *)calloc(ng, sizeof(quadfloat));
	ye   = (quadfloat *)calloc(ng, sizeof(quadfloat));
	ze   = (quadfloat *)calloc(ng, sizeof(quadfloat));
	a0   = (quadfloat *)calloc(ng, sizeof(quadfloat));
	ae   = (quadfloat *)calloc(ng, sizeof(quadfloat));
	si   = (quadfloat *)calloc(ng, sizeof(quadfloat));
	cp   = (quadfloat *)calloc(ng, sizeof(quadfloat));
	cq   = (quadfloat *)calloc(ng, sizeof(quadfloat));
	dKar = (quadfloat *)calloc(ng, sizeof(quadfloat));
	x    = (quadfloat *)calloc(ng, sizeof(quadfloat));
	y    = (quadfloat *)calloc(ng, sizeof(quadfloat));
	z    = (quadfloat *)calloc(ng, sizeof(quadfloat));
	az   = (quadfloat *)calloc(ng, sizeof(quadfloat));
	s    = (quadfloat *)calloc(ng, sizeof(quadfloat));
	dis  = (quadfloat *)calloc(ng, sizeof(quadfloat));
	dx   = (quadfloat *)calloc(ng, sizeof(quadfloat));
	dy   = (quadfloat *)calloc(ng, sizeof(quadfloat));
	dz   = (quadfloat *)calloc(ng, sizeof(quadfloat));
	da   = (quadfloat *)calloc(ng, sizeof(quadfloat));
	dr   = (quadfloat *)calloc(ng, sizeof(quadfloat));

#pragma omp parallel for shared(Res, big_unified)
	for (i = 0; i < ns; i++)
		Res[i] = &big_unified[i * MATR_DIM];

	f = one/finv;
	e2 = f*(2.*one - f);
	me2 = (one - e2);
	a2 = a*a;
	b2 = a2*me2;

	gmaxdC = -1.;
	gmaxCC = -1.;

	memset(thread_gmaxCC, 0, sizeof(thread_gmaxCC));
	memset(thread_gmaxdC, 0, sizeof(thread_gmaxdC));
	memset(thread_dc,	  0, sizeof(thread_dc));
	memset(thread_dd,	  0, sizeof(thread_dd));

	ifstream dat("data.txt");
	if (!dat.is_open())  return 1;
	ofstream des("direct_resq-opt.txt");
	ofstream com("direct_compq-opt.txt");
	ofstream sta("direct_statq-opt.txt");
	des<<endl;
	com<<endl;
	sta<<scientific<<setprecision(4)<<endl;

	char x_s[STR_SIZE], y_s[STR_SIZE], z_s[STR_SIZE], w_s[STR_SIZE];
	memset(x_s, 0, STR_SIZE * sizeof(char));
	memset(y_s, 0, STR_SIZE * sizeof(char));
	memset(z_s, 0, STR_SIZE * sizeof(char));
	memset(w_s, 0, STR_SIZE * sizeof(char));

	t1 = clock();

//  Phase 1  -  Data input

	for (i = 0; i < ng ; i++)  {
		dat>>id[i]>>x_s>>y_s>>z_s>>w_s;
		x0[i] = strtoflt128(x_s, NULL);
		y0[i] = strtoflt128(y_s, NULL);
		z0[i] = strtoflt128(z_s, NULL);
		a0[i] = strtoflt128(w_s, NULL);

		dat>>x_s>>y_s>>z_s>>w_s;
		xe[i] = strtoflt128(x_s, NULL);
		ye[i] = strtoflt128(y_s, NULL);
		ze[i] = strtoflt128(z_s, NULL);
		ae[i] = strtoflt128(w_s, NULL);

		dat>>x_s>>y_s>>z_s>>w_s;
		si[i]   = strtoflt128(x_s, NULL);
		cp[i]   = strtoflt128(y_s, NULL);
		cq[i]   = strtoflt128(z_s, NULL);
		dKar[i] = strtoflt128(w_s, NULL);
	}
	dat.close();

	t2=clock();
    t1=t2-t1;
    sta<<"\n Input  time for  "<<ng<<"  geodesics is  "<<float(t1)/cps<<"  seconds"<<endl;
	sta<<"\n Number of intervals =  "<<sub<<endl;
//  End of input

//  Phase 2  -  Computation of geodesic line
//  Parallel section (split total number of lines) - each independent

#pragma omp parallel for firstprivate(maxdC, maxCC, C0, CC, Res, delc)
	for (i = 0 ; i < ng ; i++)  {
		maxdC = -1.;
		maxCC = -1.;
		C0 = zer;

		direct(Res,x0[i],y0[i],z0[i],a0[i],si[i],x[i],y[i],z[i],az[i],s[i],C0,C,dC,maxdC,CC,maxCC);

		dx[i] = x[i] - xe[i];
		dy[i] = y[i] - ye[i];
		dz[i] = z[i] - ze[i];
		da[i] = az[i] -ae[i];
		delc  = C - cq[i];
		dr[i] = sqrtq(dx[i]*dx[i] + dy[i]*dy[i] + dz[i]*dz[i]);
		dis[i]= s[i] - si[i];

		if (fabsq(CC) > thread_gmaxCC[tID]) {
			thread_gmaxCC[tID] = fabsq(CC);
			thread_dc[tID] = id[i];
		}
		if (fabsq(dC) > thread_gmaxdC[tID]) {
			thread_gmaxdC[tID] = fabsq(dC);
			thread_dd[tID] = id[i];
		}
    }
//  End of parallel section

//  Now find which were the max values and the corresponding id's

	for (i = 0; i < threads; i++) {
		if (thread_gmaxCC[i] > gmaxCC) {
			gmaxCC = thread_gmaxCC[i];
			dc = thread_dc[i];
		}
		if (thread_gmaxdC[i] > gmaxdC) {
			gmaxdC = thread_gmaxdC[i];
			dd = thread_dd[i];
		}
	}

//  End of computation
  	t1 = clock();
    t2 = t1 - t2;
    sta<<"\n Computing  time for  "<<ng<<"  geodesics is  "<<float(t2)/cps<<"  seconds"<<endl;

//  Phase 3  -  Output of results

    for (i = 0 ; i < ng ; i++)  {
	    des<<setw(9)<<id[i];
		quadmath_snprintf(x_s, sizeof(x_s), "%.12Qe", x[i]);
		des<<setw(30)<<x_s;
		quadmath_snprintf(y_s, sizeof(y_s), "%.12Qe", y[i]);
		des<<setw(30)<<y_s;
		quadmath_snprintf(z_s, sizeof(z_s), "%.12Qe", z[i]);
		des<<setw(30)<<z_s;
		quadmath_snprintf(w_s, sizeof(w_s), "%.12Qe", az[i]);
		des<<setw(30)<<w_s;
		quadmath_snprintf(x_s, sizeof(x_s), "%.12Qe", s[i]);
		des<<setw(30)<<x_s;

/*
		quadmath_snprintf(x_s, sizeof(x_s), "%.12Qe", C0);
		des<<setw(30)<<x_s;
		quadmath_snprintf(y_s, sizeof(y_s), "%.12Qe", C);
		des<<setw(30)<<y_s;
		quadmath_snprintf(z_s, sizeof(z_s), "%.12Qe", dC);
		des<<setw(30)<<z_s;
		quadmath_snprintf(w_s, sizeof(w_s), "%.12Qe", maxdC);
		des<<setw(30)<<w_s;
		quadmath_snprintf(x_s, sizeof(x_s), "%.12Qe", CC);
		des<<setw(30)<<x_s;
		quadmath_snprintf(y_s, sizeof(y_s), "%.12Qe", maxCC);
		des<<setw(30)<<y_s<<endl;
*/

		com<<setw(9)<<id[i];
		quadmath_snprintf(x_s, sizeof(x_s), "%.12Qe", dx[i]);
		com<<setw(30)<<x_s;
		quadmath_snprintf(y_s, sizeof(y_s), "%.12Qe", dy[i]);
		com<<setw(30)<<y_s;
		quadmath_snprintf(z_s, sizeof(z_s), "%.12Qe", dz[i]);
		com<<setw(30)<<z_s;
		quadmath_snprintf(w_s, sizeof(w_s), "%.12Qe", da[i]*ras);
		com<<setw(30)<<w_s;
 		quadmath_snprintf(z_s, sizeof(z_s), "%.12Qe", dis[i]);
		com<<setw(30)<<z_s;
		quadmath_snprintf(y_s, sizeof(y_s), "%.12Qe", dr[i]);
		com<<setw(30)<<y_s;

/*
        quadmath_snprintf(x_s, sizeof(x_s), "%.12Qe", delc);
		com<<setw(30)<<x_s;
		quadmath_snprintf(w_s, sizeof(w_s), "%.12Qe", maxdC);
		com<<setw(30)<<w_s<<endl;
*/
    }
//  End of output

	t2 = clock() - t1;

	quadmath_snprintf(x_s, sizeof(x_s), "%.5Qe", gmaxdC);
	quadmath_snprintf(y_s, sizeof(y_s), "%.5Qe", gmaxCC);
	sta<<"\n Max (max|dC|) = "<<setw(15)<<x_s<<"  m , found in geodesic "<<dd;
	sta<<"\n Max (max|CC|) = "<<setw(15)<<y_s<<"    , found in geodesic "<<dc<<endl;
	sta<<"\n Output time for  "<<ng<<"  geodesics is  "<<float(t2)/cps<<"  seconds"<<endl;

	return 0;
}


void direct( quadfloat **res, quadfloat x0, quadfloat y0,
		quadfloat z0, quadfloat a0, quadfloat s0, quadfloat &xf,
		quadfloat &yf, quadfloat &zf, quadfloat &azf, quadfloat &sf,
		quadfloat &C0, quadfloat &Cf, quadfloat &dCf, quadfloat &maxdC,
		quadfloat &CCf, quadfloat &maxCC)
{
	unsigned long i;
	quadfloat v0[6],n[3],p[3],q[3];
	quadfloat s,az,C,dC,CC,xp,yp,zp;
	quadfloat x0p,y0p,z0p,h0,r0,tol=1.0e-09;
//	quadfloat res[ns][6];
	quadfloat H1,n1,n2,n3,r,PP,QQ;
	quadfloat p1,p2,p3,q1,q2,q3;

	r0  = x0*x0+y0*y0;
	h0  = sqrtq(r0+(z0*z0)/(me2*me2));
	n[0]= x0/h0;
	n[1]= y0/h0;
	n[2]= z0/(me2*h0);
	r0  = sqrtq(r0);
	s   = zer;

	if  ((r0 < tol) && (fabsq(z0 - a*sqrtq(me2)) < tol))  {
		p[0] = zer;
		p[1] = one;
		p[2] = zer;
		q[0] = -one;
		q[1] = zer;
		q[2] = zer;
	}
	else if  ((r0<tol) && (fabsq(z0+a*sqrtq(me2))<tol))  {
		p[0] = zer;
		p[1] = one;
		p[2] = zer;
		q[0] = one;
		q[1] = zer;
		q[2] = zer;
	}
	else  {
		p[0] = -y0/r0;
		p[1] = x0/r0;
		p[2] = zer;
		q[0] = -n[2]*p[1];
		q[1] = n[2]*p[0];
		q[2] = n[0]*p[1]-n[1]*p[0];
	}

	C0  = r0*sinq(a0);
	x0p = p[0]*sinq(a0)+q[0]*cosq(a0);
	y0p = p[1]*sinq(a0)+q[1]*cosq(a0);
	z0p = p[2]*sinq(a0)+q[2]*cosq(a0);

	v0[0] = x0;
	v0[1] = x0p;
	v0[2] = y0;
	v0[3] = y0p;
	v0[4] = z0;
	v0[5] = z0p;

	RK4(s,s0,v0,res);

	maxdC = -1.;
	maxCC = -1.;

	for (i = 0; i < ns ; i++) {			// may proceed in parallel (split total number in parts) OR omit entirely

		C  = res[i][0]*res[i][3] - res[i][2]*res[i][1];
		dC = C - C0;
		CC = (((res[i][0]*res[i][0] + res[i][2]*res[i][2])/a2) + ((res[i][4]*res[i][4])/b2)) - one;

		if (fabsq(dC) > maxdC)  maxdC = fabsq(dC);
		if (fabsq(CC) > maxCC)  maxCC = fabsq(CC);
	}

	xf = res[ns-1][0];
	yf = res[ns-1][2];
	zf = res[ns-1][4];
	xp = res[ns-1][1];
	yp = res[ns-1][3];
	zp = res[ns-1][5];
	sf = s;

	r  = xf*xf + yf*yf;
	H1 = r + (zf*zf)/(me2*me2);
	r  = sqrtq(r);
	H1 = sqrtq(H1);

	n1 = xf/H1;
	n2 = yf/H1;
	n3 = zf/(me2*H1);

	p1 = -yf/r;
	p2 = xf/r;
	p3 = zer;

	q1 = -n3*p2;
	q2 = n3*p1;
	q3 = n1*p2 - n2*p1;

	PP = p1*xp + p2*yp + p3*zp;
	QQ = q1*xp + q2*yp + q3*zp;

	az = angqua(QQ,PP);

	azf = az;
	Cf  = C;
	dCf = dC;
	CCf = CC;
}


void RK4( quadfloat &s, quadfloat s0, quadfloat V0[], quadfloat **res)
{
	quadfloat h,k1[6],k2[6],k3[6],k4[6];
	quadfloat w[6],v[6],H1,H2,sf;
	unsigned long i,j;

	h = s0/(ns - 1);

	for (j = 0 ; j < 6 ; j++)  { 			// initialization
		w[j] = V0[j];
		res[0][j] = V0[j];
	}

	for (i = 0 ; i < ns-1 ; i++) {		// proceed along geodesic

		for (j = 0 ; j < 6 ; j++)			// initialization
			v[j] = w[j];

//*
// ----------------  this section applies Runge-Kutta 4th order

	 // part 1
		H1 = v[0]*v[0] + v[2]*v[2] + (v[4]*v[4]/(me2*me2));
		H2 = v[1]*v[1] + v[3]*v[3] + (v[5]*v[5]/me2);
		sf = -(H2/H1);

		k1[0] = h*v[1];
		k1[1] = h*sf*v[0];
		k1[2] = h*v[3];
		k1[3] = h*sf*v[2];
		k1[4] = h*v[5];
		k1[5] = h*sf*v[4]/me2;
		for (j = 0 ; j < 6 ; j++)
			v[j] = w[j] + 0.5*k1[j];

	 // part 2
		H1 = v[0]*v[0] + v[2]*v[2] + (v[4]*v[4]/(me2*me2));
		H2 = v[1]*v[1] + v[3]*v[3] + (v[5]*v[5]/me2);
		sf = -(H2/H1);

		k2[0] = h*v[1];
		k2[1] = h*sf*v[0];
		k2[2] = h*v[3];
		k2[3] = h*sf*v[2];
		k2[4] = h*v[5];
		k2[5] = h*sf*v[4]/me2;
		for (j = 0 ; j < 6 ; j++)
			v[j] = w[j] + 0.5*k2[j];

	 // part 3
		H1 = v[0]*v[0] + v[2]*v[2] + (v[4]*v[4]/(me2*me2));
		H2 = v[1]*v[1] + v[3]*v[3] + (v[5]*v[5]/me2);
		sf = -(H2/H1);

		k3[0] = h*v[1];
		k3[1] = h*sf*v[0];
		k3[2] = h*v[3];
		k3[3] = h*sf*v[2];
		k3[4] = h*v[5];
		k3[5] = h*sf*v[4]/me2;
		for (j = 0 ; j < 6 ; j++)
			v[j] = w[j] + k3[j];

	 // part 4
		H1 = v[0]*v[0] + v[2]*v[2] + (v[4]*v[4]/(me2*me2));
		H2 = v[1]*v[1] + v[3]*v[3] + (v[5]*v[5]/me2);
		sf = -(H2/H1);

		k4[0] = h*v[1];
		k4[1] = h*sf*v[0];
		k4[2] = h*v[3];
		k4[3] = h*sf*v[2];
		k4[4] = h*v[5];
		k4[5] = h*sf*v[4]/me2;
		for (j = 0 ; j < 6 ; j++)		//  not useful any more
			v[j] = w[j] + k3[j];		//  not useful any more


	 // summation and residuals
		for (j = 0 ; j < 6 ; j++)  {
			w[j] = w[j] + (k1[j] + 2.*k2[j] + 2.*k3[j] + k4[j])/6.;
			res[i+1][j] = w[j];
		}
//  -------------------------  end of RK-4
//*/

// ----------------  this section applies Ralson



//  -------------------------  end of Ralson

	}
	s = h*(ns-1);
}


 quadfloat angqua( quadfloat x, quadfloat y)
{
	quadfloat a, g;

	if ((x == zer) && (y == zer))  g = zer;
	else  {
		if (fabsq(y) < fabsq(x)) {
			a = atanq(fabsq(y)/fabsq(x));
			if (x > zer)  {
				if (y > zer)  g = a;
				if (y < zer)  g = 2.*pi-a;
				if (y == zer) g = zer;
			}
			if (x < zer)  {
				if (y > zer)  g = pi-a;
				if (y < zer)  g = pi+a;
				if (y == zer) g = pi;
			}
		}
		if (fabsq(y) > fabsq(x)) {
			a = atanq(fabsq(x)/fabsq(y));
			if (y > zer)  {
				if (x > zer)  g = 0.5*pi-a;
				if (x < zer)  g = 0.5*pi+a;
				if (x == zer) g = 0.5*pi;
			}
			if (y < zer)  {
				if (x > zer)  g = 1.5*pi+a;
				if (x < zer)  g = 1.5*pi-a;
				if (x == zer) g = 1.5*pi;
			}
		}
	}
	return g;
}

