//
// Geodesics using cartesian representation
//
// Direct problem
//
// 03-08-2018  -  v. 9  - quad-precision
// Ralson
// No intermediate azimuths

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
typedef __float128 quadfloat;

	quadfloat const zer=0.00000000000000000000000000000;
	quadfloat const one=1.00000000000000000000000000000;
	quadfloat const a=6378137.0;
	quadfloat const finv=298.257223563;
	quadfloat const pi=4.*atan(1.);
	quadfloat const ras=648000./pi;	// Conversion from rad to arcsec
	quadfloat me2,a2,b2;
	long int sub,ns;			// number of intervals  (steps)

quadfloat angqua( quadfloat x,  quadfloat y);

//void direct( quadfloat **res, quadfloat x0, quadfloat y0, quadfloat z0, quadfloat a0, quadfloat s0, quadfloat &xf, quadfloat &yf, quadfloat &zf, quadfloat &azf, quadfloat &sf, quadfloat &C0, quadfloat &C, quadfloat &dC, quadfloat &maxdC,  quadfloat &CC, quadfloat &maxCC);

void direct( quadfloat x0, quadfloat y0, quadfloat z0, quadfloat a0, quadfloat s0, quadfloat &xf, quadfloat &yf, quadfloat &zf, quadfloat &az, quadfloat &sf, quadfloat &C0, quadfloat &dCf, quadfloat &CCf);

//void RK4( quadfloat &s, quadfloat s0, quadfloat V0[], quadfloat **res);

void RK4( quadfloat &s, quadfloat s0, quadfloat V0[], quadfloat w[]);

quadfloat huf(quadfloat v[]);

void RFk(quadfloat w[], quadfloat v[], quadfloat sf, quadfloat kc[4][6],  quadfloat h,int m);


int main()
{
	long int i,*id,ng;			// number of geodesics
	int cps=CLOCKS_PER_SEC;
	long int dd=-1,dc=-1,inda[6];
	clock_t t1, t2;
	string nam_res,nam_comp,nam_stat,fil1,fil2;

	quadfloat f,e2;
    quadfloat dab[6]; 
	quadfloat *x0,*y0,*z0,*xe,*ye,*ze,*a0,*ae,*si,*cp,*cq,C0,*dis;
	quadfloat *dx,*dy,*dz,*da,*dr,gmaxdC,gmaxCC;
//    quadfloat *big_unified,**Res;                    // big matrix Res is only useful for detailed statistics
	quadfloat *x,*y,*z,*az,*s,dC,CC;

	cout<<"\n Enter number of intervals : ";
	cin>>sub;
	ns = sub+1;

	cout<<"\n Enter number of sample geodesics : ";
	cin>>ng;
	
	i=sub/100;
	fil1=std::to_string(i);
	i=ng/100;
	fil2=std::to_string(i);	

//    big_unified = (quadfloat *)calloc((ns) * MATR_DIM, sizeof(quadfloat));
//	Res = (quadfloat **)calloc(ns, sizeof(quadfloat *));                    // big matrix Res is only useful for detailed statistics

	id   =  (long int *)calloc(ng,sizeof(long int));
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

//#pragma omp parallel for shared(Res, big_unified)                    // big matrix Res is only useful for detailed statistics
//	for (i = 0; i < ns; i++)
//		Res[i] = &big_unified[i * MATR_DIM];

	f = one/finv;
	e2 = f*(2.*one - f);
	me2 = (one - e2);
	a2 = a*a;
	b2 = a2*me2;

	gmaxdC = -1.;
	gmaxCC = -1.;

	nam_res  = "dirsh_resq-"  + fil1 + "RKp" + fil2 + ".txt";
	nam_comp = "dirsh_compq-" + fil1 + "RKp" + fil2 + ".txt";
	nam_stat = "dirsh_statq-" + fil1 + "RKp" + fil2 + ".txt";
	
	ifstream dat("data.txt");
	if (!dat.is_open())  return 1;
	ofstream des(nam_res);
	ofstream com(nam_comp);
	ofstream sta(nam_stat);
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
//		dKar[i] = strtoflt128(w_s, NULL);
	}
	dat.close();
	
	t2 = clock();
    t1 = t2 - t1;
    sta<<"\n Input  time for  "<<ng<<"  geodesics is  "<<float(t1)/cps<<"  seconds"<<endl;
	sta<<"\n Number of intervals =  "<<sub<<endl;
//  End of input

//  Phase 2  -  Computation of geodesic line
//  Parallel section (split total number of lines) - each independent

#pragma omp parallel for shared(dx,dy,dz,da,dr,s,dis)
	for (i = 0 ; i < ng ; i++)  {		

        C0 = zer;

//		direct(Res,x0[i],y0[i],z0[i],a0[i],si[i],x[i],y[i],z[i],az[i],s[i],C0,C,dC,maxdC,CC,maxCC);
		direct(x0[i],y0[i],z0[i],a0[i],si[i],x[i],y[i],z[i],az[i],s[i],C0,dC,CC);

		dx[i] = x[i] - xe[i];
		dy[i] = y[i] - ye[i];
		dz[i] = z[i] - ze[i];
		da[i] = az[i] -ae[i];
//		delc  = C - cq[i];
		dr[i] = sqrtq(dx[i]*dx[i] + dy[i]*dy[i] + dz[i]*dz[i]);
		dis[i]= s[i] - si[i];

		if (fabsq(CC) > gmaxCC) {
			gmaxCC = fabsq(CC);
			dc = id[i];
		}
		if (fabsq(dC) > gmaxdC) {
			gmaxdC = fabsq(dC);
			dd = id[i];
		}
    }
//  End of parallel section
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
		des<<endl;

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
        des.close();
        
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
		com<<endl;

/*
        quadmath_snprintf(x_s, sizeof(x_s), "%.12Qe", delc);
		com<<setw(30)<<x_s;
		quadmath_snprintf(w_s, sizeof(w_s), "%.12Qe", maxdC);
		com<<setw(30)<<w_s<<endl;
*/
    }

	com.close();
	
	quadmath_snprintf(x_s, sizeof(x_s), "%.5Qe", gmaxdC);
	quadmath_snprintf(y_s, sizeof(y_s), "%.5Qe", gmaxCC);
	sta<<"\n Max (max|dC|) = "<<setw(15)<<x_s<<"  m , found in geodesic "<<dd;
	sta<<"\n Max (max|CC|) = "<<setw(15)<<y_s<<"    , found in geodesic "<<dc<<endl;

//  Statistics section

	for (i = 0 ; i < 6 ; i++)
		dab[i] = -one;

	for (i = 0 ; i < ng ; i++)  {

		if (dab[0] < fabsq(dis[i]))  {
			dab[0] = fabsq(dis[i]);
			inda[0] = i;
		}

		if (dab[1] < fabsq(dr[i]))  {
			dab[1] = fabsq(dr[i]);
			inda[1] = i;
		}

		if (dab[2] < fabsq(dx[i]))  {
			dab[2] = fabsq(dx[i]);
			inda[2] = i;
		}

		if (dab[3] < fabsq(dy[i]))  {
			dab[3] = fabsq(dy[i]);
			inda[3] = i;
		}

		if (dab[4] < fabsq(dz[i]))  {
			dab[4] = fabsq(dz[i]);
			inda[4] = i;
		}

		if (dab[5] < fabsq(da[i]))  {
			dab[5] = fabsq(da[i]);
			inda[5] = i;
		}
	}
	
	sta<<"\n  Property            abs(Max)     (ID) \n";
	quadmath_snprintf(x_s, sizeof(x_s), "%.12Qe", dab[2]);
	sta<<"\n max(|dx|)  "<<setw(25)<<x_s<<setw(9)<<inda[2];
	quadmath_snprintf(y_s, sizeof(y_s), "%.12Qe", dab[3]);
	sta<<"\n max(|dy|)  "<<setw(25)<<y_s<<setw(9)<<inda[3];
	quadmath_snprintf(z_s, sizeof(z_s), "%.12Qe", dab[4]);
	sta<<"\n max(|dz|)  "<<setw(25)<<z_s<<setw(9)<<inda[4];
	quadmath_snprintf(w_s, sizeof(w_s), "%.12Qe", dab[5]);
	sta<<"\n max(|daz|)('')"<<setw(25)<<w_s<<setw(9)<<inda[5];
	quadmath_snprintf(z_s, sizeof(z_s), "%.12Qe", dab[0]);
	sta<<"\n max(|ds|)  "<<setw(25)<<z_s<<setw(9)<<inda[0];
	quadmath_snprintf(y_s, sizeof(y_s), "%.12Qe", dab[1]);
	sta<<"\n max(|dr|)  "<<setw(25)<<y_s<<setw(9)<<inda[1];
	sta<<endl;

	t2 = clock() - t1;
	sta<<"\n Output time for  "<<ng<<"  geodesics is  "<<float(t2)/cps<<"  seconds"<<endl;
	
//  End of output
    
    sta.close();
    
    free (id); free (x0); free (y0); free (z0); free (si);
    free (a0); free (xe); free (ye); free (ze); free (s);
    free (ae); free (x); free (y); free (z); free (cp);
    free (dis); free (dx); free (dy); free (dz); free (cq);
    free (az); free (da); free (dr);
//	free big_unified; free Res;

	return 0;
}


//void direct( quadfloat **res, quadfloat x0, quadfloat y0, quadfloat z0, quadfloat a0, quadfloat s0, quadfloat &xf, quadfloat &yf, quadfloat &zf, quadfloat &azf, quadfloat &sf, quadfloat &C0, quadfloat &Cf, quadfloat &dCf, quadfloat &maxdC,  quadfloat &CCf, quadfloat &maxCC)
void direct( quadfloat x0, quadfloat y0, quadfloat z0, quadfloat a0, quadfloat s0, quadfloat &xf, quadfloat &yf, quadfloat &zf, quadfloat &az, quadfloat &sf, quadfloat &C0, quadfloat &dCf, quadfloat &CCf)
{
	unsigned long i;
	quadfloat v0[6],n[3],p[3],q[3],pf[6];
	quadfloat s,C,dC,CC,xp,yp,zp;
	quadfloat x0p,y0p,z0p,h0,r0,tol=1.0e-09;
//	quadfloat res[ns][6];
	quadfloat H1,n1,n2,n3,r,PP,QQ;
	quadfloat p1,p2,p3,q1,q2,q3;
//    quadfloat **res;

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

//	RK4(s,s0,v0,res);                   // big matrix res is only useful for detailed statistics
	RK4(s,s0,v0,pf);                     // matrix pf is last position/velocity vector

//	maxdC = -1.;
//	maxCC = -1.;

/*    
	for (i = 0; i < ns ; i++) {			// may proceed in parallel (split total number in parts) OR omit entirely     -    big matrix res is only useful for detailed statistics

		C  = res[i][0]*res[i][3] - res[i][2]*res[i][1];
		dC = C - C0;
		CC = (((res[i][0]*res[i][0] + res[i][2]*res[i][2])/a2) + ((res[i][4]*res[i][4])/b2)) - one;

		if (fabsq(dC) > maxdC)  maxdC = fabsq(dC);
		if (fabsq(CC) > maxCC)  maxCC = fabsq(CC);
	}
*/

	xf = pf[0];
	yf = pf[2];
	zf = pf[4];
	xp = pf[1];
	yp = pf[3];
	zp = pf[5];
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

//	azf = az;
	dCf  = xf*yp-yf*xp-C0;
	CCf = (xf*xf + yf*yf)/a2 +  zf*zf/b2  - one;
;
}


//void RK4( quadfloat &s, quadfloat s0, quadfloat V0[], quadfloat **res)                   // big matrix res is only useful for detailed statistics
void RK4( quadfloat &s, quadfloat s0, quadfloat V0[], quadfloat w[])     
{
	quadfloat v[6],sf;   
	quadfloat kc[4][6],h;
	unsigned long i;
	double cof[4];
    int j,m;
    
 	h = s0/(ns - 1);
    
	 // Part 1 - initialize
    
	for (j = 0 ; j < 6 ; j++)  { 	
		w[j] = V0[j];
	}
	
	cof[0]=0.5;
	cof[1]=0.5;
	cof[2]=1.0;
	cof[3]=1.0;
	
	for (i = 0 ; i < ns-1 ; i++) {		// proceed along geodesic

       	for (j=0;j<6;j++)
			v[j]=w[j];

     //  Part 2 - evaluation section - 6 independent variables

// ----------------  this section applies Runge-Kutta 4th order
        
		for (m = 0; m < 4 ; m++)  {      
			sf=huf(v);  
			RFk(w,v,sf,kc,h,m);
            for (j=0;j<6;j++)
                v[j]=w[j]+cof[m]*kc[m][j];
		}
         
         for (j = 0 ; j < 6 ; j++)  {
 			w[j] = w[j] + (kc[0][j] + 2.*kc[1][j] + 2.*kc[2][j] + kc[3][j])/6.;
  		}	

//  -------------------------  end of RK-4
        
	}
	
	s = h*(ns-1);
}
       
//---------------------------

quadfloat huf(quadfloat v[])
{
	quadfloat H1,H2;
    H1 = v[0]*v[0] + v[2]*v[2] + (v[4]*v[4]/(me2*me2));
    H2 = v[1]*v[1] + v[3]*v[3] + (v[5]*v[5]/me2);
    return -(H2/H1);
}

void RFk(quadfloat w[], quadfloat v[], quadfloat sf, quadfloat kc[4][6],  quadfloat h,int m)	
{
	kc[m][0]=h*v[1];
	kc[m][1]=h*sf*v[0];
	kc[m][2]=h*v[3];
	kc[m][3]=h*sf*v[2];
	kc[m][4]=h*v[5];
	kc[m][5]=h*sf*v[4]/me2;
}
//------------------------------

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

