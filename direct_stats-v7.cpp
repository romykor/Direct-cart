//
// Geodesics using cartesian representation
//
// Statistics of Direct problem
//
// 25-07-2018  -  v. 7
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;

	double const pi=4.*atan(1.);
	double const ras=648000./pi;	// Conversion from rad to arcsec


int main()
{
	int j,k,ng,inda[6];
	double dab[6];

	cout<<"\n Enter number of sample geodesics : ";
	cin>>ng;
	ng=ng+1;

	double dx,dy,dz,daz,dr,ds;

	ifstream inf2("direct_compq-100opt10.txt");
	ofstream res("direct_statq-100opt10.txt",ios::app);

	for (j=0;j<6;j++)
		dab[j]=-1.;

	for (j=1;j<ng;j++)  {
		inf2>>k>>dx>>dy>>dz>>daz>>ds>>dr;

		if (dab[0]<fabs(ds))  {
			dab[0]=fabs(ds);
			inda[0]=k;
		}

		if (dab[1]<fabs(dr))  {
			dab[1]=fabs(dr);
			inda[1]=k;
		}

		if (dab[2]<fabs(dx))  {
			dab[2]=fabs(dx);
			inda[2]=k;
		}

		if (dab[3]<fabs(dy))  {
			dab[3]=fabs(dy);
			inda[3]=k;
		}

		if (dab[4]<fabs(dz))  {
			dab[4]=fabs(dz);
			inda[4]=k;
		}

		if (dab[5]<fabs(daz))  {
			dab[5]=fabs(daz);
			inda[5]=k;
		}
	}

	res<<scientific<<setprecision(4)<<endl;
	res<<"\n  Property            abs(Max)     (ID) \n";
	res<<"\n max(|dx|)  "<<setw(20)<<dab[2]<<setw(9)<<inda[2];
	res<<"\n max(|dy|)  "<<setw(20)<<dab[3]<<setw(9)<<inda[3];
	res<<"\n max(|dz|)  "<<setw(20)<<dab[4]<<setw(9)<<inda[4];
	res<<"\n max(|daz|)('')"<<setw(16)<<dab[5]<<setw(9)<<inda[5];
	res<<"\n max(|ds|)  "<<setw(20)<<dab[0]<<setw(9)<<inda[0];
	res<<"\n max(|dr|)  "<<setw(20)<<dab[1]<<setw(9)<<inda[1];
	res<<endl;

	return 0;
}
