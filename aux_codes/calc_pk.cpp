/*

pk_G(k) = (4 pi) \int r^2 dr xi_G(r) j0(kr)

where j0(kr) = sin(kr) / (kr)

*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include "../compute_pkG/spline_array.h"
#include "../compute_pkG/integration_modified.cpp"

using namespace std;

#define PI 3.14159265358979

double sinc(double x);
double int_pkG_func(double r,int dim,double params[]);

int ndata;
double *r1d,*xi1d,*xi1d2;

int main(){

	string xifname;
	cout << "Please enter the input correlation function file:" << endl;
	cin >> xifname;
	int ncolumn;
	cout << "Please enter the number of columns of the correlation function file:" << endl;
	cin >> ncolumn;
	int rc,xic;
	cout << "Please enter the column number (staring from 1) for r and xi(r) of the correlation function file:" << endl;
	cin >> rc >> xic;
	double bias;
	cout << "Please enter the linear bias (please enter 1.0 for cross power spectrum from cross correlation function):" << endl;
	cin >> bias;
	double rmax;
	cout << "Please enter rmax for integration:" << endl;
	cin >> rmax;

	double frac = bias*bias;

	ndata = 0;
	{
		double temp1;
		string temp2;
		ifstream xiin(xifname.c_str());
		while (xiin >> temp1){
			getline(xiin,temp2);
			ndata ++;
		}
		xiin.close();
	}

	r1d = new double [ndata];
	xi1d = new double [ndata];
	xi1d2 = new double [ndata];
	{
		ifstream xiin(xifname.c_str());
		double *temp = new double [ncolumn*ndata];
		for (int i=0;i<ndata;i++){
			for (int j=0;j<ncolumn;j++) xiin >> temp[j+ncolumn*i];
			r1d[i] = temp[rc-1+ncolumn*i];
			xi1d[i] = temp[xic-1+ncolumn*i]*frac;
		}
		xiin.close();
	}
	spline_arr(r1d,xi1d,ndata,1e30,1e30,xi1d2);

	integration romb(1.e-8);

	int dim = 1;
	double params[dim];

	ostringstream cr,cb,cp;
	cr << rmax;
	cb << bias;

	string pkGfname = "pk_rmax" + cr.str() + "_b" + cb.str() + ".dat";
	ofstream pkGout(pkGfname.c_str());
	pkGout.precision(10);
	for (int i=0;i<=429;i++){
		double k = 0.01*double(i)-3;
		k = pow(10.,k);
		params[0] = k;
		double pkG = romb.rombint(int_pkG_func,0,rmax,dim,params);
		pkGout << k << "\t" << 4.*PI*pkG << endl;
	}
	pkGout.close();

	return 0;
}

double int_pkG_func(double r,int dim,double params[]){

	double k = params[0];
	double xi;
	splint_arr(r1d,xi1d,xi1d2,ndata,r,&xi);

	return r*r*xi*sinc(k*r);
}

double sinc(double x){
	if (abs(x)<=1e-2){
		double x2 = x*x;
		return 1.-x2/6.+x2*x2/120.;
	}
	return sin(x)/x;
}
