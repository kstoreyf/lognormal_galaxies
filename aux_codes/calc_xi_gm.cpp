/*

xiGgm(r) = (1/2pi^2) \int k^2 dk sqrt(P_Gmm(k)*P_Ggg(k)) j0(kr)

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
double int_xiG_func(double r,int dim,double params[]);

int ndata;
double *k1d,*pk1d,*pk1d2;

int main(){
	
	cout << "This code assumes that the matter and galaxy Gaussian power spectra are defined at the same k values! Please use interpolation if this is not the case." << endl;
	string pkGfname;
	cout << "Please enter the name of galaxy Gaussian power spectrum file:" << endl;
	cin >> pkGfname;
	string mpkGfname;
	cout << "Please enter the name of matter Gaussian power spectrum file:" << endl;
	cin >> mpkGfname;
	int ncolumn;
	cout << "Please enter the number of columns of the power spectra files:" << endl;
	cin >> ncolumn;
	int rc,xic;
	cout << "Please enter the column number (staring from 1) for k and P(k) of the power spectra files:" << endl;
	cin >> rc >> xic;
	/*double bias;
	cout << "Please enter the linear bias:" << endl;
	cin >> bias;*/
	double kmax;
	cout << "Please enter kmax for integration:" << endl;
	cin >> kmax;

	ndata = 0;
	{
		double temp1;
		string temp2;
		ifstream pkin(pkGfname.c_str());
		while (pkin >> temp1){
			getline(pkin,temp2);
			ndata ++;
		}
		pkin.close();
	}

	k1d = new double [ndata];
	pk1d = new double [ndata];
	pk1d2 = new double [ndata];
	{
		ifstream pkin(pkGfname.c_str()), mpkin(mpkGfname.c_str());
		double *temp = new double [ncolumn*ndata], *temp2 = new double[ncolumn*ndata];
		for (int i=0;i<ndata;i++){
			for (int j=0;j<ncolumn;j++)
			{
				pkin >> temp[j+ncolumn*i];
				mpkin >> temp2[j+ncolumn*i];
			}
			k1d[i] = temp[rc-1+ncolumn*i];
			pk1d[i] = sqrt(temp[xic-1+ncolumn*i]*temp2[xic-1+ncolumn*i]);
		}
		delete[] temp;
		delete[] temp2;
		pkin.close();
		mpkin.close();
	}
	spline_arr(k1d,pk1d,ndata,1e30,1e30,pk1d2);

	integration romb(1.e-8);

	int dim = 1;
	double params[dim];

	ostringstream cr,cp;
	cr << kmax;
	
	double pi2 = PI*PI;
	string xiGfname = "xi_gm_kmax" + cr.str() + ".dat";
	ofstream xiGout(xiGfname.c_str());
	xiGout.precision(10);
	for (int i=0;i<=529;i++){
		double r = 0.01*double(i)-2.;
		r = pow(10.,r);
		params[0] = r;
		double xiG = romb.rombint(int_xiG_func,0,kmax,dim,params);
		xiG = xiG/2./pi2;
		double xi_gm = exp(xiG)-1.;
		xiGout << r << "\t" << xi_gm << endl;
	}
	xiGout.close();

	return 0;
}

double int_xiG_func(double k,int dim,double params[]){

	double r = params[0];
	double pk;
	splint_arr(k1d,pk1d,pk1d2,ndata,k,&pk);
	return k*k*pk*sinc(k*r);
}

double sinc(double x){
	if (abs(x)<=1e-2){
		double x2 = x*x;
		return 1.-x2/6.+x2*x2/120.;
	}
	return sin(x)/x;
}
