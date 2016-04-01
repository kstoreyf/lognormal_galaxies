#include <iostream>
#include <fstream>
#include <cmath>
#include <cassert>
#include <string>
#include "../generate_Poisson/spline.h"

using namespace std;

int main(int argc, char *argv[])
{
	string mpkfname;
	cout << "Please enter name of matter power spectrum file [1st: k in units of h/Mpc; 2nd: P(k) in units of Mpc^3/h^3]" << endl;
	cin >> mpkfname;
	int nmpk;
	cout << "Please enter the number of lines in the matter Pk file above" << endl;
	cin >> nmpk;
	string pkgmfname;
	cout << "Please enter name of cross power spectrum file [1st: k in units of h/Mpc; 2nd: P(k) in units of Mpc^3/h^3]" << endl;
	cin >> pkgmfname;
	int npkgm;
	cout << "Please enter the number of lines in the cross Pk file above" << endl;
	cin >> npkgm;
	double bias;
	cout << "Please enter the linear bias" << endl;
	cin >> bias;
	string ffname;         
	cout << "Please enter the logarithmic growth rate file [1st: k in units of h/Mpc; 2nd: f(k,z)]:" << endl;  
	cin >> ffname;       
	int ndataf;          
	cout << "Please enter the number of lines in the logarithmic growth rate file above" << endl; 
	cin >> ndataf;
	
	interpol mpk(mpkfname,nmpk,2,1,2,true,true);
	interpol cpk(pkgmfname,npkgm,2,1,2,true,true);
	interpol f(ffname, ndataf, 2, 1, 2, true, true);	

	double b2 = bias*bias;
	ofstream wf("pk_kaiser.dat");
	for (int i=0;i<=429;i++){
                double k = 0.01*double(i)-3.;
                k = pow(10.,k);
		double growth = f.value(k), growth2 = growth*growth, mpkval = mpk.value(k), cpkval = cpk.value(k);
//		growth = 0.89696, growth2 = growth*growth;
		double pk0 = mpkval*b2+2./3.*growth*cpkval+1./5.*growth2*mpkval;
		double pk2 = 4./3.*growth*cpkval+4./7.*growth2*mpkval;
		wf << k << "\t" << pk0 << "\t" << pk2 << "\n";
        }
	wf.close();
}
