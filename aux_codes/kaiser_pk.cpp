#include <iostream>
#include <fstream>
#include <cmath>
#include <cassert>
#include <string>
#include "../generate_Poisson/spline.h"

using namespace std;
int count_num_line(const string &fname);

int main(int argc, char *argv[])
{
	string mpkfname;
	cout << "Please enter name of matter power spectrum file [1st: k in units of h/Mpc; 2nd: P(k) in units of Mpc^3/h^3]" << endl;
	cin >> mpkfname;
	string pkgmfname;
	cout << "Please enter name of cross power spectrum file [1st: k in units of h/Mpc; 2nd: P(k) in units of Mpc^3/h^3]" << endl;
	cin >> pkgmfname;
	double bias;
	cout << "Please enter the linear bias" << endl;
	cin >> bias;
	string ffname;         
	cout << "Please enter the logarithmic growth rate file [1st: k in units of h/Mpc; 2nd: f(k,z)]:" << endl;  
	cin >> ffname;
	string ofname;
	cout << "Please enter the output file name" << endl;
	cin >> ofname;

	int nmpk = count_num_line(mpkfname);
	int npkgm = count_num_line(pkgmfname);
	int ndataf = count_num_line(ffname);
	
	interpol mpk(mpkfname,nmpk,2,1,2,true,true);
	interpol cpk(pkgmfname,npkgm,2,1,2,true,true);
	interpol f(ffname, ndataf, 2, 1, 2, true, true);	

	double b2 = bias*bias;
	ofstream wf; wf.open(ofname.c_str());
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

int count_num_line(const string &fname){
	int i = 0;
	ifstream ifs(fname.c_str());
	if(ifs){
		string line;
		while(true){
			getline(ifs,line);
	        if(ifs.eof()) break;
	        i++;
		}
	}
	return i;
}
