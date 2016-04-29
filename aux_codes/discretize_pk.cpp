/* 
	This program calculate the discretized power spectrum from given power spectrum.
	Aniket Agrawal
*/
#include <iostream>
#include <string>
#include <fstream>
#include <cassert>
#include <cmath>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>

using namespace std;

int main(void)
{
	string pkfname;
	cout << "Please enter the theory power spectrum!" << endl;
	cin >> pkfname;
	int npkt;
	cout << "Enter size of above file!" << endl;
	cin >> npkt;
	double Rs, b, Rs2, b2;
	cout << "Enter smoothing scale (assuming Gaussian smoothing) and bias (linear)!" << endl;
	cin >> Rs >> b;
	string ofname;
	cout << "Please enter the output file to store the power spectrum of generated file!" <<endl;
	cin >> ofname;
	double lengthx,lengthy,lengthz;
        cout << "Please enter the length of the Fourier box in x, y, and z:" << endl;
        cin >> lengthx >> lengthy >> lengthz;
        int nmax;
        cout << "Please enter the maximum size of Fourier grid (1D) !" << endl;
        cin >> nmax;

	cout << "Setting up the arrays......."<<endl;
	double max_xy = lengthx > lengthy? lengthx: lengthy; 
	double maxlength = lengthz > max_xy? lengthz: max_xy;
// -------------------------------------------------------------------------
// Allocating the array
//
// n0,n1,n2 is determined so that 
// 1. max(n0,n1,n2) = nmax
// 2. Dbin is constant for all three directions so that the Nyquist frequency
//		is the same for all three directions.
// -------------------------------------------------------------------------
	int n0 = int(ceil(nmax * lengthx / maxlength));
	int n1 = int(ceil(nmax * lengthy / maxlength));
	int n2 = int(ceil(nmax * lengthz / maxlength));

	cout << "n0,n1,n2=" << n0 << "\t" << n1 << "\t" << n2 << endl;

	double Dbin = maxlength / double(nmax); 

	double Lx = Dbin * double(n0);
	double Ly = Dbin * double(n1);
	double Lz = Dbin * double(n2);

	cout << "size of Fourier grid is (n0,n1,n2)"<< endl;
	cout << "("<<n0<<","<<n1<<","<<n2<<")"<<endl;
	cout << "Fourier resolution is " << Dbin << "[Mpc/h]" << endl;

	int cn2 = n2/2+1;				// Size of last column in complex array
	long long Ntotal = n0 * n1;
	Ntotal = Ntotal * n2;			// overall number of real cells
	long long Nc= n0 * n1;
	Nc = Nc * cn2;					// overall number of complex cells (padding)
	double pi = 3.14159265358979;
	// double sN = sqrt(N);
	double kF0 = 2. * pi / (Lx); // size of x cell in k-space
	double kF1 = 2. * pi / (Ly); // size of y cell in k-space
	double kF2 = 2. * pi / (Lz); // size of z cell in k-space
	double minL= (Lx<Ly) ? Lx : Ly;
        minL = (Ly<Lz) ? Ly : Lz;
	Rs2=Rs*Rs; b2=b*b;

	std::cout << "kF0 = " << kF0 << " kF1 = " << kF1 << " kF2 = " << kF2 << "\n";

	double *kt = new double[npkt], *pkt = new double[npkt];
	std::ifstream rf(pkfname.c_str()); assert(rf.is_open());
    
    double pk2;
	for(int i=0; i<npkt; i++){
	rf >> kt[i] >> pkt[i] >> pk2;
    }
	rf.close();
	gsl_interp_accel *acc1=gsl_interp_accel_alloc();
        gsl_spline *spline1=gsl_spline_alloc(gsl_interp_cspline, npkt);
        gsl_spline_init(spline1, kt, pkt, npkt);

	double kbin = 2.*pi/minL;
	double kmax = kF2*double(cn2-1);
	int nkmax = int(ceil(kmax/kbin));
	double *pk1d = new double [nkmax];
	double *nmodes = new double [nkmax];
	for (int i=0;i<nkmax;i++){
	        pk1d[i] = 0;
	        nmodes[i] = 0;
	}
	
	for (int i0=0;i0<n0;i0++){
		int nk0 = i0;
	        int x0 = (n0-i0)%n0;
	        if (i0>n0/2) nk0 = nk0-n0;
	        long double k0 = double(nk0) * kF0;
	        for (int i1=0;i1<n1;i1++){
	                int nk1 = i1;
	                int x1 = (n1-i1)%n1;
	                if(i1>n1/2) nk1 = nk1-n1;
	                long double k1 = double(nk1) * kF1;
	                for (int i2=0;i2<cn2;i2++){
	                        int nk2 = i2;
	                        long double k2 = double(nk2) * kF2;
	                        bool real_dof = true;
	                        if(i2==0 || ((!(n2%2)) && i2==n2/2)) real_dof = (i1+n1*i0 <= x1 + n1*x0)? true : false;
	                        if(real_dof){
	                                long double ksq = k0*k0+k1*k1+k2*k2, k = sqrt(ksq);
        	                        int n = int(round(k/kbin));
	                                if (n<nkmax && n!=0){
	                                        pk1d[n] += gsl_spline_eval(spline1, k, acc1)*b2*exp(-ksq*Rs2);
	                                        nmodes[n] ++;
	                                }
		                }
		        }
        	}
	}

	ofstream pkout(ofname.c_str());
	if (!pkout.is_open()){
	        cerr << "Output file cannot be opened!" << endl;
	        return 1;
	}
	for (int i=0;i<nkmax;i++){
	        double k = kbin*double(i);
	        pk1d[i] = pk1d[i]/nmodes[i];     //no d3r*d3r required here
	        pkout << k << "\t" << pk1d[i] << "\t" <<  nmodes[i] << "\n";
	}
	pkout.close();
	gsl_interp_accel_free(acc1);
        gsl_spline_free(spline1);
	return 0;
}
