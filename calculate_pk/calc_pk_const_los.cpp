#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <fftw3.h>
#include <omp.h>
#include <gsl/gsl_sf_legendre.h>
#include "matrix_inversion.h"

#define PI 3.14159265358979

using namespace std;

void distribute_cic(float *gal);

int n0,n1,n2;
double rbin;
double *Nijk;

int main(){

	string halofname;
	cout << "Please enter the halo file name:" << endl;
	cin >> halofname;
	int nmesh;
	cout << "Please enter the mesh number for the longest dimension:" << endl;
	cin >> nmesh;
	double aH;
	cout << "Please enter the velocity to distance factor (aH) [km/s/(Mpc/h)]:" << endl;
	cin >> aH;
	double *los = new double [3];
	cout << "Please enter the line-of-sight vector:" << endl;
	cin >> los[0] >> los[1] >> los[2];
	double kbin;
	cout << "Please enter the kbin for power spectrum:" << endl;
	cout << "(0 for fundamental frequency)" << endl;
	cin >> kbin;
	double kmax;
	cout << "Please enter the kmax for power spectrum:" << endl;
	cout << "(0 for Nyquist frequency)" << endl;
	cin >> kmax;
	int lmax;
	cout << "Please enter the maximum lmax for Legendre polynomial reconstruction:" << endl;
	cout << "(lmax>=4 and only even order is counted)" << endl;
	cin >> lmax;
	string imulfname;
	cout << "Please enter the output inverse mu-leakage matrix file:" << endl;
	cin >> imulfname;
	string pkfname;
	cout << "Please enter the output pk file:" << endl;
	cin >> pkfname;

	if (lmax<4){
		cout << "input lmax<4, so will use lmax=4!" << endl;
		lmax = 4;
	}

	int nproc = 1;
	#pragma omp parallel
	{
		nproc = omp_get_max_threads();
	}
	cout << "number of threads are using: " << nproc << endl;

	double losr = 0;
	losr = sqrt(los[0]*los[0]+los[1]*los[1]+los[2]*los[2]);
	for (int i=0;i<3;i++){
		if (losr>0) los[i] = los[i]/losr;
		else los[i] = 0;
	}

	ifstream galin(halofname.c_str());
	double Lx,Ly,Lz;
	galin.read((char*) &Lx,sizeof(double));
	galin.read((char*) &Ly,sizeof(double));
	galin.read((char*) &Lz,sizeof(double));

	double Lbox[3] = {Lx,Ly,Lz};

	double max_xy = (Lx>Ly)? Lx : Ly;
	double maxL = (Lz>max_xy)? Lz : max_xy;
	double min_xy = (Lx>Ly)? Ly : Lx;
	double minL = (Lz>min_xy)? min_xy : Lz;

	n0 = int(round(nmesh * Lx / maxL));
	n1 = int(round(nmesh * Ly / maxL));
	n2 = int(round(nmesh * Lz / maxL));

	rbin = Lx/double(n0);

	int cn2 = int(n2/2)+1;
	int cn22 = cn2*2;
	int nrtot = n0*n1*n2;
	int nktot = n0*n1*cn2;

	double kF0 = 2.*PI/Lx;
	double kF1 = 2.*PI/Ly;
	double kF2 = 2.*PI/Lz;
	double volume = Lx*Ly*Lz;

	double dr3 = pow(rbin,3);

	if (kmax>kF2*double(cn2-1)){
		cerr << "kmax is greater than the Nyquist frequency, so use kmax=kNy!" << endl;
		kmax = kF2*double(cn2-1);
	}

	if (kbin==0) kbin = 2.*PI/minL;
	if (kmax==0) kmax = kF2*double(cn2-1);
	int nmax = int(ceil(kmax/kbin));
	lmax = lmax/2+1;
	int lmax2 = lmax*lmax;

	double *nmodes = new double [nmax];
	double *imul = new double [lmax2*nmax];
	ifstream imulin(imulfname.c_str());
	if (imulin.good()){
		imulin.seekg(0);
		double Lx1,Ly1,Lz1,kbin1,kmax1;
		int nmesh1,nmax1,lmax1;
		imulin.read((char*) &Lx1,sizeof(double));
		imulin.read((char*) &Ly1,sizeof(double));
		imulin.read((char*) &Lz1,sizeof(double));
		imulin.read((char*) &nmesh1,sizeof(int));
		imulin.read((char*) &kbin1,sizeof(double));
		imulin.read((char*) &kmax1,sizeof(double));
		imulin.read((char*) &nmax1,sizeof(int));
		imulin.read((char*) &lmax1,sizeof(int));
		if (Lx1!=Lx||Ly1!=Ly||Lz1!=Lz||nmesh1!=nmesh||kbin1!=kbin||kmax1!=kmax||nmax1!=nmax||lmax1!=lmax){
			cout << "parameters in the inverse mu-leakage file is inconsistent with the input!" << endl;
			goto recompute;
		}
		imulin.read((char*) nmodes,sizeof(double)*nmax);
		imulin.read((char*) imul,sizeof(double)*nmax*lmax2);
		imulin.close();
	} else{
		cout << "cannot find the inverse mu-leakage file, so recompute!" << endl;
		recompute:
		imulin.close();
		for (int i=0;i<nmax;i++) nmodes[i] = 0;
		double *mul = new double [lmax2*nmax];
		for (int i=0;i<lmax2*nmax;i++) mul[i] = 0;
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
						long double k = sqrt(k0*k0+k1*k1+k2*k2);
						int n = int(round(k/kbin));
						if (n<nmax){
							nmodes[n] ++;
							double mu = (k==0)? 0 : k2/k;
							for (int j0=0;j0<lmax;j0++){
								int l0 = j0*2;
								double Ll0 = gsl_sf_legendre_Pl(l0,mu);
								for (int j1=j0;j1<lmax;j1++){
									int l1 = j1*2;
									double Ll1 = gsl_sf_legendre_Pl(l1,mu);
									mul[j1+lmax*j0+lmax2*n] += Ll0*Ll1;
								}
							}
						}
					}
				}
			}
		}
		double *mulk = new double [lmax2];
		double *imulk = new double [lmax2];
//		#pragma omp parallel for
		for (int n=0;n<nmax;n++){
			if (nmodes[n]>0){
				for (int j0=0;j0<lmax;j0++){
					for (int j1=j0;j1<lmax;j1++){
						int indx = j1+lmax*j0;
						mulk[indx] = mul[indx+lmax2*n];
						if (j0!=j1) mulk[j0+lmax*j1] = mulk[indx];
					}
				}
				matrix_inversion_by_LU(mulk,imulk,lmax);
				for (int j=0;j<lmax2;j++) imul[j+lmax2*n] = imulk[j];
			}
		}
		ofstream imulout(imulfname.c_str());
		if (!imulout.is_open()){
			cerr << "Output file cannot be opened!" << endl;
			exit(1);
		}
		imulout.seekp(0);
		imulout.write((char*) &Lx,sizeof(double));
		imulout.write((char*) &Ly,sizeof(double));
		imulout.write((char*) &Lz,sizeof(double));
		imulout.write((char*) &nmesh,sizeof(int));
		imulout.write((char*) &kbin,sizeof(double));
		imulout.write((char*) &kmax,sizeof(double));
		imulout.write((char*) &nmax,sizeof(int));
		imulout.write((char*) &lmax,sizeof(int));
		imulout.write((char*) nmodes,sizeof(double)*nmax);
		imulout.write((char*) imul,sizeof(double)*nmax*lmax2);
		imulout.close();
	}

	Nijk = new double [nrtot];
	#pragma omp parallel for
	for (int i=0;i<nrtot;i++) Nijk[i] = 0;

	int Ngal;
	galin.read((char*) &Ngal,sizeof(int));
	float *gal = new float [Ngal*6];
	galin.read((char*) gal,sizeof(float)*Ngal*6);
	galin.close();

	cout << "distributing particles..." << endl;
	for (int i=0;i<Ngal;i++){
		double xx = gal[  6*i];
		double yy = gal[1+6*i];
		double zz = gal[2+6*i];
		double vx = gal[3+6*i];
		double vy = gal[4+6*i];
		double vz = gal[5+6*i];

		double vel_los = (vx*los[0]+vy*los[1]+vz*los[2])/aH;

		float *galtmp = new float [3];
		galtmp[0] = xx+vel_los*los[0];
		galtmp[1] = yy+vel_los*los[1];
		galtmp[2] = zz+vel_los*los[2];

		for (int j=0;j<3;j++){
			if (galtmp[j]<0) galtmp[j] = galtmp[j]+Lbox[j];
			if (galtmp[j]>Lbox[j]) galtmp[j] = galtmp[j]-Lbox[j];
		}

		distribute_cic(galtmp);
	}

	double Nbar = double(Ngal)/double(nrtot);

	fftw_init_threads();
	fftw_plan_with_nthreads(nproc);

	fftw_complex *deltak;
	deltak = (fftw_complex*) fftw_malloc(nktot * sizeof(fftw_complex));
	double *deltar;
	deltar = (double *) deltak;
	fftw_plan plan_delta;
	plan_delta = fftw_plan_dft_r2c_3d(n0,n1,n2,deltar,deltak,FFTW_ESTIMATE);

	#pragma omp parallel for
	for (int i0=0;i0<n0;i0++){
		for (int i1=0;i1<n1;i1++){
			for (int i2=0;i2<n2;i2++){
				int indx = i2+n2*(i1+n1*i0);
				deltar[i2+cn22*(i1+n1*i0)] = Nijk[indx]/Nbar-1.;
			}
		}
	}
	delete[] Nijk;

	cout << "FFTing..." << endl;
	fftw_execute(plan_delta);

	double *pkl = new double [nmax*lmax];
	for (int i=0;i<nmax*lmax;i++) pkl[i] = 0;
	cout << "estimating the raw power spectrum..." << endl;
	for (int i0=0;i0<n0;i0++){
		int nk0 = i0;
		int x0 = (n0-i0)%n0;
		if (i0>n0/2) nk0 = nk0-n0;
		long double k0 = double(nk0) * kF0;
		double xcic = double(fabs(nk0))*PI/double(n0);
//		long double Wf0 = sinc(xcic);
		long double Cf0 = 1.-2./3.*pow(sin(xcic),2);
		for (int i1=0;i1<n1;i1++){
			int nk1 = i1;
			int x1 = (n1-i1)%n1;
			if(i1>n1/2) nk1 = nk1-n1;
			long double k1 = double(nk1) * kF1;
			xcic = double(fabs(nk1))*PI/double(n1);
//			long double Wf1 = sinc(xcic);
			long double Cf1 = 1.-2./3.*pow(sin(xcic),2);
			for (int i2=0;i2<cn2;i2++){
				int nk2 = i2;
				long double k2 = double(nk2) * kF2;
				xcic = double(fabs(nk2))*PI/double(n2);
//				long double Wf2 = sinc(xcic);
				long double Cf2 = 1.-2./3.*pow(sin(xcic),2);

				bool real_dof = true;
				if(i2==0 || ((!(n2%2)) && i2==n2/2)) real_dof = (i1+n1*i0 <= x1 + n1*x0)? true : false;
				if(real_dof){
					int ijk = i2+cn2*(i1+n1*i0);
					double deltakr = deltak[ijk][0];
					double deltaki = deltak[ijk][1];
					long double Pkadd = deltakr*deltakr+deltaki*deltaki;
//					long double Pkadd_W = Pkadd/pow(Wf0*Wf1*Wf2,4);
					long double Pkadd_C = Pkadd/(Cf0*Cf1*Cf2);

					long double k = sqrt(k0*k0+k1*k1+k2*k2);
					int n = int(round(k/kbin));
					if (n<nmax){
						double mu = (k==0)? 0 : k2/k;
						for (int j0=0;j0<lmax;j0++){
							int l0 = j0*2;
							pkl[j0+lmax*n] += Pkadd_C*gsl_sf_legendre_Pl(l0,mu);
						}
					}
				}
			}
		}
	}

	fftw_free(deltak);

	cout << "estimating the true power spectrum (correcting the mu leakage)..." << endl;
	double *pkt = new double [nmax*lmax];
	double *imulk = new double [lmax2];
//	#pragma omp parallel for
	for (int n=0;n<nmax;n++){
		if (nmodes[n]>0){
			for (int j=0;j<lmax2;j++) imulk[j] = imul[j+lmax2*n];
			for (int j0=0;j0<lmax;j0++){
				pkt[j0+lmax*n] = 0;
				for (int j1=0;j1<lmax;j1++) pkt[j0+lmax*n] += imulk[j1+lmax*j0]*pkl[j1+lmax*n];
			}
		}
	}

	delete[] pkl;
	delete[] imulk;
	delete[] imul;

	ofstream pkout(pkfname.c_str());
	if (!pkout.is_open()){
		cerr << "Output file cannot be opened!" << endl;
		return 1;
	}

	double Pshot = volume/double(Ngal);
	for (int n=0;n<nmax;n++){
		double k = kbin*double(n);
		if (nmodes[n]>0){
			double pk0 = pkt[0+lmax*n]/volume*dr3*dr3-Pshot;
			double pk2 = pkt[1+lmax*n]/volume*dr3*dr3;
			double pk4 = pkt[2+lmax*n]/volume*dr3*dr3;
			pkout << k << "\t" << pk0 << "\t" << pk2 << "\t" << pk4 << "\t" << nmodes[n] << endl;
		}
	}
	pkout.close();

	delete[] pkt;
	delete[] nmodes;

	return 0;
}

void distribute_cic(float *gal){

	float xx = gal[0];
	float yy = gal[1];
	float zz = gal[2];

	int i0 = int(floor(xx/rbin-0.5));
	int i1 = int(floor(yy/rbin-0.5));
	int i2 = int(floor(zz/rbin-0.5));

	long double w0 = xx/rbin-0.5-i0;
	long double w1 = yy/rbin-0.5-i1;
	long double w2 = zz/rbin-0.5-i2;

	int i0p = i0+1;
	int i1p = i1+1;
	int i2p = i2+1;

	if(i0p == n0) i0p = -1;
	if(i1p == n1) i1p = -1;
	if(i2p == n2) i2p = -1;

//	periodic boundary condition
	if(i0==-1) i0 = n0-1;
	if(i1==-1) i1 = n1-1;
	if(i2==-1) i2 = n2-1;
	if(i0p==-1) i0p = 0;
	if(i1p==-1) i1p = 0;
	if(i2p==-1) i2p = 0;

	int ijk;
	long double weight;
	if(i0 > -1 && i1 > -1 && i2 > -1){
		ijk = i2+n2*(i1+n1*i0);
		weight = (1.-w0)*(1.-w1)*(1.-w2);
		Nijk[ijk] += weight;
	}
	if(i0p > -1 && i1 > -1 && i2 > -1){
		ijk = i2+n2*(i1+n1*i0p);
		weight =    (w0)*(1.-w1)*(1.-w2);
		Nijk[ijk] += weight;
	}
	if(i0 > -1 && i1p > -1 && i2 > -1){
		ijk = i2+n2*(i1p+n1*i0);
		weight = (1.-w0)*   (w1)*(1.-w2);
		Nijk[ijk] += weight;
	}
	if(i0 > -1 && i1 > -1 && i2p > -1){
		ijk = i2p+n2*(i1+n1*i0);
		weight = (1.-w0)*(1.-w1)*   (w2);
		Nijk[ijk] += weight;
	}
	if(i0 > -1 && i1p > -1 && i2p > -1){
		ijk = i2p+n2*(i1p+n1*i0);
		weight = (1.-w0)*   (w1)*   (w2);
		Nijk[ijk] += weight;
	}
	if(i0p > -1 && i1 > -1 && i2p > -1){
		ijk = i2p+n2*(i1+n1*i0p);
		weight =    (w0)*(1.-w1)*   (w2);
		Nijk[ijk] += weight;
	}
	if(i0p > -1 && i1p > -1 && i2 > -1){
		ijk = i2+n2*(i1p+n1*i0p);
		weight =    (w0)*   (w1)*(1.-w2);
		Nijk[ijk] += weight;
	}
	if(i0p > -1 && i1p > -1 && i2p > -1){
		ijk = i2p+n2*(i1p+n1*i0p);
		weight =    (w0)*   (w1)*   (w2);
		Nijk[ijk] += weight;
	}
}
