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


int main(){

	string halo1fname;
	string halo2fname;
	cout << "Please enter the halo1 file name:" << endl;
	cin >> halo1fname;
	cout << "Please enter the halo2 file name:" << endl;
	cin >> halo2fname;
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

	ifstream galin1(halo1fname.c_str());
	ifstream galin2(halo2fname.c_str());
	double Lx,Ly,Lz;
	galin1.read((char*) &Lx,sizeof(double));
	galin1.read((char*) &Ly,sizeof(double));
	galin1.read((char*) &Lz,sizeof(double));
   cout <<"Lx, Ly, Lz From file 1:"<<endl;
	cout << Lx << Ly << Lz << endl;

	double Lx2,Ly2,Lz2;
	galin2.read((char*) &Lx2,sizeof(double));
	galin2.read((char*) &Ly2,sizeof(double));
	galin2.read((char*) &Lz2,sizeof(double));
   cout <<"Lx,Ly,Lz From file 2:"<<endl;
	cout << Lx2 << Ly2 << Lz2 << endl;

   if(fabs(Lx-Lx2)>1.e-6 || fabs(Ly-Ly2)>1.e-6 || fabs(Lz-Lz2)>1.e-6){
		cerr<< "box size not matching!"<<endl;
		exit(1);
	}

	double max_xy = (Lx>Ly)? Lx : Ly;
	double maxL = (Lz>max_xy)? Lz : max_xy;
	double min_xy = (Lx>Ly)? Ly : Lx;
	double minL = (Lz>min_xy)? min_xy : Lz;

	int n0 = int(round(nmesh * Lx / maxL));
	int n1 = int(round(nmesh * Ly / maxL));
	int n2 = int(round(nmesh * Lz / maxL));

	double rbin = Lx/double(n0);

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

	double *Nijk1 = new double [nrtot];
	double *Nijk2 = new double [nrtot];
	#pragma omp parallel for
	for (int i=0;i<nrtot;i++){
		Nijk1[i] = 0;
		Nijk2[i] = 0;
	}

	int Ngal1;
	galin1.read((char*) &Ngal1,sizeof(int));
	float *gal1 = new float [Ngal1*6];
	galin1.read((char*) gal1,sizeof(float)*Ngal1*6);
	galin1.close();

	int Ngal2;
	galin2.read((char*) &Ngal2,sizeof(int));
	float *gal2 = new float [Ngal2*6];
	galin2.read((char*) gal2,sizeof(float)*Ngal2*6);
	galin2.close();

	cout << "distributing particles..." << endl;

	float *galtmp = new float [3];	// modified by kayo	
	for (int i=0;i<Ngal1;i++){
		double xx = gal1[  6*i];
		double yy = gal1[1+6*i];
		double zz = gal1[2+6*i];
		double vx = gal1[3+6*i];
		double vy = gal1[4+6*i];
		double vz = gal1[5+6*i];

		double vel_los = (vx*los[0]+vy*los[1]+vz*los[2])/aH;

		galtmp[0] = xx+vel_los*los[0];
		galtmp[1] = yy+vel_los*los[1];
		galtmp[2] = zz+vel_los*los[2];

		if(galtmp[0] < 0.0) {cerr << "galtemp[0]: " << galtmp[0] << endl; galtmp[0]+=Lx;} // modified by kayo, not checked yet
		if(galtmp[1] < 0.0) {cerr << "galtemp[1]: " << galtmp[1] << endl; galtmp[1]+=Ly;} // modified by kayo, not checked yet
		if(galtmp[2] < 0.0) {cerr << "galtemp[2]: " << galtmp[2] << endl; galtmp[2]+=Lz;} // modified by kayo, not checked yet
						
		int ix = int(galtmp[0]/rbin);
		int iy = int(galtmp[1]/rbin);
		int iz = int(galtmp[2]/rbin);

		if(ix >= n0) {cerr << "ix: " << ix << endl; ix-=n0;} // modified by kayo
		if(iy >= n1) {cerr << "iy: " << iy << endl; iy-=n1;} // modified by kayo
		if(iz >= n2) {cerr << "iz: " << iz << endl; iz-=n2;} // modified by kayo

		Nijk1[iz+n2*(iy+n1*ix)] += 1.;
	}

	for (int i=0;i<Ngal2;i++){
		double xx = gal2[  6*i];
		double yy = gal2[1+6*i];
		double zz = gal2[2+6*i];
		double vx = gal2[3+6*i];
		double vy = gal2[4+6*i];
		double vz = gal2[5+6*i];

		double vel_los = (vx*los[0]+vy*los[1]+vz*los[2])/aH;

		galtmp[0] = xx+vel_los*los[0];
		galtmp[1] = yy+vel_los*los[1];
		galtmp[2] = zz+vel_los*los[2];

		if(galtmp[0] < 0.0) {cerr << "galtemp[0]: " << galtmp[0] << endl; galtmp[0]+=Lx;} // modified by kayo, not checked yet
		if(galtmp[1] < 0.0) {cerr << "galtemp[1]: " << galtmp[1] << endl; galtmp[1]+=Ly;} // modified by kayo, not checked yet
		if(galtmp[2] < 0.0) {cerr << "galtemp[2]: " << galtmp[2] << endl; galtmp[2]+=Lz;} // modified by kayo, not checked yet
						
		int ix = int(galtmp[0]/rbin);
		int iy = int(galtmp[1]/rbin);
		int iz = int(galtmp[2]/rbin);

		if(ix >= n0) {cerr << "ix: " << ix << endl; ix-=n0;} // modified by kayo
		if(iy >= n1) {cerr << "iy: " << iy << endl; iy-=n1;} // modified by kayo
		if(iz >= n2) {cerr << "iz: " << iz << endl; iz-=n2;} // modified by kayo

		Nijk2[iz+n2*(iy+n1*ix)] += 1.;
	}
	delete[] galtmp; // modified by kayo

	double Nbar1 = double(Ngal2)/double(nrtot);
	double Nbar2 = double(Ngal2)/double(nrtot);

	fftw_init_threads();
	fftw_plan_with_nthreads(nproc);

	fftw_complex *deltak1;
	deltak1 = (fftw_complex*) fftw_malloc(nktot * sizeof(fftw_complex));
	double *deltar1;
	deltar1 = (double *) deltak1;
	fftw_plan plan_delta1;
	plan_delta1 = fftw_plan_dft_r2c_3d(n0,n1,n2,deltar1,deltak1,FFTW_ESTIMATE);

	fftw_complex *deltak2;
	deltak2 = (fftw_complex*) fftw_malloc(nktot * sizeof(fftw_complex));
	double *deltar2;
	deltar2 = (double *) deltak2;
	fftw_plan plan_delta2;
	plan_delta2 = fftw_plan_dft_r2c_3d(n0,n1,n2,deltar2,deltak2,FFTW_ESTIMATE);

	#pragma omp parallel for
	for (int i0=0;i0<n0;i0++){
		for (int i1=0;i1<n1;i1++){
			for (int i2=0;i2<n2;i2++){
				int indx = i2+n2*(i1+n1*i0);
				deltar1[i2+cn22*(i1+n1*i0)] = Nijk1[indx]/Nbar1-1.;
				deltar2[i2+cn22*(i1+n1*i0)] = Nijk2[indx]/Nbar2-1.;
			}
		}
	}
	delete[] Nijk1;
	delete[] Nijk2;

	cout << "FFTing..." << endl;
	fftw_execute(plan_delta1);
	fftw_execute(plan_delta2);


	double *pkl = new double [nmax*lmax];
	for (int i=0;i<nmax*lmax;i++) pkl[i] = 0;
	cout << "estimating the raw power spectrum..." << endl;
	for (int i0=0;i0<n0;i0++){
		int nk0 = i0;
		int x0 = (n0-i0)%n0;
		if (i0>n0/2) nk0 = nk0-n0;
		long double k0 = double(nk0) * kF0;
//		double xcic = double(fabs(nk0))*PI/double(n0);
//		long double Wf0 = sinc(xcic);
//		long double Cf0 = 1.-2./3.*pow(sin(xcic),2);
		for (int i1=0;i1<n1;i1++){
			int nk1 = i1;
			int x1 = (n1-i1)%n1;
			if(i1>n1/2) nk1 = nk1-n1;
			long double k1 = double(nk1) * kF1;
//			xcic = double(fabs(nk1))*PI/double(n1);
//			long double Wf1 = sinc(xcic);
//			long double Cf1 = 1.-2./3.*pow(sin(xcic),2);
			for (int i2=0;i2<cn2;i2++){
				int nk2 = i2;
				long double k2 = double(nk2) * kF2;
//				xcic = double(fabs(nk2))*PI/double(n2);
//				long double Wf2 = sinc(xcic);
//				long double Cf2 = 1.-2./3.*pow(sin(xcic),2);

				bool real_dof = true;
				if(i2==0 || ((!(n2%2)) && i2==n2/2)) real_dof = (i1+n1*i0 <= x1 + n1*x0)? true : false;
				if(real_dof){
					int ijk = i2+cn2*(i1+n1*i0);
					double deltak1r = deltak1[ijk][0];
					double deltak1i = deltak1[ijk][1];
					double deltak2r = deltak2[ijk][0];
					double deltak2i = deltak2[ijk][1];

//					long double Pkadd = deltakr*deltakr+deltaki*deltaki;
					long double Pkadd = deltak1r*deltak2r+deltak1i*deltak2i;

					long double k = sqrt(k0*k0+k1*k1+k2*k2);
					int n = int(round(k/kbin));
					if (n<nmax){
						double mu = (k==0)? 0 : k2/k;
						for (int j0=0;j0<lmax;j0++){
							int l0 = j0*2;
							pkl[j0+lmax*n] += Pkadd*gsl_sf_legendre_Pl(l0,mu);;
						}
					}
				}
			}
		}
	}

	fftw_free(deltak1);
	fftw_free(deltak2);

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


//	double Pshot = volume/double(Ngal);
	for (int n=0;n<nmax;n++){
		double k = kbin*double(n);
		if (nmodes[n]>0){
			double pk0 = pkt[0+lmax*n]/volume*dr3*dr3;//-Pshot;
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

