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
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#define PI 3.14159265358979

using namespace std;

void ReadParams(int argc, char* argv[]);

void distribute_ngp(float *gal, double *Nijk, int i_interlace);
void distribute_cic(float *gal, double *Nijk, int i_interlacei);
double sinc(double x);

int n0,n1,n2;
int nrtot;
double rbin;
double *Nijk, *Nijk2;
 
string halo1fname;
string halo2fname;
int nmesh;
double aH;
double *los = new double [3];
double kbin;
double kmax;
int lmax;
string imulfname;
string pkfname;
int calc_mode_pk;
int flag_mass_assign;
int num_interlacing;

int main(int argc, char* argv[]){

	ReadParams(argc, argv);

	if (lmax<4){
		cout << "input lmax<4, so will use lmax=4!" << endl;
		lmax = 4;
	}
	if (flag_mass_assign > 0 or num_interlacing > 1){
		cout << "currently calc cross only works for NGP and no interlacing mode!" << endl;
        exit(1);
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
	double Lx,Ly,Lz;
	galin1.read((char*) &Lx,sizeof(double));
	galin1.read((char*) &Ly,sizeof(double));
	galin1.read((char*) &Lz,sizeof(double));
    cout <<"Lx, Ly, Lz From file 1:"<<endl;
	cout << Lx << Ly << Lz << endl;

	double max_xy = (Lx>Ly)? Lx : Ly;
	double maxL = (Lz>max_xy)? Lz : max_xy;
	double min_xy = (Lx>Ly)? Ly : Lx;
	double minL = (Lz>min_xy)? min_xy : Lz;

	n0 = int(round(nmesh * Lx / maxL));
	n1 = int(round(nmesh * Ly / maxL));
	n2 = int(round(nmesh * Lz / maxL));

	rbin = Lx/double(n0);

	int cn2 = int(n2/2)+1;
	nrtot = n0*n1*n2;

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

	// set variables for cubic mode
	int n0_tmp = n0, n1_tmp = n1, n2_tmp = n2;
	double kF0_tmp = kF0, kF1_tmp = kF1, kF2_tmp = kF2;
	int n012_max = max(max(n0,n1),n2); 
	double kFmin = min(min(kF0,kF1),kF2);
	if(calc_mode_pk == 1){
		n0_tmp = n012_max; n1_tmp = n012_max; n2_tmp = n012_max;
		kF0_tmp = kFmin; kF1_tmp = kFmin; kF2_tmp = kFmin;
	}
	int cn2_tmp = int(n2_tmp/2)+1;
	//////////////////////////////

	double *nmodes = new double [nmax];
	double *imul = new double [lmax2*nmax];
	ifstream imulin(imulfname.c_str());
	if (imulin.good()){
		imulin.seekg(0);
		double Lx1,Ly1,Lz1,kbin1,kmax1;
		int nmesh1,nmax1,lmax1, calc_mode_pk1;
		imulin.read((char*) &Lx1,sizeof(double));
		imulin.read((char*) &Ly1,sizeof(double));
		imulin.read((char*) &Lz1,sizeof(double));
		imulin.read((char*) &nmesh1,sizeof(int));
		imulin.read((char*) &kbin1,sizeof(double));
		imulin.read((char*) &kmax1,sizeof(double));
		imulin.read((char*) &nmax1,sizeof(int));
		imulin.read((char*) &lmax1,sizeof(int));
		imulin.read((char*) &calc_mode_pk1,sizeof(int));
		if (Lx1!=Lx||Ly1!=Ly||Lz1!=Lz||nmesh1!=nmesh||kbin1!=kbin||kmax1!=kmax||nmax1!=nmax||lmax1!=lmax||calc_mode_pk1!=calc_mode_pk){
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

		for (int i0=0;i0<n0_tmp;i0++){
			int nk0 = i0;
			int x0 = (n0_tmp-i0)%n0_tmp;
			if (i0>n0_tmp/2) nk0 = nk0-n0_tmp;
			long double k0 = double(nk0) * kF0_tmp;
			for (int i1=0;i1<n1_tmp;i1++){
				int nk1 = i1;
				int x1 = (n1_tmp-i1)%n1_tmp;
				if(i1>n1_tmp/2) nk1 = nk1-n1_tmp;
				long double k1 = double(nk1) * kF1_tmp;
				for (int i2=0;i2<cn2_tmp;i2++){
					int nk2 = i2;
					long double k2 = double(nk2) * kF2_tmp;

					bool real_dof = true;
					if(i2==0 || ((!(n2_tmp%2)) && i2==n2_tmp/2)) real_dof = (i1+n1_tmp*i0 <= x1 + n1_tmp*x0)? true : false;
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
		imulout.write((char*) &calc_mode_pk,sizeof(int));
		imulout.write((char*) nmodes,sizeof(double)*nmax);
		imulout.write((char*) imul,sizeof(double)*nmax*lmax2);
		imulout.close();
	}

	double *Nijk1 = new double [nrtot];
	#pragma omp parallel for
	for (int i=0;i<nrtot;i++){
		Nijk1[i] = 0;
	}

	int Ngal1;
	galin1.read((char*) &Ngal1,sizeof(int));
	float *gal1 = new float [Ngal1*6];
	galin1.read((char*) gal1,sizeof(float)*Ngal1*6);
	galin1.close();

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

		for (int j=0; j<num_interlacing; j++){
			// for interlacing		
			galtmp[0] = galtmp[0]+rbin*(1./double(num_interlacing))*j;
			galtmp[1] = galtmp[1]+rbin*(1./double(num_interlacing))*j;
			galtmp[2] = galtmp[2]+rbin*(1./double(num_interlacing))*j;

			// periodic boundary conditions
			if(galtmp[0] < 0.0) galtmp[0]+=Lx;
			if(galtmp[1] < 0.0) galtmp[1]+=Ly;
			if(galtmp[2] < 0.0) galtmp[2]+=Lz;
			if(galtmp[0] > Lx) galtmp[0]-=Lx;
			if(galtmp[1] > Ly) galtmp[1]-=Ly;
			if(galtmp[2] > Lz) galtmp[2]-=Lz;

			if(flag_mass_assign == 0){
				distribute_ngp(galtmp, Nijk1, j);
			}else if(flag_mass_assign == 1){
				distribute_cic(galtmp, Nijk1, j);
			}
		}
	}
	delete[] galtmp; 
	double Nbar1 = double(Ngal1)/double(nrtot);

	fftw_init_threads();
	fftw_plan_with_nthreads(nproc);

	fftw_complex *deltak1[num_interlacing];
	double *deltar1[num_interlacing];
	fftw_plan plan_delta1[num_interlacing];
	for (int i=0; i<num_interlacing; i++){
		deltak1[i] = (fftw_complex*) fftw_malloc(n0_tmp*n1_tmp*cn2_tmp * sizeof(fftw_complex));
		deltar1[i] = (double *) deltak1[i];
		plan_delta1[i] = fftw_plan_dft_r2c_3d(n0_tmp,n1_tmp,n2_tmp,deltar1[i],deltak1[i],FFTW_ESTIMATE);
	}
	fftw_complex *deltak2[num_interlacing];
	double *deltar2[num_interlacing];
	fftw_plan plan_delta2[num_interlacing];
	for (int i=0; i<num_interlacing; i++){
		deltak2[i] = (fftw_complex*) fftw_malloc(n0_tmp*n1_tmp*cn2_tmp * sizeof(fftw_complex));
		deltar2[i] = (double *) deltak2[i];
		plan_delta2[i] = fftw_plan_dft_r2c_3d(n0_tmp,n1_tmp,n2_tmp,deltar2[i],deltak2[i],FFTW_ESTIMATE);
	}

	#pragma omp parallel for
	for (int i0=0;i0<n0;i0++){
		for (int i1=0;i1<n1;i1++){
			for (int i2=0;i2<n2;i2++){
				int indx = i2+n2*(i1+n1*i0);
				int indx2 = i2+cn2_tmp*2*(i1+n1_tmp*i0);
				for (int i=0; i<num_interlacing; i++){
					deltar1[i][indx2] = Nijk1[indx+i*nrtot]/Nbar1-1.;
//					deltar2[i][indx2] = Nijk1[indx+i*nrtot]/Nbar1-1.;
				}
			}
		}
	}
	delete[] Nijk1;

	// reading matter density file
	ifstream galin2(halo2fname.c_str());
	double *tmp = new double[n0*n1*cn2*2];
	galin2.read((char*) tmp, sizeof(double)*n0*n1*cn2*2);
	for (int i0=0;i0<n0;i0++){
		for (int i1=0;i1<n1;i1++){
			for (int i2=0;i2<n2;i2++){
				int indx2 = i2+cn2_tmp*2*(i1+n1_tmp*i0);
				for (int i=0; i<num_interlacing; i++){
//					deltar1[i][indx2] = tmp[indx2];
					deltar2[i][indx2] = tmp[indx2];
				}
			}
		}
	}
	galin2.close();

	cout << "FFTing..." << endl;
	for (int i=0; i<num_interlacing; i++){ 
		fftw_execute(plan_delta1[i]);
		fftw_execute(plan_delta2[i]);
	}

	double Pshot = 0.0; //test
	//double Pshot = volume/double(Ngal1); //test
	double *pkl = new double [nmax*lmax];
	for (int i=0;i<nmax*lmax;i++) pkl[i] = 0;

	cout << "estimating the raw power spectrum..." << endl;
	long double Pk_add = 0.0;
	for (int i0=0;i0<n0_tmp;i0++){
		int nk0 = i0;
		int x0 = (n0_tmp-i0)%n0_tmp;
		if (i0>n0_tmp/2) nk0 = nk0-n0_tmp;
		long double k0 = double(nk0) * kF0_tmp;
		double xcic = double(fabs(nk0))*PI/double(n0_tmp);
		long double Wf0 = sinc(xcic), If0 = cos(xcic);
		long double Hf0 = If0*(1.-sin(xcic)*sin(xcic)/6.);
		long double Cf0 = 1.-2./3.*pow(sin(xcic),2);

		for (int i1=0; i1<n1_tmp; i1++){
			int nk1 = i1;
			int x1 = (n1_tmp-i1)%n1_tmp;
			if(i1>n1_tmp/2) nk1 = nk1-n1_tmp;
			long double k1 = double(nk1) * kF1_tmp;
			xcic = double(fabs(nk1))*PI/double(n1_tmp);
			long double Wf1 = sinc(xcic), If1 = cos(xcic);
			long double Hf1 = If1*(1.-sin(xcic)*sin(xcic)/6.);
			long double Cf1 = 1.-2./3.*pow(sin(xcic),2);

			for (int i2=0;i2<cn2_tmp;i2++){
				int nk2 = i2;
				long double k2 = double(nk2) * kF2_tmp;
				xcic = double(fabs(nk2))*PI/double(n2_tmp);
				long double Wf2 = sinc(xcic), If2 = cos(xcic);
				long double Hf2 = If2*(1.-sin(xcic)*sin(xcic)/6.);
				long double Cf2 = 1.-2./3.*pow(sin(xcic),2);
				long double Cf012 = Cf0*Cf1*Cf2;
				long double W_fac = Wf0*Wf1*Wf2;

				bool real_dof = true;
				if(i2==0 || ((!(n2_tmp%2)) && i2==n2_tmp/2)) real_dof = (i1+n1_tmp*i0 <= x1 + n1_tmp*x0)? true : false;
				if(real_dof){
					int ijk = i2+cn2_tmp*(i1+n1_tmp*i0);

//					if(flag_mass_assign == 0){
//						//NGP
//						if(num_interlacing == 1){
							W_fac = 1.0;
	 						Pk_add = (deltak1[0][ijk][0]*deltak2[0][ijk][0]
                                     +deltak1[0][ijk][1]*deltak2[0][ijk][1])
                                     /pow(W_fac,2);
							Pk_add = Pk_add*dr3*dr3/volume-Pshot/pow(W_fac,2);
//						}else if(num_interlacing == 2){
//							double phase = (k0+k1+k2)*rbin*0.5;
//							double deltakr1 = deltak[0][ijk][0], deltakr2 = deltak[1][ijk][0];
//							double deltaki1 = deltak[0][ijk][1], deltaki2 = deltak[1][ijk][1];
//							double deltakr = 0.5*(deltakr1+cos(phase)*deltakr2-sin(phase)*deltaki2);
//							double deltaki = 0.5*(deltaki1+cos(phase)*deltaki2+sin(phase)*deltakr2);
//							Pk_add = (deltakr*deltakr+deltaki*deltaki)/pow(W_fac,2);
//							Pk_add = Pk_add*dr3*dr3/volume-Pshot/pow(W_fac,2.)*(0.5+0.5*If0*If1*If2);
//						}
//					}else if(flag_mass_assign == 1){
//						//CIC
//						if(num_interlacing == 1){
//							Pk_add = (pow(deltak[0][ijk][0],2)+pow(deltak[0][ijk][1],2))/pow(W_fac,4);
//							Pk_add = Pk_add*dr3*dr3/volume-Pshot*Cf012/pow(W_fac,4);
//						}else if(num_interlacing == 2){
//							double phase = (k0+k1+k2)*rbin*0.5;
//							double deltakr1 = deltak[0][ijk][0], deltakr2 = deltak[1][ijk][0];
//							double deltaki1 = deltak[0][ijk][1], deltaki2 = deltak[1][ijk][1];
//							double deltakr = 0.5*(deltakr1+cos(phase)*deltakr2-sin(phase)*deltaki2);
//							double deltaki = 0.5*(deltaki1+cos(phase)*deltaki2+sin(phase)*deltakr2);	
//							Pk_add = (deltakr*deltakr+deltaki*deltaki)/pow(W_fac,4);
//							Pk_add = Pk_add*dr3*dr3/volume-Pshot/pow(W_fac,4)*0.5*(Cf012+Hf0*Hf1*Hf2);
//						}
//					}

					long double k = sqrt(k0*k0+k1*k1+k2*k2);
					int n = int(round(k/kbin));
					if (n<nmax){
						double mu = (k==0)? 0 : k2/k;
						for (int j0=0;j0<lmax;j0++){
							int l0 = j0*2;
							pkl[j0+lmax*n] += Pk_add*gsl_sf_legendre_Pl(l0,mu);
						}
					}
				}
			}
		}
	}

	for (int i=0; i<num_interlacing; i++){
		fftw_free(deltak1[i]);
		fftw_free(deltak2[i]);
	}

	cout << "estimating the true power spectrum (correcting the mu leakage)..." << endl;
	double *pkt = new double [nmax*lmax];
	double *imulk = new double [lmax2];
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

	for (int n=0;n<nmax;n++){
		double k = kbin*double(n);
		if (nmodes[n]>0){
			double pk0 = pkt[0+lmax*n];
			double pk2 = pkt[1+lmax*n];
			double pk4 = pkt[2+lmax*n];
			pkout << k << "\t" << pk0 << "\t" << pk2 << "\t" << pk4 << "\t" << nmodes[n] << endl;
		}
	}
	pkout.close();

	delete[] pkt;
	delete[] nmodes;

	return 0;
}

bool is_int(const string &str){
    return str.find_first_not_of("0123456789") == string::npos;
}

bool is_float(const string &str){
	signed int dec_point = str.find_first_of(".");
	if ((dec_point >= 0) and (str.find(".",dec_point+1) == string::npos)){
		   return str.find_first_not_of("0123456789.") == string::npos;
	}else{
		return 0;
	}
}

bool check_int(const string &str, const string &param_name){
	if(is_int(str)){
		return 1;
	}else{
		cout << param_name <<" should be integer!"<<endl;
		exit(1);
	}
}

bool check_float(const string &str, const string &param_name){
	if(is_float(str)){
		return 1;
	}else{
		cout << param_name <<" should be float!"<<endl;
		exit(1);
	}
}

void ReadParams(int argc, char* argv[]){
	if (argc == 1){
    	cout << "Please enter the halo1 file name:" << endl;
    	cin >> halo1fname;
    	cout << "Please enter the halo2 file name:" << endl;
    	cin >> halo2fname;
		cout << "Please enter the mesh number for the longest dimension:" << endl;
    	cin >> nmesh;
    	cout << "Please enter the velocity to distance factor (aH) [km/s/(Mpc/h)]:" << endl;
    	cin >> aH;
    	cout << "Please enter the line-of-sight vector:" << endl;
    	cin >> los[0] >> los[1] >> los[2];
    	cout << "Please enter the kbin for power spectrum:" << endl;
    	cout << "(0 for fundamental frequency)" << endl;
    	cin >> kbin;
    	cout << "Please enter the kmax for power spectrum:" << endl;
    	cout << "(0 for Nyquist frequency)" << endl;
    	cin >> kmax;
    	cout << "Please enter the maximum lmax for Legendre polynomial reconstruction:" << endl;
    	cout << "(lmax>=4 and only even order is counted)" << endl;
    	cin >> lmax;
    	cout << "Please enter the output inverse mu-leakage matrix file:" << endl;
    	cin >> imulfname;
    	cout << "Please enter the output pk file:" << endl;
    	cin >> pkfname;
		cout << "Please enter the calculation mode:" << endl;
		cout << "(0 for rectangular box, 1 for cubic box)" << endl;
		cin >> calc_mode_pk;
		cout << "Mass assignment scheme?(0: NGC, 1:CIC):" << endl;
		cin >> flag_mass_assign;
		cout << "number of interlacing?(1 for no interlacing. num_interlacing > 2 doesn't supported yet):" << endl;
		cin >> num_interlacing;
	}else if (argc == 16){
		halo1fname = argv[1];
		halo2fname = argv[2];
		if (check_int(argv[3],"nmesh")) nmesh = atoi(argv[3]);
		if (check_float(argv[4],"aH")) aH = atof(argv[4]);
		if (check_float(argv[5],"los[0]")) los[0] = atof(argv[5]);
		if (check_float(argv[6],"los[1]")) los[1] = atof(argv[6]);
		if (check_float(argv[7],"los[2]")) los[2] = atof(argv[7]);
		if (check_float(argv[8],"kbin")) kbin = atof(argv[8]);
		if (check_float(argv[9],"kmax")) kmax = atof(argv[9]);
		if (check_int(argv[10],"lmax")) lmax = atoi(argv[10]);
		imulfname = argv[11];
		pkfname = argv[12];
		if (check_int(argv[13],"calc_mode_pk")) calc_mode_pk = atoi(argv[13]);
		if (check_int(argv[14],"flag_mass_assign")) flag_mass_assign = atoi(argv[14]);
		if (check_int(argv[15],"num_interlacing")) num_interlacing = atoi(argv[15]);
	}else{
		cout << "number of arguments should be 0 or 15!!!" << endl;
		exit(1);
	}

	if(calc_mode_pk != 0 and calc_mode_pk != 1){
		cout << "calc_mode_pk must be 0 or 1!!!" << endl;
		exit(1);
	}
	if(flag_mass_assign != 0 and flag_mass_assign != 1){
		cout << "flag_mass_assign must be 0 or 1!!!" << endl;
		exit(1);
	}
	if(num_interlacing != 1 and num_interlacing != 2){
		cout << "num_interlacing must be 1 or 2!!!" << endl;
		exit(1);
	}

}

double sinc(double x){
	double x2 = x*x, x4 = x2*x2;
	if(abs(x)<=1e-2)
		return 1.-x2/6.+x4/120.;
	else
		return sin(x)/x;
}

void distribute_ngp(float *gal, double *Nijk, int i_interlace){

		int ix = int(gal[0]/rbin);
		int iy = int(gal[1]/rbin);
		int iz = int(gal[2]/rbin);

		// periodic boundary conditions
		if(ix >= n0) ix-=n0;
		if(iy >= n1) iy-=n1;
		if(iz >= n2) iz-=n2;

		int index = iz+n2*(iy+n1*ix)+nrtot*i_interlace;
		Nijk[index] += 1.;
}

void distribute_cic(float *gal, double *Nijk, int i_interlace){

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
		ijk = i2+n2*(i1+n1*i0)+nrtot*i_interlace;
		weight = (1.-w0)*(1.-w1)*(1.-w2);
		Nijk[ijk] += weight;
	}
	if(i0p > -1 && i1 > -1 && i2 > -1){
		ijk = i2+n2*(i1+n1*i0p)+nrtot*i_interlace;
		weight =    (w0)*(1.-w1)*(1.-w2);
		Nijk[ijk] += weight;
	}
	if(i0 > -1 && i1p > -1 && i2 > -1){
		ijk = i2+n2*(i1p+n1*i0)+nrtot*i_interlace;
		weight = (1.-w0)*   (w1)*(1.-w2);
		Nijk[ijk] += weight;
	}
	if(i0 > -1 && i1 > -1 && i2p > -1){
		ijk = i2p+n2*(i1+n1*i0)+nrtot*i_interlace;
		weight = (1.-w0)*(1.-w1)*   (w2);
		Nijk[ijk] += weight;
	}
	if(i0 > -1 && i1p > -1 && i2p > -1){
		ijk = i2p+n2*(i1p+n1*i0)+nrtot*i_interlace;
		weight = (1.-w0)*   (w1)*   (w2);
		Nijk[ijk] += weight;
	}
	if(i0p > -1 && i1 > -1 && i2p > -1){
		ijk = i2p+n2*(i1+n1*i0p)+nrtot*i_interlace;
		weight =    (w0)*(1.-w1)*   (w2);
		Nijk[ijk] += weight;
	}
	if(i0p > -1 && i1p > -1 && i2 > -1){
		ijk = i2+n2*(i1p+n1*i0p)+nrtot*i_interlace;
		weight =    (w0)*   (w1)*(1.-w2);
		Nijk[ijk] += weight;
	}
	if(i0p > -1 && i1p > -1 && i2p > -1){
		ijk = i2p+n2*(i1p+n1*i0p)+nrtot*i_interlace;
		weight =    (w0)*   (w1)*   (w2);
		Nijk[ijk] += weight;
	}
}
