/* 
	This program calculate the poisson sample from given power spectrum.
	[1]
	A Gaussian random field is generated in Fourier space for given
	LogNormal power spectrum, then real space density for each grid is 
   calculated by the Fourier + Log transform.
	[2]
	Algorithm for poisson sampling for total number of galaxies of Ngals
   is explained in HETDEX_mock.pdf

	19 May 2012
	Donghui Jeong
	[v2] Generate velocities from the matter density, rather than the galaxy density divided by the linear bias, by Aniket Agrawal, March 2016
*/
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <fftw3.h>
#include "spline.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <omp.h>

using namespace std;

int main(void)
{
	string pkfname;
	cout << "Please enter the name of the Log-normal Pk file [1st: k in units of h/Mpc; 2nd: P(k) in units of Mpc^3/h^3]" << endl;
	cin >> pkfname;
	int ndatapk;
	cout << "Please enter the number of lines in the Pk file above" << endl;
	cin >> ndatapk;
	string mpkfname;	//modified by Aniket
	cout << "Please enter the name of the Log-normal matter Pk file [1st: k in units of h/Mpc; 2nd: P(k) in units of Mpc^3/h^3]" << endl;	//modified by Aniket
	cin >> mpkfname;	//modified by Aniket
	int ndatampk;		//modified by Aniket
	cout << "Please enter the number of lines in the matter Pk file above" << endl;
	cin >> ndatampk;
	double lengthx,lengthy,lengthz;
	cout << "Please enter the length of the Fourier box in x, y, and z:" << endl;
	cin >> lengthx >> lengthy >> lengthz;
	int nmax;
	cout << "Please enter the maximum size of Fourier grid (1D) !" << endl;
	cin >> nmax;
	int Ngalaxies;
	cout << "Please enter the total number of galaxies to sample!" << endl;
	cin >> Ngalaxies;
	double aH;
	cout << "Please enter the velocity to distance factor (aH) [km/s/(Mpc/h)]:" << endl;
	cin >> aH;
	string ffname;		//modified by Aniket
	cout << "Please enter the logarithmic growth rate file [1st: k in units of h/Mpc; 2nd: f(k,z)]:" << endl;	//modified by Aniket
	cin >> ffname;		//modified by Aniket
	int ndataf;		//modified by Aniket
	cout << "Please enter the number of lines in the logarithmic growth rate file above" << endl;	//modified by Aniket
	cin >> ndataf;		//modified by Aniket
	double bias;
	cout << "Please enter the linear bias:" << endl;
	cin >> bias;
	int seed;
	cout << "Please enter the random seed for Fourier space sampling!" << endl;
	cin >> seed;
	int Pseed;
	cout << "Please enter the random seed for position of galaxy!" << endl;
	cin >> Pseed;
	int useed;
	cout << "Please enter the random seed for selecting galaxies!" << endl;
	cin >> useed;
	string Poissonfname;
	cout << "Please enter the output file to store the poisson sample!" <<endl;
	cin >> Poissonfname;
	string Densityfname;	//modified by Aniket
	cout << "Please enter the output file to store density!" << endl;	//modified by Aniket
	cin >> Densityfname;	//modified by Aniket

// -------------------------------------------------------------------------
// gsl variables
// -------------------------------------------------------------------------
	// initialize the gsl random generator
	const gsl_rng_type *rt;
	// random number in k-space sampling
	gsl_rng *r;			// with 'seed'
	// random number for position
	gsl_rng *rpos;		// with 'Pseed'
	// random number for selecting galaxies
	gsl_rng *rselect;	// with 'useed'
	// random number for flux

	gsl_rng_env_setup();
	rt = gsl_rng_mt19937; // use mt19937 generator

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
	double volume = Lx*Ly*Lz;

	double d3r = Dbin*Dbin*Dbin;
	double d3k = 1./volume;

/*-------------------------------------------------------------------------*/
/* allocate arrays and plans for Fourier transformation                    */
/*-------------------------------------------------------------------------*/

	int nproc = 1;
	#pragma omp parallel
	{
		nproc = omp_get_max_threads();
	}
	fftw_init_threads();
	fftw_plan_with_nthreads(nproc);

	fftw_complex *deltak, *deltamk;// complex array in k-space	//modified by Aniket
	double *deltar, *deltamr; 	// storage arrays for positions and velocities	//modified by Aniket

	deltak = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nc);
	deltamk = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nc);	//modified by Aniket
	// aliasing to the double pointer deltar for in-place tranform
	deltar = (double *)deltak;
	deltamr = (double *)deltamk;	//modified by Aniket

	fftw_plan Pfftw_c2r, mPfftw_c2r;	//modified by Aniket
	Pfftw_c2r = fftw_plan_dft_c2r_3d(n0, n1, n2, deltak, deltar, FFTW_ESTIMATE);
	mPfftw_c2r = fftw_plan_dft_c2r_3d(n0, n1, n2, deltamk, deltamr, FFTW_ESTIMATE);	//modified by Aniket
//	fftw_plan Pfftw_r2c; // r2c is used only for the testing purpose.
//	Pfftw_r2c = fftw_plan_dft_r2c_3d(n0, n1, n2, deltar, deltak, FFTW_ESTIMATE);

//	declare velocity variables
	fftw_complex *vxk;
	fftw_complex *vyk;
	fftw_complex *vzk;
	vxk = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nc);
	vyk = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nc);
	vzk = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nc);

	double *vxr;
	double *vyr;
	double *vzr;
	vxr = (double *) vxk;
	vyr = (double *) vyk;
	vzr = (double *) vzk;

	fftw_plan vx_r2c;
	vx_r2c = fftw_plan_dft_r2c_3d(n0, n1, n2, vxr, vxk, FFTW_ESTIMATE);
	fftw_plan vx_c2r;
	vx_c2r = fftw_plan_dft_c2r_3d(n0, n1, n2, vxk, vxr, FFTW_ESTIMATE);
	fftw_plan vy_c2r;
	vy_c2r = fftw_plan_dft_c2r_3d(n0, n1, n2, vyk, vyr, FFTW_ESTIMATE);
	fftw_plan vz_c2r;
	vz_c2r = fftw_plan_dft_c2r_3d(n0, n1, n2, vzk, vzr, FFTW_ESTIMATE);

/*-----------------------------------------------------------------------------------------------------------------------------*/
/* opening and reading transfer function file, matter power spectrum and logarithimic growth rate                              */
/*-----------------------------------------------------------------------------------------------------------------------------*/
	interpol pk(pkfname,ndatapk,2,1,2,true,true);
	interpol mpk(mpkfname, ndatampk, 2,1,2,true,true);	//modified by Aniket
	interpol f(ffname, ndataf, 2,1,2,true,true);		//modified by Aniket

/*-------------------------------------------------------------------------*/
/* Generating the mock density field in k-space                            */
/*-------------------------------------------------------------------------*/
	cout<<"Generating mock density field in k-space"<<endl;

	r  = gsl_rng_alloc(rt);
	gsl_rng_set(r,seed);

	// First, generate the mock density field for 
	// deltak(0:n0-1,0:n1-1,1:cn2-2) when n2 is even
	// deltak(0:n0-1,0:n1-1,1:cn2-1) when n2 is odd
	// CAUTION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// 	Do NOT parallelize this loop. 
	// 	Order of random numbers will not be reproducable.
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// Note that we use the periodicity of the density field 
	// in Fourier space :: deltak[i0,i1,i2] = deltak[i0-n0,i1-n1,i2]
	int cn2max = (n2%2==1? cn2 : cn2-1);
	for(int i0 = 0; i0 < n0; i0++){
		double k0 = kF0 * double((i0 < n0/2) ? i0 : i0-n0);
		for(int i1 = 0; i1 < n1; i1++){
			double k1 = kF1 * double((i1 < n1/2) ? i1 : i1-n1);
			for(int i2 = 1; i2 < cn2max; i2++){
				double k2 = kF2 * double(i2);
				double kvalue = sqrt(k0*k0 + k1*k1 + k2*k2);

				double gauss = gsl_ran_gaussian(r,1.0);
				deltak[i2 + cn2*(i1 + n1*(i0))][0]
				=gauss*sqrt(0.5*pk.value(kvalue)*volume);
				deltamk[i2+cn2*(i1 + n1*(i0))][0]
				=gauss*sqrt(0.5*mpk.value(kvalue)*volume);	//modified by Aniket

				gauss = gsl_ran_gaussian(r,1.0);
				deltak[i2 + cn2*(i1 + n1*(i0))][1]
				=gauss*sqrt(0.5*pk.value(kvalue)*volume);
				deltamk[i2+cn2*(i1 + n1*(i0))][1]
				=gauss*sqrt(0.5*mpk.value(kvalue)*volume);	//modified by Aniket
			}
		}
	}

	// Generating deltak for boundary plane.
	// i2 = 0 && i2 = cn2-1 for even n2
	// i2 = 0               for odd  n2

	/* For those boundary plane, we have to impose the Hermitian condition
		explicitly. That is, we have to impose
		for n2=0 plane,

			deltak[i0,i1,0] = deltak*[-i0,-i1,0] = deltak*[n0-i0,n1-i1,0]

		and, in addition, for even n2,

			deltak[i0,i1,n2/2] = deltak*[-i0,-i1,-n2/2]
                            = deltak*[n0-i0,n1-i1,n2/2]

		<<Pseudo code>>
		do i0=0,n0-1
			x0 = mod(n0-i0,n0)
			do i1=0,n1-1
				x1 = mod(n1-i1,n0)
				if(x0==i0 && x1==i1)
					if(i0==0 && i1==0 && i2==0) then delta[0]=0
					else
					assign real number (Nyquest)
				if(i1+n1*i0 < x1+n1*x0) 
					then assign C to deltak[i0,i1] and C* to deltak[x0,x1]
				else
					Do not assign any value since delta[i0,i1] is already assigned.
			enddo
		enddo
	*/
	// For odd n2, fill only i2=0 plane
	// For even n2, fill i2=0, i2=n2/2 plane

	for(int i0 = 0; i0 < n0; i0++){
		int x0=(n0-i0)%n0;
		double k0 = kF0 * double((i0 < n0/2) ? i0 : i0-n0);
		for(int i1 = 0; i1 < n1; i1++){
			int x1=(n1-i1)%n1;
			double k1 = kF1 * double((i1 < n1/2) ? i1 : i1-n1);
			double kvalue0 = sqrt(k0*k0 + k1*k1);
			double kvalue1 = sqrt(k0*k0 + k1*k1 + pow(kF2*n2/2.,2.));
			if(i0==x0 && i1==x1){
				if(i0==0 && i1==0){
					deltak[0][0] = 0.;	// DC mode = 0 'cause <delta>=0
					deltak[0][1] = 0.;
					deltamk[0][0] = 0.;	//modified by Aniket
					deltamk[0][1] = 0.;	//modified by Aniket
					// For even n2, generate the real number thing deltak(0,0,n2/2)
					if(!(n2%2)){
						double gauss = gsl_ran_gaussian(r,1.0);
						deltak[n2/2+cn2*(i1 + n1*(i0))][0]
						=gauss*sqrt(pk.value(kvalue1)*volume);
						deltamk[n2/2+cn2*(i1 + n1*(i0))][0]
                                                =gauss*sqrt(mpk.value(kvalue1)*volume);		//modified by Aniket
						deltak[n2/2+cn2*(i1 + n1*(i0))][1] = 0.;
						deltamk[n2/2+cn2*(i1 + n1*(i0))][1] = 0.;	//modified by Aniket
					}
				}
				else{ // Nyquist frequency modes are real.
					double gauss = gsl_ran_gaussian(r,1.0);
					deltak[cn2*(i1 + n1*(i0))][0]
					=gauss*sqrt(pk.value(kvalue0)*volume);
					deltamk[cn2*(i1+n1*(i0))][0]
                                        =gauss*sqrt(mpk.value(kvalue0)*volume);			//modified by Aniket
					deltak[cn2*(i1 + n1*(i0))][1] = 0.;
					deltamk[cn2*(i1+n1*(i0))][1] = 0.;			//modified by Aniket
					// For even n2, repeat the same procedure for i2=n2/2 
					if(!(n2%2)){
						gauss = gsl_ran_gaussian(r,1.0);
						deltak[n2/2+cn2*(i1 + n1*(i0))][0]
						=gauss*sqrt(pk.value(kvalue1)*volume);
						deltamk[n2/2+cn2*(i1 + n1*(i0))][0]
                                                =gauss*sqrt(mpk.value(kvalue1)*volume);		//modified by Aniket
						deltak[n2/2+cn2*(i1 + n1*(i0))][1] = 0.;
						deltamk[n2/2+cn2*(i1 + n1*(i0))][1] = 0.;	//modified by Aniket
					}
				}
			}
			else{
				if(i1+n1*i0 < x1+n1*x0){
					// generate the gaussian random field for
					// deltak(i0,i1,0)
					double gauss = gsl_ran_gaussian(r,1.0);
					deltak[cn2*(i1 + n1*(i0))][0]
					=gauss*sqrt(0.5*pk.value(kvalue0)*volume);
					deltamk[cn2*(i1 + n1*(i0))][0]
					=gauss*sqrt(0.5*mpk.value(kvalue0)*volume);		//modified by Aniket
					gauss = gsl_ran_gaussian(r,1.0);
					deltak[cn2*(i1 + n1*(i0))][1]
					=gauss*sqrt(0.5*pk.value(kvalue0)*volume);
					deltamk[cn2*(i1 + n1*(i0))][1]
					=gauss*sqrt(0.5*mpk.value(kvalue0)*volume);		//modified by Aniket
					// assign the complex conjugate to
					// deltak(x0,x1,0)
					deltak[cn2*(x1 + n1*(x0))][0]=deltak[cn2*(i1 + n1*(i0))][0];
					deltamk[cn2*(x1 + n1*(x0))][0]=deltamk[cn2*(i1 + n1*(i0))][0];	//modified by Aniket
					deltak[cn2*(x1 + n1*(x0))][1]=-deltak[cn2*(i1 + n1*(i0))][1];
					deltamk[cn2*(x1 + n1*(x0))][1]=-deltamk[cn2*(i1 + n1*(i0))][1];	//modified by Aniket
					// For even n2, repeat the same procedure for i2=n2/2 
					if(!(n2%2)){
						gauss = gsl_ran_gaussian(r,1.0);
						deltak[n2/2+cn2*(i1 + n1*(i0))][0]
						=gauss*sqrt(0.5*pk.value(kvalue1)*volume);
						deltamk[n2/2+cn2*(i1 + n1*(i0))][0]
						=gauss*sqrt(0.5*mpk.value(kvalue1)*volume);	//modified by Aniket

						gauss = gsl_ran_gaussian(r,1.0);
						deltak[n2/2+cn2*(i1 + n1*(i0))][1]
						=gauss*sqrt(0.5*pk.value(kvalue1)*volume);
						deltamk[n2/2+cn2*(i1 + n1*(i0))][1]
						=gauss*sqrt(0.5*mpk.value(kvalue1)*volume);	//modified by Aniket

						deltak[n2/2+cn2*(x1 + n1*(x0))][0]
						=deltak[n2/2+cn2*(i1 + n1*(i0))][0];
						deltamk[n2/2+cn2*(x1 + n1*(x0))][0]
						=deltamk[n2/2+cn2*(i1 + n1*(i0))][0];		//modified by Aniket
						deltak[n2/2+cn2*(x1 + n1*(x0))][1]
						=-deltak[n2/2+cn2*(i1 + n1*(i0))][1];
						deltamk[n2/2+cn2*(x1 + n1*(x0))][1]
						=-deltamk[n2/2+cn2*(i1 + n1*(i0))][1];		//modified by Aniket
					}
				}
			}
		}
	}

	// free random generator
	gsl_rng_free(r);
	cout << "Finished generating mock density field." << endl;

/*-------------------------------------------------------------------------*/
/* Testing the mock density field                                          */
/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/
/* Store deltak into the temporary array                                   */
/*-------------------------------------------------------------------------*/
/*	fftw_complex* temp = new fftw_complex[Nc];

	int i0=0;
   int chunk = n0 / 10;
#pragma omp parallel private(i0)
	{
#pragma omp for schedule(dynamic, chunk)
		for(i0 = 0; i0 < n0; i0++){
			for(int i1 = 0; i1 < n1; i1++){
				for(int i2 = 0; i2 < cn2; i2++){
					temp[i2 + cn2*(i1 + n1*(i0))][0]
					=deltak[i2 + cn2*(i1 + n1*(i0))][0];
					temp[i2 + cn2*(i1 + n1*(i0))][1]
					=deltak[i2 + cn2*(i1 + n1*(i0))][1];
				}
			}
		}
	}
*/
/*-------------------------------------------------------------------------*/
/* Check if r2c(c2r(deltak)) = n0*n1*n2*deltak                             */
/* This test fails if deltak is not Hermitian.                             */
/*-------------------------------------------------------------------------*/
/*	fftw_execute(Pfftw_c2r);
	fftw_execute(Pfftw_r2c);

	for(i0 = 0; i0 < n0; i0++){
		for(int i1 = 0; i1 < n1; i1++){
			for(int i2=0; i2 < cn2; i2++){
//				cout<<deltak[i2 + cn2*(i1 + n1*(i0))][0]<<endl;
				double ratio1=deltak[i2 + cn2*(i1 + n1*(i0))][0]/temp[i2 + cn2*(i1 + n1*(i0))][0];
				double ratio2=deltak[i2 + cn2*(i1 + n1*(i0))][1]/temp[i2 + cn2*(i1 + n1*(i0))][1];
				if( (long(ratio1)-N) > .1 || (long(ratio2) - N) > .1 ){
					cout<<i0<<" "<<i1<<" "<<i2<<" "<<ratio1<<" "<<ratio2<<" "<<n0*n1*n2<<endl;
				}
			}
		}
	}
	delete[]	temp;
*/
/*-------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------*/
/* calculate density map and store it to the file                          */
/*-------------------------------------------------------------------------*/
	// c2r Fourier transform to generate the real space density map
	cout << "Doing FFT for density field." << endl;
	fftw_execute(Pfftw_c2r);
	fftw_execute(mPfftw_c2r);	//modified by Aniket
//	cout << "Finished." << endl;
//	fftw_free(deltak);

	// Physical size of the last column of deltar is cn2*2.
	// However, the last 2 (even n2) or 1 (odd n2) elements are junks,
	// which is the result of the in-place fftw padding.
	int cn22=cn2*2;

	// Normalize the Fourier output to get the density field
	// also calculate the density variance
	long double averageG = 0., averagemG = 0.;	//modified by Aniket
	long double sigmasqG = 0., sigmasqmG = 0.;	//modified by Aniket
	for(int i0 = 0; i0 < n0; i0++){
		for(int i1 = 0; i1 < n1; i1++){
			for(int i2 = 0; i2 < n2; i2++){
				int ijk = i2 + cn22*(i1 + n1*(i0));
				// normalizing deltar
				deltar[ijk] = deltar[ijk]*d3k;
				deltamr[ijk] = deltamr[ijk]*d3k;			//modified by Aniket
				averageG += (long double)deltar[ijk];
				sigmasqG += (long double)deltar[ijk]*deltar[ijk];
				averagemG += (long double)deltamr[ijk];			//modified by Aniket
				sigmasqmG += (long double)deltamr[ijk]*deltamr[ijk];	//modified by Aniket
			}
		}
	}
	cout<<"Average of Log-normal density field  :"<<averageG<<endl;
	cout<<"Variance of Log-normal density field :"<<sqrt(sigmasqG)<<endl;
	averageG = averageG/(long double)(Ntotal);
	sigmasqG = sigmasqG/(long double)(Ntotal) - averageG*averageG;
	cout<<"Average of Log-normal density field  :"<<averageG<<endl;
	cout<<"Variance of Log-normal density field :"<<sqrt(sigmasqG)<<endl;
	
	cout<<"Average of matter Log-normal density field :"<<averagemG<<endl;		//modified by Aniket
	cout<<"Variance of matter Log-normal density field :"<<sqrt(sigmasqmG)<<endl;	//modified by Aniket
	averagemG = averagemG/(long double)(Ntotal);					//modified by Aniket
	sigmasqmG = sigmasqmG/(long double)(Ntotal) - averagemG*averagemG;		//modified by Aniket
	cout<<"Average of matter Log-normal density field :"<<averagemG<<endl;		//modified by Aniket
	cout<<"Variance of matter Log-normal density field :"<<sqrt(sigmasqmG)<<endl;	//modified by Aniket
	

	long double GtoLN = exp(-0.5*sigmasqG), mGtoLN = exp(-0.5*sigmasqmG);		//modified by Aniket
	double deltar_max=0.;
	double deltar_min=100.;
	long double avg_dr = 0;
	long double var_dr = 0;

//	int *deltar_hist = new int [2010];
//	for (int i=0;i<2010;i++) deltar_hist[i] = 0;

	for(int i0 = 0; i0 < n0; i0++){
		for(int i1 = 0; i1 < n1; i1++){
			for(int i2 = 0; i2 < n2; i2++){
				int ijk = i2 + cn22*(i1 + n1*(i0));
				// Transform to the density field
				deltar[ijk] = GtoLN*exp(deltar[ijk]) - 1.;
				deltamr[ijk] = mGtoLN*exp(deltamr[ijk]) - 1.;		//modified by Aniket
				long double deltar_l = (long double)(deltar[ijk]);
				avg_dr += deltar_l;
				var_dr += deltar_l*deltar_l;
				// calculate the maximum deltar
				deltar_max = ((deltar_max > deltar[ijk])? deltar_max:deltar[ijk]);
				deltar_min = ((deltar_min < deltar[ijk])? deltar_min:deltar[ijk]);
				vxr[ijk] = deltamr[ijk];				//modified by Aniket
//				vxr[ijk] = deltar[ijk];
//				int indx = int((deltar[ijk]+1.)/0.1);
//				if (indx>2009) indx = 2009;
//				deltar_hist[indx] ++;
			}
		}
	}

//	ofstream deltarhout("deltar_histogram.dat");
//	for (int i=0;i<2010;i++){
//		double deltar_bin = -1.+double(i)*0.1;
//		deltarhout << deltar_bin << "\t" << deltar_hist[i] << endl;
//	}
//	deltarhout.close();

	cout << "density maximum = " << deltar_max << endl;
	cout << "density minimum = " << deltar_min << endl;
	avg_dr = avg_dr/(long double)(Ntotal);
	var_dr = var_dr/(long double)(Ntotal)-avg_dr*avg_dr;
	cout << "Average of density field: " << avg_dr << endl;
	cout << "Variance of density field: " << var_dr << endl;

	cout << "Doing FFT for the density field." << endl;
	fftw_execute(vx_r2c);

	cout << "Calculating velocity field in Fourier space..." << endl;
	#pragma omp parallel for
	for(int i0 = 0; i0 < n0; i0++){
		double k0 = kF0 * double((i0 < n0/2) ? i0 : i0-n0);
		for(int i1 = 0; i1 < n1; i1++){
			double k1 = kF1 * double((i1 < n1/2) ? i1 : i1-n1);
			for(int i2 = 0; i2 < cn2; i2++){
				double k2 = kF2 * double(i2);
				double k = sqrt(k0*k0 + k1*k1 + k2*k2);
				double const_fac = (k>0)? aH/k/k*d3r*(f.value(k)) : 0;	//modified by Aniket - smoothing doesn't seem to work

				int indx = i2+cn2*(i1+n1*i0);
				double deltakr = vxk[indx][0]*const_fac;
				double deltaki = vxk[indx][1]*const_fac;
				vxk[indx][0] = -k0*deltaki;
				vyk[indx][0] = -k1*deltaki;
				vzk[indx][0] = -k2*deltaki;
				vxk[indx][1] =  k0*deltakr;
				vyk[indx][1] =  k1*deltakr;
				vzk[indx][1] =  k2*deltakr;
			}
		}
	}

	cout << "Doing FFT for the vx field." << endl;
	fftw_execute(vx_c2r);
	cout << "Doing FFT for the vy field." << endl;
	fftw_execute(vy_c2r);
	cout << "Doing FFT for the vz field." << endl;
	fftw_execute(vz_c2r);

/*-------------------------------------------------------------------------*/
// Output density grid to get the density power spectrum
/*-------------------------------------------------------------------------*/	
//	string densityfname = "mock_density_grid.bin";
	ofstream densityout(Densityfname.c_str(),ios::binary);			//modified by Aniket
	if(!densityout.is_open()){
		cerr<<"Output file cannot be opened!"<<endl;
	}
	// output to file
	densityout.seekp(0);							//modified by Aniket
	densityout.write((char*) deltar, sizeof(double) * n0 * n1 * cn22);	//modified by Aniket
	densityout.close();							//modified by Aniket

	// renormalize the density field as
	// density = (1+density)/(1+density_max)
//	for(int i0 = 0; i0 < n0; i0++){
//		for(int i1 = 0; i1 < n1; i1++){
//			for(int i2 = 0; i2 < n2; i2++){
//				int ijk = i2 + cn22*(i1 + n1*(i0));
//				deltar[ijk] = (1.+deltar[ijk])/(1.+deltar_max);
//			}
//		}
//	}

// -------------------------------------------------------------------------
// Initialize the gsl random generator
// -------------------------------------------------------------------------
	cout << "Initializing random generater.." << endl;

	rpos  = gsl_rng_alloc(rt);
	gsl_rng_set(rpos,Pseed);
	rselect = gsl_rng_alloc(rt);
	gsl_rng_set(rselect,useed);

// -------------------------------------------------------------------------
// generating random samples for given total number of galaxies
// -------------------------------------------------------------------------	
	float xx,yy,zz;
	long double vx,vy,vz;
	float* galdata = new float[6*Ngalaxies*2];
	double urand;

	cout << "Generating Poisson particles..........."<<endl;

	long double ngalbar = (long double)(Ngalaxies)/(long double)(Ntotal);

	long double avg_vx = 0;
	long double var_vx = 0;
	double vx_min = 1e5;
	double vx_max = -1e5;
	long double avg_vy = 0;
	long double var_vy = 0;
	double vy_min = 1e5;
	double vy_max = -1e5;
	long double avg_vz = 0;
	long double var_vz = 0;
	double vz_min = 1e5;
	double vz_max = -1e5;

//	int *vz_hist = new int [1020*2];
//	for (int i=0;i<1020*2;i++) vz_hist[i] = 0;

	int nPoisson = 0;
	for(int i0 = 0; i0 < n0; i0++){
		for(int i1 = 0; i1 < n1; i1++){
			for(int i2 = 0; i2 < n2; i2++){
				int ijk = i2 + cn22*(i1 + n1*(i0));

				double xmin_thiscell = i0 * Dbin;
				double ymin_thiscell = i1 * Dbin;
				double zmin_thiscell = i2 * Dbin;

				double Nthiscell = ngalbar*(1.+deltar[ijk]);
				int Nthiscell_int = gsl_ran_poisson(rselect,Nthiscell);
//				int Nthiscell_int = (int) floor(Nthiscell);
//				Nthiscell = Nthiscell - Nthiscell_int;

				vx = vxr[ijk]*d3k;	//modified by Aniket
				vy = vyr[ijk]*d3k;	//modified by Aniket
				vz = vzr[ijk]*d3k;	//modified by Aniket

/*				avg_vx += vx;
				avg_vy += vy;
				avg_vz += vz;
				var_vx += vx*vx;
				var_vy += vy*vy;
				var_vz += vz*vz;
				if (vx<=vx_min) vx_min = vx;
				if (vy<=vy_min) vy_min = vy;
				if (vz<=vz_min) vz_min = vz;
				if (vx>=vx_max) vx_max = vx;
				if (vy>=vy_max) vy_max = vy;
				if (vz>=vz_max) vz_max = vz;*/

//				int indx = int(abs(vz)/50.);
//				vz_hist[indx*2] ++;

				// random number for selecting galaxy
//				urand = gsl_rng_uniform(rselect);
//				if(urand < Nthiscell) Nthiscell_int = Nthiscell_int+1;

				for(int ithiscell=0;ithiscell<Nthiscell_int;ithiscell++){
					urand = gsl_rng_uniform(rpos);
					xx = xmin_thiscell + urand * Dbin;
					urand = gsl_rng_uniform(rpos);
					yy = ymin_thiscell + urand * Dbin;
					urand = gsl_rng_uniform(rpos);
					zz = zmin_thiscell + urand * Dbin;
					// store position of this galaxy into array
					galdata[6*nPoisson  ] = xx;
					galdata[6*nPoisson+1] = yy;
					galdata[6*nPoisson+2] = zz;
					galdata[6*nPoisson+3] = vx;
					galdata[6*nPoisson+4] = vy;
					galdata[6*nPoisson+5] = vz;
//					vz_hist[indx*2+1] ++;
					nPoisson++;

					avg_vx += vx;
					avg_vy += vy;
					avg_vz += vz;
					var_vx += vx*vx;
					var_vy += vy*vy;
					var_vz += vz*vz;
					if (vx<=vx_min) vx_min = vx;
					if (vy<=vy_min) vy_min = vy;
					if (vz<=vz_min) vz_min = vz;
					if (vx>=vx_max) vx_max = vx;
					if (vy>=vy_max) vy_max = vy;
					if (vz>=vz_max) vz_max = vz;
				}
			}
		}
	}

	cout << "Total number of " << nPoisson << " galaxies are generated!" << endl;

	fftw_free(deltak);
	fftw_free(deltamk);	//modified by Aniket
	fftw_free(vxk);
	fftw_free(vyk);
	fftw_free(vzk);

/*	ostringstream cRs;
	cRs << Rs;
	string cfilter = (ifilter==0)? "tophat_" : "gaussian_";
	if (Rs==0) cfilter = "";
	string vzfname = "vz_histogram_" + cfilter + "Rs=" + cRs.str() + ".dat";

	ofstream vzhout(vzfname.c_str());
	for (int i=0;i<1020;i++){
		double vz_bin = double(i)*50.;
		vzhout << vz_bin << "\t" << vz_hist[i*2] << "\t" << vz_hist[i*2+1] << endl;
	}
	vzhout.close();*/

/*	avg_vx = avg_vx/(long double)(Ntotal);
	var_vx = var_vx/(long double)(Ntotal)-avg_vx*avg_vx;
	avg_vy = avg_vy/(long double)(Ntotal);
	var_vy = var_vy/(long double)(Ntotal)-avg_vy*avg_vy;
	avg_vz = avg_vz/(long double)(Ntotal);
	var_vz = var_vz/(long double)(Ntotal)-avg_vz*avg_vz;*/
	avg_vx = avg_vx/(long double)(nPoisson);
	var_vx = var_vx/(long double)(nPoisson)-avg_vx*avg_vx;
	avg_vy = avg_vy/(long double)(nPoisson);
	var_vy = var_vy/(long double)(nPoisson)-avg_vy*avg_vy;
	avg_vz = avg_vz/(long double)(nPoisson);
	var_vz = var_vz/(long double)(nPoisson)-avg_vz*avg_vz;

	cout << "min[vx] = " << vx_min << "\tmax[vx] = " << vx_max << endl;
	cout << "avg[vx] = " << avg_vx << "\tvar[vx] = " << var_vx << endl;
	cout << "min[vy] = " << vy_min << "\tmax[vy] = " << vy_max << endl;
	cout << "avg[vy] = " << avg_vy << "\tvar[vy] = " << var_vy << endl;
	cout << "min[vz] = " << vz_min << "\tmax[vz] = " << vz_max << endl;
	cout << "avg[vz] = " << avg_vz << "\tvar[vz] = " << var_vz << endl;

	ofstream fout(Poissonfname.c_str(),ios::binary);
	if(!fout.is_open()){
		cerr << "Output file cannot be opened!" << endl;
	}
	// output to file
	fout.seekp(0);
	fout.write((char*) &Lx, sizeof(double));
	fout.write((char*) &Ly, sizeof(double));
	fout.write((char*) &Lz, sizeof(double));
	fout.write((char*) &nPoisson, sizeof(int));
	fout.write((char*) galdata, sizeof(float)*6*nPoisson);
	fout.close();

	delete [] galdata;

	fftw_destroy_plan(Pfftw_c2r);
	fftw_destroy_plan(mPfftw_c2r);	

	fftw_destroy_plan(vx_r2c);
	fftw_destroy_plan(vx_c2r);
	fftw_destroy_plan(vy_c2r);
	fftw_destroy_plan(vz_c2r);

	return 0;
}
