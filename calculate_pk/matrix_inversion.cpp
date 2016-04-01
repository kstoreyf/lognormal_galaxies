//	In this code, we apply the LU decomposition to calculate the inversion of a matrix

#include <iostream>
#include <cmath>

using namespace std;

void matrix_inversion_by_LU(double *Am,double *iAm,int dim){

	long double *Lm = new long double [dim*dim];
	long double *Um = new long double [dim*dim];

//	initialize matrix L and U
	for (int i=0;i<dim;i++){
		for (int j=0;j<dim;j++){
			int indx = j+i*dim;
			if (i==j){
				Lm[indx] = 1.;
			} else if (j>i){
				Lm[indx] = 0.;
			} else{
				Um[indx] = 0.;
			}
		}
	}

//	decompose A = L * U
	for (int j=0;j<dim;j++){
		for (int i=0;i<=j;i++){
			int indx = j+i*dim;
			Um[indx] = Am[indx];
			for (int k=0;k<=i-1;k++){
				Um[indx] -= Lm[k+i*dim]*Um[j+k*dim];
			}
		}
		for (int i=j+1;i<dim;i++){
			int indx = j+i*dim;
			Lm[indx] = Am[indx];
			for (int k=0;k<=j-1;k++){
				Lm[indx] -= Lm[k+i*dim]*Um[j+k*dim];
			}
			Lm[indx] = Lm[indx]/Um[j+j*dim];
		}
	}

	long double *iLm = new long double [dim*dim];
	for (int j=0;j<dim;j++){
		for (int i=0;i<=j;i++){
			int indx = j+i*dim;
			if (i==j){
				iLm[indx] = 1./Lm[indx];
			} else{
				iLm[indx] = 0.;
			}
		}
		for (int i=j+1;i<dim;i++){
			int indx = j+i*dim;
			iLm[indx] = 0.;
			for (int k=0;k<=i-1;k++){
				iLm[indx] -= Lm[k+i*dim]*iLm[j+k*dim];
			}
			iLm[indx] = iLm[indx]/Lm[i+i*dim];
		}
	}

	delete[] Lm;

//	compute inverse U
	long double detUm = 1.;
	for (int i=0;i<dim;i++) detUm = detUm*Um[i+i*dim];
	if (detUm==0){
		cerr << "Matrix cannot be inverted by LU decomposition!" << endl;
	}

	long double *iUm = new long double [dim*dim];
	for (int i=0;i<dim;i++){
		for (int j=0;j<=i;j++){
			int indx = j+i*dim;
			if (i==j){
				iUm[indx] = 1./Um[indx];
			} else{
				iUm[indx] = 0.;
			}
		}
		for (int j=i+1;j<dim;j++){
			int indx = j+i*dim;
			iUm[indx] = 0.;
			for (int k=0;k<=j-1;k++){
				iUm[indx] -= iUm[k+i*dim]*Um[j+k*dim];
			}
			iUm[indx] = iUm[indx]/Um[j+j*dim];
		}
	}

	delete[] Um;

//	iA = iU * iL
	for (int i=0;i<dim;i++){
		for (int j=0;j<dim;j++){
			int indx = j+i*dim;
			iAm[indx] = 0.;
			for (int k=0;k<dim;k++){
				iAm[indx] += iUm[k+i*dim]*iLm[j+k*dim];
			}
		}
	}

	delete[] iLm;
	delete[] iUm;
}

void matrix_inversion_by_Cholesky(double *Am,double *iAm,int dim){

	double *Lm = new double [dim*dim];

	for (int i=0;i<dim;i++) for (int j=0;j<dim;j++) Lm[j+dim*i] = 0;

	for (int i=0;i<dim;i++){
		for (int j=0;j<(i+1);j++){
			long double sum = 0;
			for (int k=0;k<j;k++) sum += Lm[i*dim+k]*Lm[j*dim+k];
			Lm[i*dim+j] = (i == j) ? sqrt(Am[i*dim+i]-sum) : (1.0/Lm[j*dim+j]*(Am[i*dim+j]-sum));
		}
	}

	long double detLm = 1.;
	for (int i=0;i<dim;i++) detLm = detLm*Lm[i+i*dim];
	if (detLm==0) cerr << "Matrix cannot be inverted by Cholesky decomposition because of Lm!" << endl;

	double *iLm = new double [dim*dim];
	for (int j=0;j<dim;j++){
		for (int i=0;i<=j;i++){
			int indx = j+i*dim;
			if (i==j){
				iLm[indx] = 1./Lm[indx];
			} else{
				iLm[indx] = 0.;
			}
		}
		for (int i=j+1;i<dim;i++){
			int indx = j+i*dim;
			iLm[indx] = 0.;
			for (int k=0;k<=i-1;k++) iLm[indx] -= Lm[k+i*dim]*iLm[j+k*dim];
			iLm[indx] = iLm[indx]/Lm[i+i*dim];
		}
	}

	delete[] Lm;

	double *iUm = new double [dim*dim];
	for (int i=0;i<dim;i++) for (int j=0;j<dim;j++) iUm[j+dim*i] = iLm[i+dim*j];

	for (int i=0;i<dim;i++){
		for (int j=0;j<dim;j++){
			int indx = j+dim*i;
			iAm[indx] = 0;
			for (int k=0;k<dim;k++) iAm[indx] += iUm[k+dim*i]*iLm[j+dim*k];
		}
	}

	delete[] iLm;
	delete[] iUm;
}

