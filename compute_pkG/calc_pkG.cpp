/*

pk_G(k) = (4 pi) \int r^2 dr xi_G(r) j0(kr)

where j0(kr) = sin(kr) / (kr)

*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include "spline_array.h"
#include "integration_modified.cpp"

using namespace std;

#define PI 3.14159265358979

double sinc(double x);
double int_pkG_func(double r,int dim,double params[]);

int ndata;
double *r1d,*xi1d,*xi1d2;

void ReadParams(int argc, char* argv[]);
string ofile_prefix;
string xifname;
int ncolumn;
double bias;
double rmax;

int main(int argc, char* argv[]){

	ReadParams(argc, argv);

	double frac = bias*bias;
	//the column number (staring from 1) for r and xi(r) of the correlation function file
	int rc = 1; int xic = 2;

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

	string pkGfname = "./data/inputs/"+ofile_prefix+"_pkG_b" + cb.str() + ".dat";
	ofstream pkGout(pkGfname.c_str());
	pkGout.precision(10);
	for (int i=0;i<=450;i++){
		double k = 0.01*double(i)-3.5;
		k = pow(10.,k);
		params[0] = k;
		double pkG = romb.rombint(int_pkG_func,0,rmax,dim,params);
		if( pkG > 0.0){
			pkGout << k << "\t" << 4.*PI*pkG << endl;
		}
	}
	pkGout.close();

	return 0;
}

double int_pkG_func(double r,int dim,double params[]){

	double k = params[0];
	double xi;
	splint_arr(r1d,xi1d,xi1d2,ndata,r,&xi);
	double xiG = log(1.+xi);

	return r*r*xiG*sinc(k*r);
}

double sinc(double x){
	if (abs(x)<=1e-2){
		double x2 = x*x;
		return 1.-x2/6.+x2*x2/120.;
	}
	return sin(x)/x;
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
		cout << "Please enter the prefix of output file name:" << endl;
		cin >> ofile_prefix;
		cout << "Please enter the input correlation function file [1st column: r; 2nd column: xi(r)]:" << endl;
		cin >> xifname;
		cout << "Please enter the number of columns of the correlation function file:" << endl;
		cin >> ncolumn;
		cout << "Please enter the linear bias:" << endl;
		cin >> bias;
		cout << "Please enter rmax for integration:" << endl;
		cin >> rmax;
	}else if (argc == 6){
		ofile_prefix = argv[1];
		xifname = argv[2];
		if(check_int(argv[3],"ncolumn")) ncolumn = atoi(argv[3]);
		if(check_float(argv[4],"bias")) bias = atof(argv[4]);
		if(check_float(argv[5],"rmax")) rmax = atof(argv[5]);
	}else{
		cout << "number of arguments should be 0 or 5!!!" << endl;
		exit(1);
	}
}
