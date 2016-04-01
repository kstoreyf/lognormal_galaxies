#ifndef SPLINE_H_
#define SPLINE_H_
#include<iostream>
#include<string>
#include<fstream>
#include<sstream>
#include<cmath>

using namespace std;

class interpol{
private:
	double *xx;		// pointer for x array
	double *yy;		// pointer for y array
	double *yy2;	// pointer for 2nd derivative of y array
	string fname;	// string variable for the file name
	int ndata;		// number of data points (the same as the length of the file)
	int ncolumn,nxx,nyy; // ncolumn 	= number of column
								// nxx		= the column for x
								// nyy		= the column for y
	bool xlog,ylog;		// whether interpolation is done in log interval or not
	void hunt(const double x,int &jlo);
	void spline(const double yp1,const double ypn);
public:
	interpol(const string &afname,const int andata,int ancolumn,int anxx,int anyy,
			bool axlog,bool aylog); // Constructor.
	interpol(const interpol &Other);
	~interpol(){
		delete [] xx;
		delete [] yy;
		delete [] yy2;
	}
	double value(const double x);
	double deriv(const double x);
	double get_xx(int n){
		return xx[n];
	}
	double get_yy(int n){
		return yy[n];
	}
};

#endif
