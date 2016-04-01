#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include "spline.h"
#include <stdlib.h>
using namespace std;
/*----------------------------------------------------------------------*/
// Spline function to calculate the 2nd derivative of yy w.r.t. xx
/*----------------------------------------------------------------------*/
void interpol::spline(const double yp1,const double ypn){

	double u[ndata-1];
	if(yp1 > 0.99e30)
		yy2[0]=u[0]=0.0;
	else{
		yy2[0]=-0.5;
		u[0]=(3.0/(xx[1]-xx[0]))*((yy[1]-yy[0])/(xx[1]-xx[0])-yp1);
	}

	double sig,p;
	for(int i=1;i<ndata-1;i++){
		sig=(xx[i]-xx[i-1])/(xx[i+1]-xx[i-1]);
		p=sig*yy2[i-1]+2.0;
		yy2[i]=(sig-1.0)/p;
		u[i]=(yy[i+1]-yy[i])/(xx[i+1]-xx[i])-(yy[i]-yy[i-1])/(xx[i]-xx[i-1]);
		u[i]=(6.0*u[i]/(xx[i+1]-xx[i-1])-sig*u[i-1])/p;
	}
	double qn,un;
	if(ypn > 0.99e30)
		qn=un=0.0;
	else{
		qn=0.5;
		un=(3.0/(xx[ndata-1]-xx[ndata-2]))
			*(ypn-(yy[ndata-1]-yy[ndata-2])/(xx[ndata-1]-xx[ndata-2]));
	}
	yy2[ndata-1]=(un-qn*u[ndata-2])/(qn*yy2[ndata-2]+1.0);
	for(int k=ndata-2;k>=0;k--)
		yy2[k]=yy2[k]*yy2[k+1]+u[k];
}
//  (C) Copr. 1986-92 Numerical Recipes Software D041&0(9p#3.
//     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

/*----------------------------------------------------------------------*/
// Class constructor
/*----------------------------------------------------------------------*/
interpol::interpol(const string &afname,const int andata,int ancolumn,int anxx,int anyy,
		bool axlog=true,bool aylog=true){
//assigning the arguments to the private values 
	fname=afname;
	ncolumn=ancolumn;
	nxx=anxx-1;
	nyy=anyy-1;
	xlog=axlog;
	ylog=aylog;
// opening the file
	ifstream infile(fname.c_str());
	if(!infile.is_open()){
		cout<<"Input file "<<fname<<" could not be opened."<<endl;
		exit(1);
	}

	string temp;
// if andata is nonzero, then, it is "ndata"
	if(andata){
		ndata=andata;
	}
// if andata is zero, then, calculate it from the file
	else{
// get length of the file
		int fstart, fend, oneline, length;
// get a size of one line
		getline(infile,temp);
		oneline=temp.size()+1; // size of one line = size of character + "\n"
		infile.seekg(0,ios::end);
		fend=infile.tellg();
		infile.seekg(0,ios::beg);
		fstart=infile.tellg();
		length=(fend-fstart)/oneline;
		ndata=length;
		cout<<"ndata of the file is "<<ndata<<endl;
	}
// read from the file
	cout<<"reading from file "<<fname<<endl;

// array storing the data
	xx=new double[ndata];
	yy=new double[ndata];

//	double readdata;
	double readdata;
	for(int dindx=0;dindx<ndata;dindx++){
		getline(infile,temp);
	   istringstream iss(temp,istringstream::in);
		for(int cindx=0;cindx<ncolumn;cindx++){
			iss>>readdata;
			if(cindx==nxx){
				if(xlog) readdata=log(readdata);
				xx[dindx] = readdata;
			}
			if(cindx==nyy){
				if(ylog) readdata=log(readdata);
				yy[dindx] = readdata;
			}
		}
/*// to print out the input file reading
		cout.width(10);
		cout.precision(4);
		cout<<xx[dindx];
		cout.width(15);
		cout.precision(4);
		cout<<yy[dindx]<<endl;
*/
	}
	infile.close();
// array storing the 2nd derivative data
	yy2=new double[ndata];
	interpol::spline(1.e30,1.e30);
}
/*----------------------------------------------------------------------*/
// Copy constructor 
/*----------------------------------------------------------------------*/
interpol::interpol(const interpol &Other){
	ndata = Other.ndata;
	xlog  = Other.xlog;
	ylog  = Other.ylog;
	
	xx = new double[ndata];
	yy = new double[ndata];
	yy2 = new double[ndata];

	for(int indx=0; indx<ndata; indx++){
		xx[indx]  = Other.xx[indx];
		yy[indx]  = Other.yy[indx];
		yy2[indx] = Other.yy2[indx];
	}
}

/*----------------------------------------------------------------------*/
// Hunt function to calculate the maximum xx which does not exceed x
/*----------------------------------------------------------------------*/
void interpol::hunt(const double x,int &jlo){
	int jm,jhi,inc;
	bool ascnd=(xx[ndata-1]>=xx[0]);

	if(jlo <0 || jlo > ndata+1){
		jlo=-1;
		jhi=ndata;
	}
	else{
		inc=1;
		if( (x>=xx[jlo]) == ascnd ){
			if(jlo==ndata-1) return;
			jhi=jlo+1;
			while( (x>=xx[jhi]) == ascnd ){
				jlo=jhi;
				inc+=inc;
				jhi=jlo+inc;
				if(jhi>ndata-1){
					jhi=ndata;
					break;
				}
			}
		}
		else{
			if(jlo==0){
				jlo-=1;
				return;
			}
			jhi=jlo--;
			while( (x<xx[jlo]) == ascnd ){
				jhi=jlo;
				inc<<=1;
				if(inc>=jhi){
					jlo=-1;
					break;
				}
				else jlo=jhi-inc;
			}
		}
	}
	while(jhi-jlo !=1){
		jm=(jhi+jlo)>>1;
		if( (x>=xx[jm]) == ascnd )
			jlo=jm;
		else
			jhi=jm;
	}
	if(x==xx[ndata-1]) jlo=ndata-2;
	if(x==xx[0]) jlo=0;
}
//  (C) Copr. 1986-92 Numerical Recipes Software D041&0(9p#3.
//     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

/*----------------------------------------------------------------------*/
// function to calculate the interpolation
/*----------------------------------------------------------------------*/
double interpol::value(const double ax){
	double x;
	if(xlog){
		x=log(ax);
	}else{
		x=ax;
	}
	double a,b,h,y;
	int jlo=0;
	interpol::hunt(x,jlo);
// When x is outside of the table
	if(jlo==-1 || jlo== ndata-1){
		cout<<"argument outside of the bound "<<ax<<endl;
		exit(1);
	}else{
		h=xx[jlo+1]-xx[jlo];
		a=(xx[jlo+1]-x)/h;
		b=(x-xx[jlo])/h;
      y=a*yy[jlo]+b*yy[jlo+1]+((a*a*a-a)*yy2[jlo]+(b*b*b-b)*yy2[jlo+1])*(h*h)/6.;
	}
	if(ylog)	y=exp(y);
	return y;
}
/*----------------------------------------------------------------------*/
// function to calculate the derivative of the tabulated function
/*----------------------------------------------------------------------*/
double interpol::deriv(const double ax){
	double x;
	if(xlog){
		x=log(ax);
	}else{
		x=ax;
	}
	double a,b,h,dydx;
	int jlo=0;
	interpol::hunt(x,jlo);
// When x is outside of the table
	if(jlo==-1 || jlo== ndata-1){
		cout<<"argument outside of the bound"<<endl;
		exit(1);
	}else{
		h=xx[jlo+1]-xx[jlo];
		a=(xx[jlo+1]-x)/h;
		b=(x-xx[jlo])/h;
      dydx=(yy[jlo+1]-yy[jlo])/h+(-(3.*a*a-1.)*yy2[jlo]+(3.*b*b-1.)*yy2[jlo+1])*h/6.;
	}
	return dydx;
}

/*
//This is the test routine for spline.cpp
int main(){
	interpol pk("wmap_5yr_baosn_matterpower_z=0.dat",2,1,2,false,true);
	int nhunt;
	pk.hunt(1.e3,nhunt);
	cout<<"nhunt="<<nhunt<<endl;
	pk.get_xx(nhunt-1);
	pk.get_xx(nhunt);
	pk.get_xx(nhunt+1);

	cout<<"interpolation value is "<<pk.value(1.e3)<<endl;
	cout<<"has to be between "<<pk.get_yy(nhunt)<<" "<<pk.get_yy(nhunt+1)<<endl;

	ofstream ofile("test_spline.dat");
	double k=1.e-4;
	while(k<1.e3){
		ofile.width(10);
		ofile.precision(5);
		ofile<<k;
		ofile.width(15);
		ofile.precision(5);
		ofile<<pk.value(k)<<endl;
		k*=1.1e0;		
	}
	ofile.close();


	return 0;
}*/
