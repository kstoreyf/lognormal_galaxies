#include <iostream>
#include <cmath>
#include "integration_modified.h"
#include <stdlib.h>

using namespace std;

/* rombint returns the integral from a to b of func using Romberg integration.
	The method converges provided that f(x) is continuous in (a,b).
	tol indicates the desired relative accuracy in the integral. 
*/
template<class T>
double integration::rombint(T &func, const double a, const double b,int dim,double params[]){
	const int MAXITER = 40;
	const int MAXJ = 5;
	double g[MAXJ+1];

	int nint,i,j,k,jmax;
	double h,gmax,g0,fourj,g1,error;

	h=0.5*(b-a);
	gmax=h*(func(a,dim,params)+func(b,dim,params));
	g[0]=gmax;
	nint=1;
	error=1.e20;
	g0=0.;

	for(i=0; !((i > MAXITER) || ((i > 5) && (fabs(error) < tol))); i++){
	/*     Calculate next trapezoidal rule approximation to integral. */
		g0=0.;
		for (k=1; k<=nint; k++){
			g0+=func(a+(k+k-1)*h,dim,params);
		}
   
		g0=0.5*g[0]+h*g0;
		h=0.5*h;
		nint*=2;
		jmax=((i<MAXJ)? i : MAXJ);
		fourj=1.;
	/*     Use Richardson extrapolation. */
		for (j=1; j<=jmax; j++){
			fourj*=4.0;
			g1=g0+(g0-g[j-1])/(fourj-1.0);
			g[j-1]=g0;
			g0=g1;
		}
		if (fabs(g0) > tol){
			error=1.0-gmax/g0;
		}
      else{ 
			error=gmax;
		}
		gmax=g0;
		g[jmax]=g0;
	}
	if ((i > MAXITER) && (fabs(error) > tol)){
		cout << "rombint failed to converge; integral=" << g0 << ", error=" << error << endl;
	}
	return g0;
}

/* rk5int returns the integral from a to b of func using 
   5th order Runge-Kutta integration.
	tol indicates the desired relative accuracy in the integral. 
*/
template<class T>
void integration::rkck
(const double y,const double dydx,const double x,const double h,
double &yout,double &yerr,T &func,int dim,double params[]){
	double ak2,ak3,ak4,ak5,ak6;
   static const double A2=0.2,A3=0.3,A4=0.6,A5=1.0,
      A6=0.875,B21=0.2,B31=3.0/40.0,B32=9.0/40.0,
      B41=0.3,B42=-0.9,B43=1.2,B51=-11.0/54.0,
      B52=2.5,B53=-70.0/27.0,B54=35.0/27.0,
      B61=1631.0/55296.0,B62=175.0/512.0,
      B63=575.0/13824.0,B64=44275.0/110592.0,
      B65=253.0/4096.0,C1=37.0/378.0,
      C3=250.0/621.0,C4=125.0/594.0,
      C6=512.0/1771.0,DC1=C1-2825.0/27648.0,
      DC3=C3-18575.0/48384.0,DC4=C4-13525.0/55296.0,
      DC5=-277.0/14336.0,DC6=C6-0.25;

    ak2=func(x+A2*h,dim,params);
    ak3=func(x+A3*h,dim,params);
    ak4=func(x+A4*h,dim,params);
    ak5=func(x+A5*h,dim,params);
    ak6=func(x+A6*h,dim,params);

    yout=y+h*(C1*dydx+C3*ak3+C4*ak4+C6*ak6);
    yerr=h*(DC1*dydx+DC3*ak3+DC4*ak4+DC5*ak5+DC6*ak6);
}

template<class T>
void integration::rkqs
(double &y,const double dydx,double &x,const double htry,const double yscal,
double &hdid,double &hnext,T &func,int dim,double params[]){
	double errmax,h,htemp,xnew;
	double yerr,ytemp;
	static const double SAFETY=0.9, PGROW=-0.2, PSHRNK=-0.25, ERRCON=1.89e-4;

    h=htry;

	while(1){
		rkck(y,dydx,x,h,ytemp,yerr,func,dim,params);
		errmax=fabs(yerr/yscal)/tol;
		if (errmax <= 1.0) break;
		htemp=SAFETY*h*pow(errmax,PSHRNK);
		h=((h>=0.0)? 
			((htemp>0.1*h)? htemp: 0.1*h): ((htemp>0.1*h)? 0.1*h: htemp));
		xnew=x+h;
		if (xnew == x){
			cout<<"stepsize underflow in rkqs!"<<endl;
		}
	}
	if (errmax > ERRCON){
		hnext=SAFETY*h*pow(errmax,PGROW);
	}
	else{
		hnext=5.0*h;
	}
	hdid=h;
	x=x+h;
	y=ytemp;
}
template<class T>
double integration::rk5int(T &func, const double a, const double b,int dim,double params[]){
	const double h1=0.01;
	const double TINY=1.e-20;
	const double hmin=0.;
	const int MAXSTEP=100000;
	int nstp;
	double result;
	double h,hdid,hnext,x;
	double dydx,y,yscal;

	x=a;
	h=h1*((b > a) ? 1. : -1.);
	result=0.;
	y=result;

	for(nstp=1; nstp<=MAXSTEP; nstp++){
		dydx=func(x,dim,params);
		yscal=fabs(y)+fabs(h*dydx)+TINY;
		if((x+h-b)*(x+h-a) > 0.) h=b-x;
		rkqs(y,dydx,x,h,yscal,hdid,hnext,func,dim,params);
		if((x-b)*(b-a)>=0.){
			result=y;
			return result;
		}
		if(fabs(hnext) < hmin){
			cout<<"stepsize is smaller than minimum in rk_integral"<<endl;
		}
		h=hnext;
	}
	cout<<"too many steps in rk5int"<<endl;
	exit(1);
}

/*
//What follows are the testing code for integration class.

class printpower{
private:
	struct npower{
		double n;
		double operator()(const double xx){
//			return pow(xx,n);
			return sin(xx*n);
		}
	};
	npower power;
public:
	printpower(const double an){
		power.n = an;
	}
	void print_int(const double a, const double b){
		integration romb(1.e-10);
		double result=romb.rombint(power,a,b);
		cout<<"result is from romb "<<result<<endl;	
		integration rk5(1.e-10);
		result=rk5.rk5int(power,a,b);
		cout<<"result is from rk5 "<<result<<endl;	
	}
};

double test(double x){
	return sin(x*3.);
}
int main(){
//	npower square(2.);
	printpower square(3.);
	square.print_int(3.,7.);

	integration romb(1.e-10);
	double result=romb.rombint(test,1.,2.);
	cout<<"result is "<<result<<endl;
	return 1;
}
*/
