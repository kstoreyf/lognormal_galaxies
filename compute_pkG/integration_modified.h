#ifndef INTEGRATION_H_
#define INTEGRATION_H_
#include <iostream>
#include <cmath>

// This class implement the two different integration methods:
// 1. qromb
// 2. rk5int

using namespace std;

class integration{
private:
	double tol;		// error tolerance

	template<class T>
	void rkck(const double y,const double dydx,const double x,const double h,
				double &yout,double &yerr,T &func,int dim,double params[]);
	template<class T>
	void rkqs(double &y,const double dydx,double &x,const double htry,
				const double yscal,double &hdid,double &hnext,T &func,int dim,double params[]);
public:
	integration(double atol):tol(atol){}
	template<class T>
	double rombint(T &func, const double a, const double b,int dim,double params[]);
	template<class T>
	double rk5int(T &func, const double a, const double b,int dim,double params[]);
};

#endif
