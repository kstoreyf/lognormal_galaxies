#!/usr/bin/env python26
#============================================================================
# This is the Python script of calculating the power spectrum from an input
# correlation function. 
# The code uses a Romberg integration to calculate the power spectrum 
# from correlation function as
# P(k)=4*pi*int r^2*xi(r)*dr.
# The error in integration can be changed by changing the parameter in the
# integration romb(1.e-8) line in the C++ code. The default value is 1e-8.
#
#
# You need to provide the input correlation function, the no. of columns
# in the correlation function file, the column no. of r, xi(r), the linear
# bias to use and the value of rmax. The output file is named 
# pk_rmax%d_b%f.dat.
#
#
# 15 Mar 2016
# by Aniket Agrawal
#============================================================================

import os
import random
class executable:
	"""A class for running executables"""
	def __init__(self,exename):
		self.exename=exename
#		self.pfname='params.'+self.exename
	def __call__(self,params):
#		pfile=open(self.pfname,'a')
#		for param in params:
#			try:
#				pfile.write(param+'\n')
#			except:
#				if type(param)==int:
#					pfile.write('%d\n' % param)
#				else:
#					pfile.write('%g\n' % param)
#		pfile.close()
#		cmd = 'time ./'+self.exename+' <'+self.pfname
		nparams = len(params)
		paramfname = params[nparams-1]
		pfile=open(paramfname,'w')
		for iparam in range(0,nparams-1):
			try:
				pfile.write(params[iparam]+'\n')
			except:
				if type(params[iparam])==int:
					pfile.write('%d\n' % params[iparam])
				else:
					pfile.write('%g\n' % params[iparam])
		pfile.close()
		cmd = 'time ./'+self.exename+' <'+paramfname
		print cmd
		os.system(cmd)

xirfname = 'xi_gm_kmax30.dat'			#name of input correlation function file
ncol = 2                                        #no. of columns in power spectra files
rcol = 1                                        #column no. containing r (starting from 1)
xicol = 2                                       #column no. containing xi(r) (starting from 1)	
bias = 1.0			
rmax = 10000					#rmax used for integration
p1fname = 'params_calc_pk.dat'
executable('aux_codes/calc_pk')([xirfname, ncol, rcol, xicol, bias, rmax, p1fname])    
