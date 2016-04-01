#!/usr/bin/env python26
#============================================================================
# This is the Python script of calculating the discretized power spectrum. 
# The discretized power spectrum is the power spectrum calculated at the 
# exact same values as are used while generating the fields. When we measure
# the power spectrum we average over a spherical shell. This code does the 
# same for a theoretical power spectrum. The resulting discrete power 
# spectrum provides a more accurate match to the measured power spectrum.
#
# You need to provide the name of the theory power spectrum, (the length is
# calculated in this code itself), the smoothing scale (if used to smooth 
# the power spectrum), linear bias, the name of the output file, the length
# of the box used and no. of points used to generate the mesh while 
# generating density fields. 
# NOTE that this code assumes a standard k binning with fundamental
# frequency spacing in each direction. To use a non standard grid, please
# modify the code.  
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

#-----------------------------------------------------------------------------
# function to read number of lines in the file
#-----------------------------------------------------------------------------
def file_len(fname):
	f=open(fname)
	for i, l in enumerate(f):
		pass
	f.close()
	return i + 1

pkfname = 'wavenumber_pk.txt'			#name of theory P(k) file
npk = file_len(pkfname)				#length of above file
print 'number of data in '+pkfname+' is '+str(npk)
Rs = 0.1                                        #smoothing radius if used for smoothing power spectrum
bias = 1.455                                   
ofname = 'discrete_'+pkfname                   	
Lx = 3666					
Ly = 1486
Lz = 734
nmesh = 1024					#please use the same as used for measuring P(k)
p1fname = 'params_discretize_pk.dat'
executable('aux_codes/discretize_pk')([pkfname, npk, Rs, bias, ofname, Lx, Ly, Lz, nmesh, p1fname])    