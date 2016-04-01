#!/usr/bin/env python26
#============================================================================
# This is the Python script of calculating the power spectrum from the
# mock galaxy catalogue generated in /genmock folder.
#
# First parts calculate the number of galaxyes N_g, and the number of 
# random samples N_s inside of a grid.
# [1] calc_Fourier_box::
#		calculate the box size, and 
# [2] galaxy_to_grid::
#		put galaxies into the Fourier grid
# [3] calc_Luminosity_fn::
#		code to estimate n(L) from given galaxy sample.
#
# the resulting files are as follow:
# 1) size of the Fourier box
# ../Grids/fourier_dimension_[selection name]_[survey name].dat 
# 2) galaxy numbers in the grid
# ../Grids/gpicc_[survey name]_[selection name]_nmesn=[nmesh].bin 
# 3) galaxy luminosity function n(L)
# ../Grids/gpicc_[survey name]_nmesn=[nmesh].bin 
# 
# The second part generates the random sample inside of the box, 
# applies exactly the same selection as LAEs, and assign the numbers to
# the Fourier grid.
# [1] rsample_to_grid::
#		Select random samples for given selection scheme
#		selection=0	::	no selection
#		selection=1	::	geometry
#		selection=2	::	shot only
#		selection=3	::	shot+mask
#		selection=4	::	radial only
#		selection=5	::	shot+mask+radial
#
# the resulting files are:
# 1) number of random smaples in the grid
# ../Grids/random_[survey name]_[selection name]_rspace.bin
#
# In the third part, we estimate the power spectrum of LAEs from the 
# number of LAEs/random samples in grids.
#
# 24 July 2011
# by Donghui Jeong
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

#-----------------------------------------------------------------------------
# how many realisation? initial seed?
#-----------------------------------------------------------------------------
#nrealisation = int(raw_input('N_realisation: '))
#iseed = int(raw_input('iseed: '))
tmp = raw_input('N_realisation & iseed [separated by space]: ')
nrealisation = int(tmp.split()[0])
iseed = int(tmp.split()[1])

print 'the number of realisation = ',nrealisation
print 'the initial seed = ',iseed

#-----------------------------------------------------------------------------
# set default values 
#-----------------------------------------------------------------------------
# cosmological parameters
Omegam=0.272		# Omega matter for computing the velocity field
Omegade=1-Omegam
zz=1.3			# redshift for computing the velocity field
aHz=100.*pow(Omegam*pow(1.+zz,3)+Omegade,0.5)/(1.+zz)

# parameters for lognormal realization
bias = 1.455		# linear bias for computing the matter fluctuation and then velocity field (need to be consistent with input P(k))

Lx = 3666	# box size in x [h^-1 Mpc]
Ly = 1486	# box size in y [h^-1 Mpc]
Lz = 734	# box size in z [h^-1 Mpc]
#Pnmax = 2500	# mesh number of the longest dimension for lognormal realization
Pnmax = 1024


# file names
inputpkfname = 'pkG_rmax10000_b'+str(bias)+'.dat'	# input log-transformed power spectrum (format from calc_pkG)
num_data_inputpk = file_len(inputpkfname)
print 'number of data in '+inputpkfname+' is '+str(num_data_inputpk)

nmax = Pnmax	# mesh number of the longest dimension for computing power spectrum
kbin = 0.005	# kbin for the power spectrum (0 for the largest fundamental frequency)
kmax = 0	# kmax for the output power spectrum (0 for Nyquist frequency)
losx = 0	# line-of-sight direction of x (this is a vector, and enter (0,0,0) for real-space calculation)
losy = 0	# line-of-sight direction of y
losz = 0	# line-of-sight direction of z

lmax=4          # lmax for reconstructing the leakage (the output pk file will contain k, pmono, pquad, and phexa, so lmax>=4)
imulfname='coupling/imul_Lx'+str(Lx)+'_Ly'+str(Ly)+'_Lz'+str(Lz)+'_nmax'+str(nmax)+'_kbin'+str(kbin)+'_kmax'+str(kmax)+'_lmax'+str(lmax)+'.bin'     # input inverse mu-leakage matrix file; if it does not exist, the code will compute and output it

random.seed(iseed)
print 'iseed=',iseed
for i in range(0,nrealisation):
	seed1 = random.randint(1,100000)
	seed2 = random.randint(1,100000)
	seed3 = random.randint(1,100000)
        p2fname = 'pk/params.calc_pk_const_los_is'+str(iseed)+'_rlz'+str(i)
	Poissonfname = 'lognormal/lognormal_is'+str(iseed)+'_rlz'+str(i)+'_Pnmax'+str(Pnmax)+'_b'+str(bias)+'.bin'	#output lognormal galaxy catalog in (x,y,z,vx,vy,vz) [h^-1 Mpc] and [km/s]
        pkfname='pk/pk_ngp_is'+str(iseed)+'_rlz'+str(i)+'_z'+str(zz)+'_b'+str(bias)+'_los'+str(losx)+'.'+str(losy)+'.'+str(losz)+'.dat'
        executable('calculate_pk/calc_pk_const_los_ngp')([Poissonfname,nmax,aHz,losx,losy,losz,kbin,kmax,lmax,imulfname,pkfname,p2fname])     # same as before, but with two more inputs (lmax and imulfname)
