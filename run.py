#!/usr/bin/env python
#============================================================================
# This is the Python script of running the following steps of the simulation:
# [1] generating input files (power spectrum, 2-point correaltion function and 
#     the power-spectrum of the log-normal field)
# [2] generating mock galaxy catalog from the input files calculated in [1]
# [3] calculating the power spectrum from the mock galaxy catalog generated in [2]
# 
# To run the code, type:
#
# > ./run.py example.ini
#
# where exmaple.ini is the exmaple of input file.
# It contains the input parameters of simulation e.g. cosmological parameters,  
# survey geometries, etc.
#
# The resulting files are as follows:
# 1) input files
# ./data/inputs/example_pk.txt
# ./data/inputs/exmaple_pkG_b1.5.txt
# ...
# 2) mock galaxy catalog
# For each realization, the code generates the followng two files:
# ./data/lognormal/example_lognormal_rlz{# of realization}.bin (for galaxy field)
# ./data/lognormal/example_density_lognormal_rlz{# of realization}.bin (for density field)
# 3) Power spectrum
# ./data/pk/exmaple_pk_rlz{# of realization}.dat
# The column 1, 2, 3 and 4 are k [Mpc/h], P_0, P_2 and P_4 [Mpc^3/h^3], respectively.
#
# For the detailes of the calculations of each step, please see readme.txt.
#   
# 04 Sep 2016
# by Ryu Makiya
#============================================================================

import sys
import random
import numpy as np
from multiprocessing import Pool
from multiprocessing import Process
from run_module import *

# read parameters from .ini file
ini_file_name = sys.argv[1]
params = read_params(ini_file_name)

# check whether output directory exists or not
check_dir(params)

# create seeds for random
random.seed(params['seed'])
seed1 = []
seed2 = []
seed3 = []

for i in range(0,params['Nrealization']):
    seed1.append(random.randint(1,100000))
    seed2.append(random.randint(1,100000))
    seed3.append(random.randint(1,100000))

exe = executable("exe")

# generate input files
if params['gen_inputs']:
    # input power spectrum
    params['ofile_eh'] = params['out_dir']+'/inputs/'+params['ofile_prefix']
    args = ['ofile_eh','om0','ode0','ob0','h0','w','ns','run','As','mnu','z'] # do not change the order
    exe.run('eisensteinhubaonu/compute_pk',args,params)

    # xi
    params['ofile_xi'] = params['out_dir']+'/inputs/'+params['ofile_prefix']
    params['len_inp_pk'] = sum(1 for line in open(params['inp_pk_fname']))
    args = ['ofile_xi','inp_pk_fname','len_inp_pk'] # do not change the order
    exe.run('compute_xi/compute_xi',args,params)

    # galaxy pkG
    params['ncol'] = np.size(np.loadtxt(params['xi_fname'])[0,:])
    args = ['pkg_fname','xi_fname','ncol','bias','rmax'] # do not change the order
    exe.run('compute_pkG/calc_pkG',args,params)

    # matter pkG
    args = ['mpkg_fname','xi_fname','ncol','bias_mpkG','rmax'] # do not change the order
    exe.run('compute_pkG/calc_pkG',args,params)

    # galaxy-matter pkG
    if params['use_cpkG'] == 1: 
        params['bias_cpkG'] = np.sqrt(params['bias_cpkG'])
        args = ['cpkg_fname','xi_fname','ncol','bias_cpkG','rmax'] # do not change the order
    exe.run('compute_pkG/calc_pkG',args,params) #idk about indentation?
    params['bias_cpkG'] = (params['bias_cpkG'])**2.0
    print("DONE HERE")
else:
    print 'skip: generating input files'

# generate Poisson catalog
if params['run_lognormal']:
    print("RUN")
    def gen_Poisson(i):
        params_tmp = params
        params_tmp['seed1'] = seed1[i]
        params_tmp['seed2'] = seed2[i]
        params_tmp['seed3'] = seed3[i]
        
        # output file names
        params_tmp['Poissonfname'] = params['out_dir']+'/lognormal/'+params['ofile_prefix']\
                                    +'_lognormal_rlz'+str(i)+'.bin'
        params_tmp['Densityfname'] = params['out_dir']+'/lognormal/'+params['ofile_prefix']\
                                    +'_density_lognormal_rlz'+str(i)+'.bin'

        if os.path.exists(params_tmp['Poissonfname']):
            print(params_tmp['Poissonfname'], "already exists!")
            return

        args = ['pkg_fname','mpkg_fname','use_cpkG','cpkg_fname','Lx','Ly','Lz','Pnmax',\
                        'Ngalaxies','aH','f_fname','bias','seed1','seed2','seed3','Poissonfname',\
                        'Densityfname','output_matter','output_gal'] # do not change the order
#       args = ['pkg_fname','mpkg_fname','use_cpkG','cpkg_fname','Lx','Ly','Lz','Pnmax',\
#                        'Ngalaxies','aH','f_fname','bias','seed1','seed2','seed3','Poissonfname',\
#                        'Densityfname'] # do not change the order
        print("running")
        exe.run('generate_Poisson/gen_Poisson_mock_LogNormal',args,params_tmp)
    p = Pool(params['num_para']) 
    run = p.map(gen_Poisson,range(params['Nrealization']))
else:
    print 'skip: generating log-normal catalog'

# calculate Pk
if params['calc_pk']:
    def calc_Pk(i):
        params_tmp = params

        # input file names
        if (params_tmp['halofname_prefix']) == '':
            params_tmp['halofname'] = params['out_dir']+'/lognormal/'+params['ofile_prefix']+'_lognormal_rlz'+str(i)+'.bin'
        else:
            params_tmp['halofname'] = params['out_dir']+'/lognormal/'+params['halofname_prefix']+'_lognormal_rlz'+str(i)+'.bin'
        
        if (params_tmp['imul_fname']) == '':
            params_tmp['imul_fname'] = params['out_dir']+'/coupling/'+params['ofile_prefix']+'_coupling.bin'
        params_tmp['pk_fname'] = params['out_dir']+'/pk/'+params['ofile_prefix']+'_pk_rlz'+str(i)+'.dat'

        args = ['halofname','Pnmax','aH','losx','losy','losz','kbin','kmax','lmax','imul_fname','pk_fname','calc_mode_pk']
        exe.run('calculate_pk/calc_pk_const_los_ngp',args,params_tmp)

    p = Pool(params['num_para']) 
    run = p.map(calc_Pk,range(params['Nrealization']))
else:
    print 'skip: calculate Pk'

# calculate Pk
if params['calc_cpk']:
    def calc_cPk(i):
        params_tmp = params

        # input file names
        if (params_tmp['halofname_prefix']) == '':
            params_tmp['halofname1'] = params['out_dir']+'/lognormal/'+params['ofile_prefix']+'_lognormal_rlz'+str(i)+'.bin'
        else:
            params_tmp['halofname1'] = params['out_dir']+'/lognormal/'+params['halofname_prefix']+'_lognormal_rlz'+str(i)+'.bin'
        params_tmp['halofname2'] = params['out_dir']+'/lognormal/'+params['ofile_prefix']+'_density_lognormal_rlz'+str(i)+'.bin'
    
        if (params_tmp['imul_fname']) == '':
            params_tmp['imul_fname'] = params['out_dir']+'/coupling/'+params['ofile_prefix']+'_coupling.bin'
        params_tmp['cpk_fname'] = params['out_dir']+'/pk/'+params['ofile_prefix']+'_cpk_rlz'+str(i)+'.dat'

        params_tmp['tmp1'] = 0
        params_tmp['tmp2'] = 1
        
        args = ['halofname1','halofname2','Pnmax','aH','losx','losy','losz','kbin','kmax','lmax','imul_fname','cpk_fname','calc_mode_pk','tmp1','tmp2']
        exe.run('calculate_cross/calc_cpk_const_los_v2',args,params_tmp)

    p = Pool(params['num_para']) 
    run = p.map(calc_cPk,range(params['Nrealization']))
else:
    print 'skip: calculate cPk'
