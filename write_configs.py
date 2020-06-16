

def main():

    L = 750
    N = 125000
    seed = 0
    fn = f'configs/cat_lognormal_L{L}_N{int(N/1000)}k.cfg'
    params = set_params(L=L, N=N, seed=seed)
    write(params, fn)

def set_params(L=None, N=None, seed=None):
    params = {}
    if seed:
        params['seed'] = seed
        params['Nrealization'] = 1
    if L:
        params['Lx'] = L    
        params['Ly'] = L
        params['Lz'] = L
    if N:
        params['Ngalaxies'] = N
    return params

def write(params, fn):
    
    final_params = default_params.copy()

    for p, val in params.items():
        assert p in default_params, f"Parameter {p} not recognized!"
        final_params[p] = val
    
    with open(fn, "w") as f:
        for p, val in final_params.items():
            if val==None:
                continue
            f.write(f"{p} = {val}")
            f.write("\n")

    return            



default_params = {'out_dir': '../catalogs',
        'z': 0.0, #survey
        'bias': 1.5,
        'Ngalaxies': 125000,
        'Lx': 750,
        'Ly': 750,
        'Lz': 750,
        'Nrealization': 0, # number or realization
        'seed': 0,
        'oc0h2': 0.11880, # \Omega_c h^2 # PLANCK
        'mnu': 0.00,   # \Sigma m_{\nu} total neutrino mass
        'ns': 0.9667,
        'lnAs': 3.064,
        'ob0h2': 0.02230, #\Omega_baryon h^2
        'h0': 0.6774, # H0/100
        'w': -1.0,
        'run': 0.0,
        'inp_pk_fname': None,
        gen_inputs = True
run_lognormal = True
calc_pk = True
calc_mode_pk = 0 # 0:rectangular, 1:cubic
Nrealization = 0 # number or realization
seed = 0 # seed for random numbers 
Pnmax = 512 # mesh number
losx = 0
losy = 0
losz = 0
kbin = 0.01 # bin width of Pk
kmax = 0 # maximum k. if 0, automatically set to be nyquist frequency
lmax = 4         
num_para = 24 # parallel number

halofname_prefix = 

# set 1 if you want to use cross pk as an input
use_cpkG = 0
# should be less than ~90% of galaxy bias
bias_cpkG = 1.35 #?? what if galaxy bias is 1?

output_gal = 1 # 0: don't output gal catalog, 1: output gal catalog
output_matter = 0 # 0: don't output matter density file, 1: output matter density file
}

if __name__=='__main__':
    main()