function compute_pk {
REDSHIFT=$1
eisensteinhubaonu/compute_pk << eof
$REDSHIFT
eof
}

function compute_pk2 {
REDSHIFT=$1
eisensteinhubaonu/compute_pk2 << eof
$REDSHIFT
eof
}

function compute_xi {
compute_xi/compute_xi
}

function compute_pkG {
BIAS=$1
compute_pkG/calc_pkG << eof
Rh_xi.txt
2
1 2
$BIAS
10000
eof
}

# generate the linear matter power spectrum, P(k), using Eisenstein&Hu and the logarithmic growth rate f(k)
compute_pk2 1.3		# specify the redshift 
# compute 2-point correlation function, xi(r), from P(k)
compute_xi
# compute the power spectrum of log-normal field from xi(r) 
compute_pkG 1.455	# specify the linear galaxy bias
compute_pkG 1.0		# compute the power spectrum of matter field to generate velocities
