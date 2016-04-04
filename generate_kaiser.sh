function calc_xi_gm {
KMAX=$1
aux_codes/calc_xi_gm << eof
pkG_rmax10000_b1.455.dat
pkG_rmax10000_b1.dat
2
1 2
$KMAX
eof
}

function calc_pk {					
RMAX=$1					
aux_codes/calc_pk << eof		
xi_gm_kmax19.dat
2
1 2
1.0
$RMAX					
eof						
}						

function discretize_pk {
LORDER=$1
aux_codes/discretize_pk << eof
pk_kaiser.dat
430
3
1 $LORDER
0.1 1.455
pkt_kaiser_$LORDER.dat
3666
1486
734
1024
eof
}

function kaiser_pk {
BIAS=$1
aux_codes/kaiser_pk << eof
wavenumber_pk.txt
1037
pk_rmax1500_b1.dat
430
$BIAS
wavenumber_fnu.txt
1037
eof
}

# generate the cross correlation function 
calc_xi_gm 19		# specify the redshift 	#modified by Aniket
# compute P(k) from xi(r)
calc_pk 1500
# compute Kaiser prediction
kaiser_pk 1.455
# discretize Kaiser P(k) monopole
discretize_pk 2 
# discretize Kaiser P(k) quadrupole
discretize_pk 3
