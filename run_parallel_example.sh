#!/bin/sh
#$ -cwd
#$ -j y
## -M komatsu@mpa-garching.mpg.de
## -m e
#$ -N calc_pk
#$ -l h_rt=24:00:00
#$ -pe impi_hydra 16

PATH=/bin:/sbin:/usr/bin:/usr/sbin:/usr/local/bin:/usr/local/sbin:~/bin
export PATH

module load intel

echo "50 1" | python run_genPoissonmock_iseed.py
