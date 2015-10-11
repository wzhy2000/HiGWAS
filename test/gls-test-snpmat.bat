#!/bin/bash 

#----------------------------------------------------
# Generic SLURM script 
#----------------------------------------------------

#SBATCH -J gls.test.snpmat # Job name
#SBATCH -o gls.test.snpmat.%j.out # stdout; %j expands to jobid
#SBATCH -e gls.test.snpmat.%j.err # stderr; skip to combine stdout and stderr
#SBATCH -p normal # queue
#SBATCH -N 1 -n 8 # one node and one task
#SBATCH -t 48:00:00 # max time
#SBATCH --mail-user=zw355@cornell.edu
#SBATCH --mail-type=ALL

$R --vanilla --quiet  < gls-test-snpmat.r > gls-test-snpmat.out 
