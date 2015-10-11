#!/bin/bash 

#----------------------------------------------------
# Generic SLURM script 
#----------------------------------------------------

#SBATCH -J bls.test.simple # Job name
#SBATCH -o bls.test.simple.%j.out # stdout; %j expands to jobid
#SBATCH -e bls.test.simple.%j.err # stderr; skip to combine stdout and stderr
#SBATCH -p normal # queue
#SBATCH -N 1 -n 8 # one node and one task
#SBATCH -t 48:00:00 # max time
#SBATCH --mail-user=zw355@cornell.edu
#SBATCH --mail-type=ALL

$R --vanilla --quiet  < bls-test-simple.r > bls-test-simple.out 

