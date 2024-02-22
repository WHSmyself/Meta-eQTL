#!/bin/bash -xe
#PBS -N gene_cluster
#PBS -e /job_error/
#PBS -o /job_out/
#PBS -q middle
#PBS -l nodes=1:ppn=1
#PBS -l walltime=999:00:00
#PBS -t 0
# Kill script if any commands fail
Rscript gene_cluster.R -i data.txt -g LYC00211 -t 0.5 -e 1 -o /data1/test/ 
