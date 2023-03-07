#!/bin/bash
#PBS -N SIDM_test2
#PBS -l nodes=1:ppn=16
#PBS -j oe
#PBS -o outputs/multi_SIDM_test2_16mergerTree.log
#PBS -V
#PBS -m abe
#PBS -M nahvazi@carnegiescience.edu

cd $PBS_O_WORKDIR
export OMP_THREAD_NUM=16
/usr/bin/time -v Galacticus.exe parameters/myParams/multi_SIDM_test_16mergerTree.xml


