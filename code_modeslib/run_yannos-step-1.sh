#!/bin/sh
# Inputs are Earth Model: PREMQL6 and yannos file: yannos.dat
# Use the full path as per your computer of yannos or yannos_MPI
NMS_bin_path=path_of_modeslib_bin_in_your_system
# for prallel computation
#time mpiexec -n 6 ${NMS_bin_path}/yannos_MPI
# for non-parallel computation
time ${NMS_bin_path}/yannos_MPI
