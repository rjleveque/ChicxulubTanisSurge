#! /bin/sh

# Use the full path as per your computer of nms
NMS_bin_path=path_of_modeslib_bin_in_your_system

# input file: nms.dat,receivers.dat & sources.dat

${NMS_bin_path}/nms<< !
nms.dat
receivers.dat
sources.dat
!
