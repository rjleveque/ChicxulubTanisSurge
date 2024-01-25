Code provided by Satish Maura

Programs to compute synthetic seismograms with normal mode summation
in spherically symmetric Earth models.

Details of the Program can be found here, which includes installation
instruction and background of the input parameters.

https://gitlab.univ-nantes.fr/capdeville-y/modeslib/

After installation successfully done.

- Remember to replace `/path/to/modeslib/bin` with the actual path where the modeslib binaries are located on your system. 
Users should follow the installation instructions provided in the link before executing these scripts.

1. run_yannos-step-1.sh
2. run_nms_step-2.sh


###########
Script-1: run_yannos-step-1.sh
###########

# Set the path to the modeslib binary on your system after successful installation
MODESLIB_BIN_PATH=/path/to/modeslib/bin

# For parallel computation
#time mpiexec -n 6 ${MODESLIB_BIN_PATH}/yannos_MPI

# For non-parallel computation
time ${MODESLIB_BIN_PATH}/yannos_MPI


###########
Script-2: run_nms_step-2.sh
###########

# Set the path to the modeslib binary on your system after successful installation
MODESLIB_BIN_PATH=/path/to/modeslib/bin

# Specify input files: nms.dat, receivers.dat, and sources.dat
${MODESLIB_BIN_PATH}/nms << !
nms.dat
receivers.dat
sources.dat
!

############
Output: UR_R032 UZ_R032 UT_R032 (synthetic seismogram for Radial, Vertical, and Transversal components)
############


- Maintained consistent formatting especially for receivers.dat and sources.dat file.

