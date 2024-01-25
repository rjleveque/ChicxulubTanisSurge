Code to run GeoClaw tsunami simulations, provided by Randy LeVeque

Commit cc5634177 of the clawpack/geoclaw master branch was used for these
simulations.   At the time done, the 1D GeoClaw code was not yet in the
latest release v5.9.2 of GeoClaw, but should appear soon in v5.10.0, at
which point these codes should run with that release.

Standard GeoClaw codes to set up the problem (see the documentation):

    Makefile
    claw1ez.f
    b4step1.f90
    src1.f90
    setprob_module.f90
    qinit.f90

    make_celledges.py
    setrun.py
    setplot.py
    b4run.py

The routine src1.f90 adds horiztonal acceleration as a source term,
based on the seismic signal in:

    Radial_accel_R032_HD60s_source-depth_10km_filter_T40-500s.txt

The geometry of the basin and depth of the river are specified in
make_celledges.py. The code was run 4 times after re-executing this
script for the 4 cases shown in the paper and SI.

Scripts to make figures and animations:

    make_plots_for_paper.py
    make_anims.py
