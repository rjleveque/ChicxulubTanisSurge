"""
Create the topography for a new test case.
"""
from pylab import *
from clawpack.geoclaw import nonuniform_grid_tools


h0 = 200.   # max depth of basin
h0r = 10. # depth of river
xrm = 0.  # location of river mouth
dxr = 10e3  # length of river
x2 = xrm + dxr  # right boundary
xshore1 = -10e3  # left shore location
xflat_start = xshore1 + 1.  # 1 meter gives nearly vertical wall
#xflat_start = xshore1 + 5000.  # for SI
slope1 = -h0/(xflat_start - xshore1)
B1 = 20  # elevation at x1
x1 = xshore1 + B1/slope1  # left boundary
xflat_end = xrm - 5e3
#slope2 = (h0r - h0) / (xflat_end - xrm)
#xflat_end = xrm - 5e3
slope2 = -h0 / (xrm - xflat_end)
B2 = -h0r 



xzpairs = [(x1, B1),                # left edge
           (xflat_start , -h0),   # start of flat
           (xflat_end , -h0),   # end of flat
           (xrm, -h0r),   # river mouth
           (x2, B2)]                # right edge

topo_fcn = nonuniform_grid_tools.make_pwlin_topo_fcn(xzpairs)

mx = 5000
hmin = 15.   # use uniform grid for shallower

nonuniform_grid_tools.make_celledges_cfl(x1,x2,mx,topo_fcn,
        hmin=hmin, fname='celledges.txt', plot_topo=True)


# compute travel time from shore to shore by integrating 1/sqrt(g*h):

celldata = loadtxt('celledges.txt',skiprows=1)
ttime = 0
for k in range(celldata.shape[0]-1):
    if celldata[k,1] < 0:
        ttime += (celldata[k+1,0]-celldata[k,0])/sqrt(-9.81*celldata[k,1])
print('Travel time across domain = %.2f second' % ttime)
