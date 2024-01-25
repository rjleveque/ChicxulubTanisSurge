"""
Controls plotting of the time frames.
Only figno=2 is used for plots in the paper.
Adjust plotdata.print_fignos to print others that may be useful to see.
"""

import os, sys
from scipy.interpolate import interp1d


from clawpack.visclaw import geoplot


from clawpack.geoclaw.nonuniform_grid_tools import make_mapc2p
import numpy

try:
    fname = '_output/fgmax.txt'
    d = numpy.loadtxt(fname)
    etamax = numpy.where(d[:,1]>1e-6, d[:,3], numpy.nan)
    xmax = d[:,0]
    jmax = numpy.where(d[:,1]>0)[0].max()
    print("run-in = %8.2f,  run-up = %8.2f" % (d[jmax,0],d[jmax,3]))
    print('Loaded hmax from ',fname)
except:
    xmax = None
    print("Failed to load fgmax.txt")
    
xmax = None
print("Not plotting fgmax.txt")



#xlimits = [-20,70]
xlimits = [-6e3,6e3]
ylimits = [-110,20]


if 1:
    fforce = 'Radial_accel_Filter-period_40-500s_R032_Chicxulub_100_-Eq-Mw11.txt'
    d = numpy.loadtxt(fforce, skiprows=1)
    xforce = interp1d(d[:,0],d[:,1],kind='linear',bounds_error=False,
                      fill_value=0.)
                      #fill_value=numpy.nan)

if 0:
    def xforce(t):
        from numpy import arctan,sin
        amp_xforce = 0.020
        period_force = 200. #130.  #1500.
        omega_force = 8.*arctan(1.) / period_force
        return amp_xforce * sin(t*omega_force)

def timeformat(t):
    from numpy import mod
    hours = int(t/3600.)
    tmin = mod(t,3600.)
    min = int(tmin/60.)
    sec = int(mod(tmin,60.))
    timestr = '%s:%s:%s' % (hours,str(min).zfill(2),str(sec).zfill(2))
    return timestr

def title_hours(current_data, name=None):
    from pylab import title
    t = current_data.t
    timestr = timeformat(t)
    if name:
        title('%s at time %s after impact' % (name,timestr))
    else:
        title('%s after impact' % timestr)


def setplot(plotdata):

    plotdata.clearfigures()

    fname_celledges = os.path.join(plotdata.outdir,'celledges.txt')
    mapc2p1, mx_edge, xp_edge = make_mapc2p(fname_celledges)
    
    try:
        fname = os.path.join(plotdata.outdir, 'salt.txt')
        d = numpy.loadtxt(fname)
        tsalt = d[:,0]
        xsalt = d[:,1]
        print('Loaded xsalt from %s, min = %.1f, max = %.1f' \
                % (fname, xsalt.min(),xsalt.max()))
    except:
        xsalt = None
        print("Failed to load salt.txt")

    def fixticks_vel(current_data):
        from pylab import ticklabel_format, grid,tight_layout,plot, where
        ticklabel_format(useOffset=False)
        grid(True)
        title_hours(current_data, 'Velocity')
        tight_layout()

    def fixticks(current_data):
        from pylab import ticklabel_format, plot,grid,gca, where,text
        ticklabel_format(useOffset=False)
        title_hours(current_data, 'Surface')
        if xmax is not None:
            plot(xmax, etamax, 'r', linewidth=0.8)
        if xsalt is not None:
            jsalt = where(tsalt <= current_data.t)[0].max()
            #print('+++ jsalt = ',jsalt)
            xs = xsalt[jsalt]
            ysalt = (-10,8)
            plot([xs,xs],ysalt,'k--')
            text(xs, ysalt[1]*1.3, 'salt/fresh\n%.0f m' % xs, 
                 horizontalalignment='center', fontsize=11,
                 backgroundcolor='yellow', alpha=1)
        grid(True)
        
    def velocity(current_data):
        from pylab import where
        q = current_data.q
        u = where(q[0,:]>1e-3, q[1,:] / q[0,:], 0.)
        return u


    plotfigure = plotdata.new_plotfigure(name='domain', figno=0)
    plotfigure.kwargs = {'figsize':(8,7)}
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(311)'
    plotaxes.xlimits = xlimits
    #plotaxes.ylimits = [-15,15]
    #plotaxes.ylimits = [-3,3]
    plotaxes.ylimits = [-40,40]
    plotaxes.title = 'Surface displacement'
    plotaxes.afteraxes = fixticks

    plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    #plotitem.show = False
    plotitem.plot_var = geoplot.surface
    plotitem.plot_var2 = geoplot.topo
    plotitem.color = [.7,.7,1]
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.surface
    plotitem.color = 'b'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    #plotitem.show = False
    plotitem.plot_var = geoplot.topo
    plotitem.color = 'g'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(312)'
    plotaxes.xlimits = xlimits
    #plotaxes.ylimits = [-10,10]
    plotaxes.ylimits = [-20,20]
    #plotaxes.title = 'Velocity'
    plotaxes.afteraxes = fixticks_vel
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.u_velocity
    plotitem.color = 'b'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(313)'
    plotaxes.xlimits = xlimits
    #plotaxes.ylimits = [-110,30]
    plotaxes.ylimits = ylimits
    plotaxes.title = 'Full depth'

    def add_force(current_data):
        from pylab import ticklabel_format, grid,tight_layout, arrow, text
        t = current_data.t
        f = xforce(t)
        fscale = (xlimits[1]-xlimits[0]) 
        #print('arrow length = ',f*fscale)
        xloc = (xlimits[0] + xlimits[1])/2.
        yloc = (ylimits[0] + ylimits[1])/2.
        #if abs(f) > 1e-6:
        if 0:
            arrow(xloc, yloc, f*fscale, 0., width=0.003, 
                  #head_width=0.006*f*fscale, 
                  #head_length=0.5*f*fscale,
                  fc='r', ec='r')
        if 1:
            text(4e3,-50,'seismic acceleration \n%.4f m/s^2' % f, 
                 color='r', horizontalalignment='center',fontsize=12)
        ticklabel_format(useOffset=False)
        title_hours(current_data, 'Full depth')
        grid(True)
        tight_layout()


    plotaxes.afteraxes = add_force
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1
    
    plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    #plotitem.show = False
    plotitem.plot_var = geoplot.surface
    plotitem.plot_var2 = geoplot.topo
    plotitem.color = [.7,.7,1]
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.surface
    plotitem.color = 'b'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    #plotitem.show = False
    plotitem.plot_var = geoplot.topo
    plotitem.color = 'g'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1
    #----------

    plotfigure = plotdata.new_plotfigure(name='water', figno=2)
    plotfigure.kwargs = {'figsize':(10,3)}
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [-10.2e3,6e3] #xlimits
    #plotaxes.ylimits = [-110,20]
    plotaxes.ylimits = [-20,20]
    plotaxes.title = 'Surface displacement'
    
    def fixticks3(current_data):
        from pylab import ticklabel_format, plot,grid,gca, where,text,title,\
             xticks,yticks,xlabel,ylabel
        ticklabel_format(useOffset=False)
        #tminutes = current_data.t / 60.
        #title('Water surface at time t = %.2f minutes' % tminutes)
        title_hours(current_data, 'Surface')
        if xmax is not None:
            plot(xmax, etamax, 'r', linewidth=0.8)
        if xsalt is not None:
            jsalt = where(tsalt <= current_data.t)[0].max()
            #print('+++ jsalt = ',jsalt)
            xs = xsalt[jsalt]
            #plot([xs,xs],[-20,10],'k--')
            #text(xs-1e3, 12, '%.0f m' % xs, fontsize=9)
            ysalt = (-10,13)
            plot([xs,xs],ysalt,'k--')
            text(xs, ysalt[1]*1., 'salt/fresh\n%.0f m' % xs, 
                 horizontalalignment='center', fontsize=11,
                  alpha=1,
                 bbox={'facecolor':'yellow','edgecolor':'k','lw':0.5})
        grid(True)
        xlabel('meters', fontsize=12)
        ylabel('meters', fontsize=12)
        xticks(fontsize=10)
        yticks(fontsize=10)
        
    plotaxes.afteraxes = fixticks3
    
    def saltwater(current_data):
        from pylab import where, nan
        if xsalt is not None:
            jsalt = where(tsalt <= current_data.t)[0].max()
            #print('+++ jsalt = ',jsalt)
            xs = xsalt[jsalt]
        else:
            xs = -1e6
        x = mapc2p1(current_data.x)
        eta = current_data.q[-1,:]
        swater = where(x<xs, eta, nan)
        #print('+++ saltwater')
        #import pdb; pdb.set_trace()
        return swater
    
    def freshwater(current_data):
        from pylab import where, nan
        if xsalt is not None:
            jsalt = where(tsalt <= current_data.t)[0].max()
            #print('+++ jsalt = ',jsalt)
            xs = xsalt[jsalt]
        else:
            xs = -1e6
        x = mapc2p1(current_data.x)
        eta = current_data.q[-1,:]
        fwater = where(x>=xs, eta, nan)
        #print('+++ freshwater')
        #import pdb; pdb.set_trace()
        return fwater    

    plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    #plotitem.show = False
    plotitem.plot_var = saltwater
    plotitem.plot_var2 = geoplot.topo
    #plotitem.color = [.6,.6,1]
    plotitem.color = [1,.6,.6]

    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    #plotitem.show = False
    plotitem.plot_var = freshwater
    plotitem.plot_var2 = geoplot.topo
    plotitem.color = [.8,.8,1]
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1
    
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.surface
    plotitem.color = 'b'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    #plotitem.show = False
    plotitem.plot_var = geoplot.topo
    plotitem.color = 'g'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    #----- figure for topo and initial eta -----


    plotfigure = plotdata.new_plotfigure(name='initial', figno=22)
    plotfigure.kwargs = {'figsize':(10,6)}
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [-10e3,6e3] #xlimits
    #plotaxes.ylimits = [-110,20]
    plotaxes.ylimits = [-220,20]
    plotaxes.title = ''
    def fixticks4(current_data):
        fixticks3(current_data)
        
    plotaxes.afteraxes = fixticks4


    plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    #plotitem.show = False
    plotitem.plot_var = saltwater
    plotitem.plot_var2 = geoplot.topo
    #plotitem.color = [.6,.6,1]
    plotitem.color = [1,.6,.6]

    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    #plotitem.show = False
    plotitem.plot_var = freshwater
    plotitem.plot_var2 = geoplot.topo
    plotitem.color = [.8,.8,1]
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1
    
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.surface
    plotitem.color = 'b'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    #plotitem.show = False
    plotitem.plot_var = geoplot.topo
    plotitem.color = 'g'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1
    
        
    #----------
    plotfigure = plotdata.new_plotfigure(name='domain2', figno=3)
    plotfigure.kwargs = {'figsize':(8,6)}
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = [-15,15]
    plotaxes.title = 'Surface displacement'
    plotaxes.afteraxes = fixticks

    plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    #plotitem.show = False
    plotitem.plot_var = saltwater
    plotitem.plot_var2 = geoplot.topo
    plotitem.color = [.5,.4,1]
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    #plotitem.show = False
    plotitem.plot_var = freshwater
    plotitem.plot_var2 = geoplot.topo
    plotitem.color = [.8,.8,1]
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1
    

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.surface
    plotitem.color = 'b'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    #plotitem.show = False
    plotitem.plot_var = geoplot.topo
    plotitem.color = 'g'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1


    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.xlimits = xlimits
    #plotaxes.ylimits = [-110,30]
    plotaxes.ylimits = ylimits
    plotaxes.title = 'Full depth'

    def add_force(current_data):
        from pylab import ticklabel_format, grid,tight_layout, arrow, text
        t = current_data.t
        f = xforce(t)
        fscale = (xlimits[1]-xlimits[0]) 
        #print('arrow length = ',f*fscale)
        xloc = (xlimits[0] + xlimits[1])/2.
        yloc = (ylimits[0] + ylimits[1])/2.
        #if abs(f) > 1e-6:
        if 0:
            arrow(xloc, yloc, f*fscale, 0., width=0.003, 
                  #head_width=0.006*f*fscale, 
                  #head_length=0.5*f*fscale,
                  fc='r', ec='r')
        if 0:
            text(4e3,-50,'seismic acceleration \n%.4f m/s^2' % f, 
                 color='r', horizontalalignment='center',fontsize=12)
        ticklabel_format(useOffset=False)
        title_hours(current_data, 'Full depth')
        grid(True)
        tight_layout()


    plotaxes.afteraxes = add_force
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1
    
    plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    #plotitem.show = False
    plotitem.plot_var = saltwater
    plotitem.plot_var2 = geoplot.topo
    plotitem.color = [.5,.4,1]
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    #plotitem.show = False
    plotitem.plot_var = freshwater
    plotitem.plot_var2 = geoplot.topo
    plotitem.color = [.8,.8,1]
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.surface
    plotitem.color = 'b'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    #plotitem.show = False
    plotitem.plot_var = geoplot.topo
    plotitem.color = 'g'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1
    #----------


    plotfigure = plotdata.new_plotfigure(name='shore', figno=1)
    plotfigure.kwargs = {'figsize':(8,7)}
    plotfigure.show = False
    

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.xlimits = [-400,400]
    plotaxes.ylimits = [-15,15]
    plotaxes.title = 'Zoom near left shore'

    plotaxes.afteraxes = fixticks

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.show = False
    plotitem.plot_var = geoplot.surface

    plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    plotitem.plot_var = geoplot.surface
    plotitem.plot_var2 = geoplot.topo
    plotitem.color = [.7,.7,1]
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.surface
    plotitem.color = 'b'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.topo
    plotitem.color = 'g'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1


    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.xlimits = [19600,20400]
    plotaxes.ylimits = [-15,15]
    plotaxes.title = 'Zoom near right shore'

    plotaxes.afteraxes = fixticks

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.show = False
    plotitem.plot_var = geoplot.surface

    plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    plotitem.plot_var = geoplot.surface
    plotitem.plot_var2 = geoplot.topo
    plotitem.color = [.7,.7,1]
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.surface
    plotitem.color = 'b'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.topo
    plotitem.color = 'g'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1



    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='q', figno=300, \
                                         type='each_gauge')
    plotfigure.clf_each_gauge = True
    plotfigure.kwargs = {'figsize': (9,5)}

    plotaxes = plotfigure.new_plotaxes()

    def addgrid(current_data):
        from pylab import title,grid,xlabel,ylabel
        if current_data.gaugeno == 1:
            title('Surface at river mouth')
        elif current_data.gaugeno == 2:
            title('Surface 2 km upstream')
        grid(True)
        xlabel('Time (seconds)')
        ylabel('Elevation (m)')

    plotaxes.afteraxes = addgrid
    plotaxes.xlimits = 'auto'
    #plotaxes.ylimits = [-1,15]
    plotaxes.ylimits = [-6,15]
    plotaxes.title = 'Water depth'
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 2  # eta
    plotitem.plotstyle = 'b-'


    plotdata.printfigs = True          # Whether to output figures
    plotdata.print_format = 'png'      # What type of output format
    plotdata.print_framenos = 'all'      # Which frames to output
    plotdata.print_fignos = [2]        # Which figures to print
    plotdata.html = True               # Whether to create HTML files
    plotdata.latex = False             # Whether to make LaTeX output
    plotdata.parallel = True

    return plotdata
