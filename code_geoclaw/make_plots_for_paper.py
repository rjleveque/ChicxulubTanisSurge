"""
Make plots of solution at selected times for the 4 cases used
in Figure 4 of the paper and Figure S3 of the SI.

"""
from pylab import *
import os

savefig_ext = '.pdf'
figdir = './figures'
os.system('mkdir -p %s' % figdir)

def save_figure(fname):
    """Save figure to figdir with desired extension"""
    full_fname = os.path.join(figdir,fname) + savefig_ext
    savefig(full_fname, bbox_inches='tight')
    print('Created %s' % full_fname)


def make_plot(jplot,frameno,outdir):

    subplot(4,2,jplot)
    fname = os.path.join(outdir, 'salt.txt')
    d = loadtxt(fname)
    tsalt = d[:,0]
    xsalt = d[:,1]
    print('Loaded xsalt from %s, min = %.1f, max = %.1f' \
            % (fname, xsalt.min(),xsalt.max()))
                    
    griddata = loadtxt(outdir + '/celledges.txt', skiprows=1)
    xgrid = griddata[:,0]
    zgrid = griddata[:,1]
    xcell = 0.5*(xgrid[:-1]+xgrid[1:])
    zcell = 0.5*(zgrid[:-1]+zgrid[1:])


    fname = outdir + '/fort.q%s' % str(frameno).zfill(4)
    q = loadtxt(fname, skiprows=6)
    eta = q[:,-1]
    fname = outdir + '/fort.t%s' % str(frameno).zfill(4)
    t = float(open(fname).readline().split()[0])
    print('t = %g' % t)
    jsalt = where(tsalt <= t)[0].max()
    xs = xsalt[jsalt]
    swater = where(xcell<xs, eta, nan)
    fwater = where(xcell>=xs, eta, nan)
    

    fill_between(xcell,zcell,swater,color=[1,.6,.6])
    fill_between(xcell,zcell,fwater,color=[.8,.8,1])
    plot(xcell,zcell,'g')
    plot(xcell, eta, 'b')
    xlim(-10.2e3,6e3)
    yticks(range(-20,21,5))
    xticks(rotation=20)
    grid(True)
    
    ysalt = (-10,13)
    plot([xs,xs],ysalt,'k--')
    text(xs, 12, '%.0f m' % xs, 
         horizontalalignment='center', fontsize=8,
          alpha=1,
         bbox={'facecolor':'yellow','edgecolor':'k','lw':0.5})
    text(-9000, 13, 't = %.0f sec' % t, fontsize=10, alpha=1,
         horizontalalignment='left', verticalalignment='top',
         bbox={'facecolor':'w','edgecolor':'k','lw':0.5})
    ylim(-15,15)


#framenos = [56,64,69,81]
framenos = [60,64,76,82]
       
# Cases 1 and 2:

fig,axs = subplots(4,2,sharex=True,sharey=True,figsize=(10,8),squeeze=True)

outdir = '_output_slope1_200m_r10'
subplot(4,2,1)
title('(a) Case 1: River depth 10 meters')
ylabel('meters')
for k,frameno in enumerate(framenos):
    jplot = 2*k+1
    make_plot(jplot,frameno,outdir)
xlabel('meters')

outdir = '_output_slope1_200m_r1'
subplot(4,2,2)
title('(b) Case 2: River depth 1 meter')
#ylabel('meters')
for k,frameno in enumerate(framenos):
    jplot = 2*k+2
    make_plot(jplot,frameno,outdir)
xlabel('meters')
    
tight_layout()
save_figure('numerical12')
close(fig)

if 1:
    # Cases 3 and 4 for SI:

    fig,axs = subplots(4,2,sharex=True,sharey=True,figsize=(10,8),squeeze=True)
    outdir = '_output_slope5000_200m_r10'
    subplot(4,2,1)
    title('(a) Case 3: River depth 10 meters')
    ylabel('meters')
    for k,frameno in enumerate(framenos):
        jplot = 2*k+1
        make_plot(jplot,frameno,outdir)
    xlabel('meters')

    outdir = '_output_slope5000_200m_r1'
    subplot(4,2,2)
    title('(b) Case 4: River depth 1 meter')
    #ylabel('meters')
    for k,frameno in enumerate(framenos):
        jplot = 2*k+2
        make_plot(jplot,frameno,outdir)
    xlabel('meters')
        
    tight_layout()
    save_figure('numerical34')
