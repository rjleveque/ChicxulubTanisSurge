
"""
Create mp4 files for 4 animations.
"""

import matplotlib
matplotlib.use('Agg')

from clawpack.visclaw import animation_tools
from clawpack.visclaw import plotclaw
import os


def make_anim(outdir):
    plotdir = outdir.replace('output','plots')
    #if not os.isdir(plotdir):
    if 1:
        plotclaw.plotclaw(outdir,plotdir)
    title = plotdir.replace('_plots_','')
    anim = animation_tools.make_anim(plotdir=plotdir,figsize=(11,4))
    file_name = title + '.mp4'
    animation_tools.make_mp4(anim, file_name=file_name)
    
for h0r in [10,1]:
    for xwallend in [1,5e3]:
        outdir = '_output_slope%.0f_200m_r%s' % (xwallend, h0r)        
        make_anim(outdir)
