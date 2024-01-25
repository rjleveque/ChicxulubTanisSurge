#!/usr/bin/env python
# coding: utf-8



from pylab import *
from matplotlib.patches import Polygon



def pwcubic(xi, zl, zr, slopel, sloper, x):
    from numpy import where, zeros

    dx = xi[1:] - xi[:-1]

    s = (zl[1:] - zr[:-1]) / dx
    c2 = (s - sloper[:-1]) / dx
    c3 = (slopel[1:] - 2.*s + sloper[:-1]) / (dx**2)

    # set to linear function for x<xi[0] or x>= xi[-1]
    z = where(x < xi[0],  zl[0] + (x-xi[0]) * slopel[0], 
                          zr[-1] + (x-xi[-1]) * sloper[-1])

    # replace by appropriate cubic in intervals:
    for i in range(len(xi)-1):
        cubic = zr[i] + sloper[i]*(x - xi[i])                     + (x-xi[i])**2 * (c2[i] + c3[i]*(x-xi[i+1]))
        z = where((x >= xi[i]) & (x < xi[i+1]), cubic, z)

    return z



hl = 20.
hr = 5.
x = linspace(0, 60e3,1001)
xi = array([-1e3, 20e3, 40e3, 60e3])
Bl = array([-hl, -hl, -hr, -hr])
Br = Bl
slopel = array([0,0,0,0])
sloper = slopel
B = pwcubic(xi, Bl, Br, slopel, sloper, x)
xmouth = 45e3


color_salt = [1,.8,.8]
color_fresh = [.8,.8,1]
#fs = 9 # fontsize
fs = 14 # fontsize

def plot_fig(eta_salt,eta_fresh):
    fill_between(x,B,eta_salt,color=color_salt)
    fill_between(x,B,eta_fresh,color=color_fresh)
    plot(x,B,'g')
    plot(x,eta_salt,'b')
    plot(x,eta_fresh,'b')
    plot([xmouth,xmouth], [-13,0], 'b--')
    text(xmouth,-16,'River mouth',horizontalalignment='center', fontsize=fs)
    #ymax = max(eta_salt.max(),eta_fresh.max())
    #ylim(-22,ymax+4)
    ylim(-22,7)
    



def sine_hump(x,x0,L,eta_max):
    eta1 = eta_max*where(logical_or(x<x0, x>x0+L), 0, sin(pi*(x-x0)/L))
    return eta1

def add_arrow(x0,L,fraca,eta_max):
    xa = x0 + (1-fraca)*L/2
    dxa = fraca*L
    ya = eta_max + 1.5
    #arrow(xa,ya,dxa,0,head_width=0.5, head_length=300)
    arrow(xa,ya,dxa,0,head_width=1.3, head_length=800)


figure(figsize=(14,16))

subplot(511)

eta0 = zeros(x.shape)
eta_salt = ma.masked_where(x>xmouth, eta0)
eta_fresh = ma.masked_where(x<=xmouth, eta0)
#plot_fig(eta_salt,eta_fresh)

    
eta_max = 2.
L = 6e3

x0 = 2e3
eta_salt = ma.masked_where(x>xmouth, sine_hump(x,x0,L,eta_max))
eta_fresh = ma.masked_where(x<=xmouth, eta0)
plot_fig(eta_salt,eta_fresh)

x1 = x0+L
plot([x1,x1], [-20,0], 'k--')
          
add_arrow(x0,L,0.5,eta_max);

fss = 12

xlab = x1 + 5e3
arrow(xlab,-19.8,0,19.6,head_width=300, head_length=1, length_includes_head=True)
arrow(xlab,-0.2,0,-19.6,head_width=300, head_length=1, length_includes_head=True)
text(xlab+600,-10,'$H_1$',fontsize=fss,ha='left')

arrow(xlab,0.2,0,1.8,head_width=300, head_length=0.3, length_includes_head=True)
arrow(xlab,2.0,0,-1.8,head_width=300, head_length=0.3, length_includes_head=True)
text(xlab+600,1.0,'$\eta_1^*$',fontsize=fss,ha='left')

arrow(x0+400,-1,L-800,0,head_width=0.5, head_length=400, length_includes_head=True)
arrow(x1-400,-1,-(L-800),0,head_width=0.5, head_length=400, length_includes_head=True)
text(x0+L/2,-3.5,'$L_1$',fontsize=fss)

text(x0+L/2,5.,'wave speed',fontsize=fss,ha='center')

text(150,-18.5,'(a)',fontsize=fs)
axis(False)


#=======================

subplot(512)

x0 = 10e3
eta_salt = ma.masked_where(x>xmouth, sine_hump(x,x0,L,eta_max))
eta_fresh = ma.masked_where(x<=xmouth, eta0)
plot_fig(eta_salt,eta_fresh)

plot([x1,x1], [-20,0], 'k--')

A1 = 2/pi * L * eta_max
dx = A1/20.
dx = 3*dx # looks better
plot([x1+dx,x1+dx], [-20,0], 'k--')
fill_between([x1,x1+dx], [-20,-20], [0,0], facecolor=color_salt, hatch='///')
fill_between(x, eta0, eta_salt, facecolor=color_salt, hatch='///')
print(A1, dx)

add_arrow(x0,L,0.5,eta_max)

#arrow(x1,-20.5,dx,0,head_width=0.5, head_length=300, length_includes_head=True)
text(150,-18.5,'(b)',fontsize=fs)

annotate("equal areas $A_1$", 
            xy=(x1+dx,-15), xycoords='data',
            xytext=(x1+4*dx, -10), textcoords='data',
            arrowprops=dict(arrowstyle="->", connectionstyle="angle3"),
            fontsize=fss)

annotate("          ",
            xy=(x1+0.8*L,-0.55), xycoords='data',
            xytext=(x1+6.0*dx, -9), textcoords='data',
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3"))

axis(False)

#=======================

subplot(513)




x0 = 40e3
L2 = 3.5e3
eta_max2 = 4.
eta_salt = ma.masked_where(x>xmouth, sine_hump(x,x0,L2,eta_max2))
eta_fresh = ma.masked_where(x<=xmouth, eta0)
plot_fig(eta_salt,eta_fresh)
add_arrow(x0,L2,0.5,eta_max2)
text(150,-18.5,'(c)',fontsize=fs)

axis(False)

#=======================

subplot(514)

A2 = 2/pi*L2*eta_max2
dxsalt = A2 / 5.
xsalt_max = xmouth + dxsalt
xsalt = xmouth + dxsalt / 2.
x0 = xsalt - L2/2
eta_salt = ma.masked_where(x>xsalt, sine_hump(x,x0,L2,eta_max2))
eta_fresh = ma.masked_where(x<=xsalt, sine_hump(x,x0,L2,eta_max2))
plot_fig(eta_salt,eta_fresh)

plot([xsalt, xsalt], [-10,eta_max2], 'k--')
text(xsalt+800, -10, 'Salt/fresh interface', fontsize=fs);
add_arrow(x0,L2,0.5,eta_max2)
text(150,-18.5,'(d)',fontsize=fs)

axis(False)

#=======================

subplot(515)



A2 = 2/pi*L2*eta_max2
dxsalt = A2 / 5.
xsalt_max = xmouth + dxsalt
x0 = 50e3
eta_salt = ma.masked_where(x>xsalt_max, eta0)
eta_fresh = ma.masked_where(x<=xsalt_max, sine_hump(x,x0,L2,eta_max2))
plot_fig(eta_salt,eta_fresh)

plot([xsalt_max, xsalt_max], [-10,0], 'k--')
text(xsalt_max+800, -10, 'Salt/fresh interface',fontsize=fs);

fill_between([xmouth,xsalt_max], [-5,-5], [0,0], facecolor=color_salt, hatch='///')
fill_between(x, eta0, eta_fresh, facecolor=color_fresh, hatch='///')
add_arrow(x0,L2,0.5,eta_max2)
text(150,-18.5,'(e)',fontsize=fs)

xlab = 58e3
arrow(xlab,-4.8,0,4.6,head_width=300, head_length=1, length_includes_head=True)
arrow(xlab,-0.2,0,-4.6,head_width=300, head_length=1, length_includes_head=True)
text(xlab+600,-3,'$H_2$',fontsize=fss,ha='left')

arrow(xlab,0.2,0,eta_max2-0.4,head_width=300, head_length=0.3, length_includes_head=True)
arrow(xlab,eta_max2-0.2,0,-eta_max2+0.4,head_width=300, head_length=0.3, length_includes_head=True)
text(xlab+600,1.4,'$\eta_2^*$',fontsize=fss,ha='left')

annotate("equal areas $A_2$", 
            xy=(44500,-2), xycoords='data',
            xytext=(35000, 3), textcoords='data',
            arrowprops=dict(arrowstyle="->", connectionstyle="angle3"),
            fontsize=fss)
annotate("          ",
            xy=(50000,2), xycoords='data',
            xytext=(40000, 3), textcoords='data',
            arrowprops=dict(arrowstyle="->", connectionstyle="angle3"))
                
axis(False)

if 1:
    fname = 'cartoon.pdf'
    savefig(fname, bbox_inches='tight')
    print('Created ',fname)
