#!/usr/bin/env python
import sys,os
sys.path.insert(1,'../../')
import numpy as np
import pylab as py
import tools
from tools import tex,plot_band,fill_between
import lhapdf
import matplotlib.gridspec as gridspec
from  matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text',usetex=True)
from master import COMPOSER,FITPACK

# prepare grid 
gs = gridspec.GridSpec(1,2) 
gs.update(left=0.11,right=0.98,wspace=0,hspace=0.0,top=0.98,bottom=0.1)

# construct (L)eft panel
axL=py.subplot(gs[0,0])
axL.spines['right'].set_visible(False)
axL.yaxis.set_ticks_position('left')
axL.semilogx()
axL.set_xticks(10**np.linspace(-4,-2,3))
#axL.set_xlim(1e-5,0.1)
axL.axvline(0.1,color='k',ls='--',alpha=0.5)
axL.tick_params(axis='both', which='major', labelsize=20)

# construct (R)ight panel
axR=py.subplot(gs[0,1])
axR.spines['left'].set_visible(False)
axR.axes.yaxis.set_ticklabels([])
axR.axes.get_yaxis().set_ticks([])
axR.set_xticks(np.arange(0.1,1.1,0.2))
axR.set_xlim(0.1,1.0)
axR.tick_params(axis='both', which='major', labelsize=17)

CJ=COMPOSER(name='CJ15_NLO_KP_AV18')
def plot(flav,color,factor=1):
  XL=10**np.linspace(-4,-1,100)
  XR=np.linspace(0.1,1,100)
  L=CJ.get_xpdf(flav,X=XL,Q2=Q2)
  R=CJ.get_xpdf(flav,X=XR,Q2=Q2)
  axL.plot(XL,L['xf0']*factor,color=color,ls='-')
  p1,=axR.plot(XR,R['xf0']*factor,color=color,ls='-',label=flav)
  fill_between(XL,
    factor*(L['xf0']-L['dxf-']*10),
    factor*(L['xf0']+L['dxf+']*10),
    ax=axL,
    facecolor=color,
    edgecolor='none',
    alpha=0.5)
  p2=fill_between(XR,
    factor*(R['xf0']-R['dxf-']*10),
    factor*(R['xf0']+R['dxf+']*10),
    ax=axR,
    facecolor=color,
    edgecolor='none',
    alpha=0.5)#,hatch=None)
  return (p2,p1)

Q2=10
LH={}
LH['u']=plot('u','r')
LH['d']=plot('d','b')
LH['db+ub']=plot('db+ub','g')
LH['db-ub']=plot('db-ub','m')
LH['g']=plot('g','y',factor=1e-1)

# set limits
ylimL = axL.get_ylim()
ylimR = axR.get_ylim()
ymin=-0.1#np.amin([ylimL[0],ylimL[1],ylimR[0],ylimR[1]])
ymax=1.6#np.amax([ylimL[0],ylimL[1],ylimR[0],ylimR[1]])
axL.set_ylim(ymin,ymax)
axR.set_ylim(ymin,ymax)

# axis labels
axL.set_ylabel('$xf(x,Q^2)$',size=25)
axL.set_xlabel('$x$',size=25)
axL.xaxis.set_label_coords(1.0,-0.04,transform=axL.transAxes)

# legend
labelmap={}
labelmap['u']='$u$'
labelmap['d']='$d$'
labelmap['db+ub']=r'$\bar{d}+\bar{u}$'
labelmap['db-ub']=r'$\bar{d}-\bar{u}$'
labelmap['g']='$g/10$'

H,L=[],[]
for k in ['u','d','db+ub','db-ub','g']:
  p12=LH[k]
  H.append(LH[k])
  L.append(labelmap[k])
axR.legend(H,L,loc=1,frameon=0,fontsize=20)

# info
axL.text(0.1,0.1,'$Q^2=%0.0f$'%(Q2)+tex('~GeV^2'),transform=axL.transAxes,size=20)

# the end
py.savefig('gallery/xpdf.pdf')
