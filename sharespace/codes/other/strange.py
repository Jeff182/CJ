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
from master import COMPOSER,FITPACK,COMPOSER4NNPDF
##################################################

#NNPDF=COMPOSER4NNPDF('NNPDF30_nlo_as_0118')
NNPDF=COMPOSER4NNPDF('NNPDF30_nlo_as_0118_nolhc_1000')

MMHT14=COMPOSER('MMHT2014nlo68cl')

ax=py.subplot(111)
Q2=1
#X=np.linspace(1e-3,1,100)
X=10**np.linspace(-3,0,100)
D=NNPDF.get_xpdf('s',X=X,Q2=Q2)
ax.plot(X,D['xf0'],'r-',label='NNPDF')
ax.plot(X,D['xf0']-D['dxf-'],'r:')
ax.plot(X,D['xf0']+D['dxf+'],'r:')

D=MMHT14.get_xpdf('s',X=X,Q2=Q2)
ax.plot(X,D['xf0'],'b-',label='MMHT')
ax.plot(X,D['xf0']-D['dxf-'],'b:')
ax.plot(X,D['xf0']+D['dxf+'],'b:')

ax.legend(frameon=0)
ax.semilogx()
ax.set_ylim(0,1)
py.savefig('gallery/strange.pdf')
