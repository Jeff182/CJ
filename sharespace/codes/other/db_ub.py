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

def main():

  T=10
  Q2=100
  X=np.linspace(1e-3,0.4)

  CJ=COMPOSER(name='CJ15_NLO_KP_AV18')
  d=CJ.get_xpdf('db-ub',X=X,Q2=Q2)

  ax=py.subplot(111)
  p1,=ax.plot(X,d['xf0']/X,color='r',ls='-')
  p2=fill_between(X,
    (d['xf0']-d['dxf-']*T)/X,
    (d['xf0']+d['dxf+']*T)/X,
    ax=ax,
    facecolor='r',
    edgecolor='none',
    alpha=0.5)

  H=[(p2,p1)]
  L=[tex('CJ15')+'\ $(T=%d)$'%T]
  ax.legend(H,L,loc=1,frameon=0,fontsize=20)
  ax.set_ylim(-0.2,1.2)
  ax.axhline(0,ls='--',color='k',alpha=0.5)

  py.savefig('gallery/db_ub.pdf')

if __name__=='__main__':

  main()





