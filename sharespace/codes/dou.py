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

def plot_dou(ax,com,Q2,X,color,ls='-',T=1,hatch=None,alpha=0.4,facecolor='none',edgecolor='none'):

  dou=com.get_dou(X=X,Q2=Q2)

  ax.plot(X,dou['central'],color=color,ls=ls)
  p1,=ax.plot(X,dou['central'],color=color,ls=ls)
  p2=fill_between(X,
    (dou['central']-dou['asym err -']*T),
    (dou['central']+dou['asym err +']*T),
    ax=ax,
    facecolor=facecolor,
    edgecolor=edgecolor,
    alpha=alpha,
    hatch=hatch)

  return (p2,p1)

CJ=COMPOSER('CJ15_NLO_KP_AV18')
HERA15=COMPOSER('HERAPDF15NLO_EIG')
MMHT14=COMPOSER('MMHT2014nlo68cl')
NNPDF=COMPOSER4NNPDF('NNPDF30_nlo_as_0118')


ax=py.subplot(111)
Q2=10
X=np.linspace(1e-2,0.9,100)

cj=plot_dou(ax,CJ,Q2,X,'r',ls='-',T=10,hatch=None,alpha=1.0,facecolor='Yellow',edgecolor='none')
hera=plot_dou(ax,HERA15,Q2,X,'k',ls='-',T=1,hatch='||',alpha=1.0,facecolor='none',edgecolor='k')
mmht=plot_dou(ax,MMHT14,Q2,X,'b',ls='-',T=1,hatch='.',alpha=0.4,facecolor='none',edgecolor='b')
nnpdf=plot_dou(ax,NNPDF,Q2,X,'g',ls='-',T=1,hatch='//',alpha=0.4,facecolor='none',edgecolor='g')

# legend
H=[cj,hera,mmht,nnpdf]
L=[tex('CJ15'),tex('HERA15'),tex('MMHT14'),tex('NNPDF')]
ax.legend(H,L,frameon=0,fontsize=20,loc=3)#,bbox_to_anchor=(0.05, 1.0))
ax.set_ylabel(r'$d/u$',size=20)
ax.set_xlabel('$x$',size=20)
ax.set_ylim(-0.2,1.2)
ax.tick_params(axis='both',labelsize=20)
py.tight_layout()
py.savefig('gallery/dou.pdf')

