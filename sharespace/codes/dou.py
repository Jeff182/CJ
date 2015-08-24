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
rc('font',**{'family':'sans-serif','sans-serif':['Times-Roman']})
rc('text',usetex=True)
from master import COMPOSER,FITPACK,COMPOSER4NNPDF

def plot_dou(ax,com,Q2,X,color,ls='-',T=1,hatch=None,alpha=0.4,facecolor='none',edgecolor='none'):
  dou=com.get_dou(X=X,Q2=Q2)
  #ax.plot(X,dou['central'],color=color,ls=ls)
  #p1,=ax.plot(X,dou['central'],color=color,ls=ls)
  p2=fill_between(X,
    (dou['central']-dou['asym err -']*T),
    (dou['central']+dou['asym err +']*T),
    ax=ax,
    facecolor=facecolor,
    edgecolor=edgecolor,
    alpha=alpha,
    hatch=hatch)
  #return (p2,p1)
  return p2

CJ15=COMPOSER('CJ15_NLO_KP_AV18')
MMHT14=COMPOSER('MMHT2014nlo68cl')
JR14=COMPOSER('JR14NLO08VF')
CT14=COMPOSER('CT14nlo')
#NNPDF=COMPOSER4NNPDF('NNPDF30_nlo_as_0118')
#ABM11=COMPOSER('abm11_4n_nlo')
#HERA15=COMPOSER('HERAPDF15NLO_EIG')

ax=py.subplot(111)
Q2=10
X=np.linspace(1e-3,0.95,1000)

ct=plot_dou(ax,CT14,Q2,X,color='g',ls='-',T=1,hatch=None,alpha=0.7,facecolor='g',edgecolor='g')
mmht=plot_dou(ax,MMHT14,Q2,X,color='y',ls='-',T=1,hatch=None,alpha=0.8,facecolor='y',edgecolor='y')
jr=plot_dou(ax,JR14,Q2,X,color='b',ls='-',T=10,hatch='..',alpha=0.7,facecolor='none',edgecolor='b')
cj=plot_dou(ax,CJ15,Q2,X,color='r',ls='-',T=10,hatch=None,alpha=0.8,facecolor='r',edgecolor='r')
#hera=plot_dou(ax,HERA15,Q2,X,color='k',ls='-',T=1,hatch='||',alpha=1.0,facecolor='none',edgecolor='k')
#nnpdf=plot_dou(ax,NNPDF,Q2,X,color='k',ls='-',T=1,hatch='//',alpha=0.4,edgecolor='k')
#abm=plot_dou(ax,ABM11,Q2,X,color='b',ls='-',T=10,hatch='...',alpha=0.7,facecolor='none',edgecolor='b')


# legend
H=[cj,mmht,ct,jr]
L=[tex('CJ15'),tex('MMHT14'),tex('CT14'),tex('JR14')]
ax.legend(H,L,frameon=0,fontsize=20,loc=3,bbox_to_anchor=(0.5,0.57))
ax.set_ylabel(r'$d/u$',size=25)
ax.set_xlabel(r'$x$',size=25)
ax.set_xlim(0,0.95)
ax.set_ylim(0,1)
ax.tick_params(axis='both',labelsize=20)
py.tight_layout()
py.savefig('gallery/du0.pdf')
