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
import pandas as pd
from master import COMPOSER,FITPACK

def plot(fname,ax,Q2,T=1,color='r',hatch='...',edgecolor='k',alpha=1):
  path='../CJ15data/FITS_dquark_series'
  DATA=FITPACK().get_PDFs(path+'/'+fname)
  D=DATA['Q2'][Q2]
  #print D.keys()
  #p1,=ax.plot(D['x'],D['d/u'],color+'-')
  p2=fill_between(D['x'],
    (D['d/u']-D['err-d/u']*T),
    (D['d/u']+D['err-d/u']*T),
    ax=ax,
    facecolor=color,
    edgecolor=edgecolor,
    alpha=alpha,
    hatch=hatch)
  #return (p2,p1)
  return p2

def main():

  T=10
  Q2=10.0
  ax=py.subplot(111)

  f1='ddat_D.pdf'
  f2='ddat_D_BNS.pdf'
  f3='ddat_D_BNS_Zrap.pdf'
  f4='ddat_D_BNS_Zrap_Lasy.pdf'
  f5='ddat_D_BNS_Zrap_Lasy_WCDF.pdf'
  f6='ddat_D_BNS_Zrap_Lasy_Wasy.pdf'
  g1='ddat_nonuk_D.pdf'
  g2='ddat_nonuk_D_BNS.pdf'
  g3='ddat_nonuk_D_BNS_Zrap.pdf'
  g4='ddat_nonuk_D_BNS_Zrap_Lasy.pdf'
  g5='ddat_nonuk_D_BNS_Zrap_Lasy_WCDF.pdf'
  g6='ddat_nonuk_D_BNS_Zrap_Lasy_Wasy.pdf'
  h1='ddat_noD.pdf'

# plot d/u fig. 1
  pf1=plot(f1,ax,Q2,T,color='y',hatch=None,edgecolor='y',alpha=1.0)
  pf2=plot(f2,ax,Q2,T,color='g',hatch=None,edgecolor='g',alpha=1.0)
  pf4=plot(f4,ax,Q2,T,color='b',hatch=None,edgecolor='b',alpha=1.0)
  pf6=plot(f6,ax,Q2,T,color='r',hatch=None,edgecolor='r',alpha=1.0)

# plot d/u fig. 2
  #pg6=plot(g6,ax,Q2,T,color='y',hatch='...',edgecolor='y',alpha=0.7)
  #ph1=plot(h1,ax,Q2,T,color='g',hatch=None,edgecolor='g',alpha=0.5)
  #pf6=plot(f6,ax,Q2,T,color='r',hatch=None,edgecolor='r',alpha=0.7)

# plot d/u fig. ?
  #pf1=plot(g1,ax,Q2,T,color='y',hatch=None,edgecolor='y',alpha=0.5)
  #pf2=plot(g2,ax,Q2,T,color='g',hatch=None,edgecolor='g',alpha=0.5)
  #pf3=plot(g3,ax,Q2,T,color='m',hatch=None,edgecolor='m',alpha=1.0)
  #pf4=plot(g4,ax,Q2,T,color='b',hatch=None,edgecolor='b',alpha=1.0)
  #pf5=plot(g5,ax,Q2,T,color='k',hatch=None,edgecolor='k',alpha=1.0)
  #pf6=plot(g6,ax,Q2,T,color='r',hatch=None,edgecolor='r',alpha=1.0)


  H,L=[],[]
# legend fig. 1
  H.append(pf1);L.append(tex('DIS\ only'))
  H.append(pf2);L.append(tex('+\ BONuS'))
  H.append(pf4);L.append(r'$+\ \ell\ {\rm asym}\ (\&\ Z\ {\rm rap})$')
  H.append(pf6);L.append(r'$+\ W\ {\rm asym}$')

# legend fig. 2
  #H.append(pg6);L.append(tex('no\ nuclear'))
  #H.append(ph1);L.append(tex('no\ deuteron'))
  #H.append(pf6);L.append(tex('CJ15'))

# legend fig. ?
  #H.append(pf1);L.append(tex('DIS\ only'))
  #H.append(pf2);L.append(tex('+\ BONuS'))
  #H.append(pf3);L.append(r'$+\ Z\ {\rm rap})$')
  #H.append(pf4);L.append(r'$+\ \ell\ {\rm asym}$')
  #H.append(pf5);L.append(r'$+\ W\ {\rm CDF}$')
  #H.append(pf6);L.append(r'$+\ W\ {\rm asym}$')

  ax.legend(H,L,loc=1,frameon=0,fontsize=20,bbox_to_anchor=(0.95,1))
  ax.set_ylim(0.0,1.0)
  ax.set_xlim(0.0,0.95)
  ax.set_xlabel('$x$',size=25)
  ax.set_ylabel('$d/u$',size=25)
  py.tick_params(axis='both',labelsize=20)
  py.savefig('gallery/du1.pdf')

if __name__=='__main__':

  main()
