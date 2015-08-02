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

def plot(fname,ax,Q2,T=1,color='r',edgecolor='none',alpha=1):
  path='../CJ15data/FITS_dquark_series'
  DATA=FITPACK().get_PDFs(path+'/'+fname)
  D=DATA['Q2'][Q2]
  print D.keys()
  p1,=ax.plot(D['x'],D['d/u'],color+'-') 
  p2=fill_between(D['x'],
    (D['d/u']-D['err-d/u']*T),
    (D['d/u']+D['err-d/u']*T),
    ax=ax,
    facecolor=color,
    edgecolor=edgecolor,
    alpha=alpha)
  return (p2,p1)

def main():

  T=10
  Q2=100.0
  ax=py.subplot(111)

  f1='ddat_D_BNS.pdf'
  f2='ddat_D_BNS_Zrap_Lasy.pdf'
  f3='ddat_D_BNS_Zrap_Lasy_Wasy.pdf'
  f4='ddat_D_BNS_Zrap_Lasy_WCDF.pdf'
  f5='ddat_D_BNS_Zrap.pdf'
  f6='ddat_D.pdf'
  f7='ddat_noD.pdf'
  f8='ddat_nonuk_D_BNS.pdf'
  f9='ddat_nonuk_D_BNS_Zrap_Lasy.pdf'
  f10='ddat_nonuk_D_BNS_Zrap_Lasy_Wasy.pdf'
  f11='ddat_nonuk_D_BNS_Zrap_Lasy_WCDF.pdf'
  f12='ddat_nonuk_D_BNS_Zrap.pdf'
  f13='ddat_nonuk_D.pdf'


  pf6=plot(f6,ax,Q2,T,color='r',alpha=0.4)
  pf1=plot(f1,ax,Q2,T,color='b',alpha=0.4)
  pf2=plot(f2,ax,Q2,T,color='g',alpha=0.4)
  pf3=plot(f3,ax,Q2,T,color='g',alpha=0.4)

  H,L=[],[]
  H.append(pf6);L.append(tex('DIS'))
  H.append(pf1);L.append(tex('BONUS'))
  H.append(pf2);L.append(tex('Lasy+Zrap'))
  H.append(pf3);L.append(tex('Wasy+Zrap'))
  ax.legend(H,L,loc=1,frameon=0,fontsize=20)
  #ax.set_ylim(-0.2,1.2)
  #ax.set_xlim(0,0.4)
  py.savefig('gallery/d_over_u.pdf')

if __name__=='__main__':

  main()


