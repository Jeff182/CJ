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
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('font',**{'family':'sans-serif','sans-serif':['Times-Roman']})
rc('text',usetex=True)
import pandas as pd
from  scipy.interpolate import  interp1d

class LASY(object):

  def __init__(self):
    self.load_data()
    self.make_plot()

  def load_data(self):
    F=open('../CJ15data/CJ15_NLO_KP_AV18/CJ15_Wasym.dat')
    L=F.readlines()
    F.close()
    L=[l.strip() for l in L]
    L=[l for l in L if l!='']
    H=L[0].split()
    L=[l.split() for l in L[1:]]
    DF=pd.DataFrame(L,columns=H)
    DF=DF.convert_objects(convert_numeric=True)

    D={}
    D['cdfLasy05'] = DF[DF.ITYPE=='cdfLasy05']
    D['d0Lasy_e15']= DF[DF.ITYPE=='d0Lasy_e15']
    D['d0Lasy13']  = DF[DF.ITYPE=='d0Lasy13']
    self.D=D

  def make_plot(self):
    D=self.D
    ax=py.subplot(111)

    Y=D['d0Lasy_e15']['Y']
    T=D['d0Lasy_e15']['THEORY']
    ET=D['d0Lasy_e15']['ERROR']

    iT= interp1d(Y,T,kind='cubic')
    iET= interp1d(Y,ET, kind='cubic')
    Y=np.linspace(np.amin(Y),np.amax(Y),100)

    T=10
    p2=fill_between(Y,iT(Y)-iET(Y)*T,iT(Y)+iET(Y)*T,
      ax=ax,
      facecolor='yellow',
      edgecolor='yellow')
    p1,=ax.plot(Y,iT(Y),'r-')
    H=[(p2,p1)]
    L=[tex('CJ15')]

    dmap={}
    dmap['cdfLasy05']  = {'color':'g','marker':'d'}
    dmap['d0Lasy_e15'] = {'color':'b','marker':'o'}
    dmap['d0Lasy13']   = {'color':'c','marker':'^'}

    for k in D.keys():
      color=dmap[k]['color']
      marker=dmap[k]['marker']
      markersize=4
      p3=ax.errorbar(D[k]['Y'],D[k]['DATA'],\
        yerr=D[k]['DERROR'],fmt=color+marker,mfc=color,mec=color,\
        markersize=markersize,zorder=1,alpha=0.9)
      H.append(p3)
      #L.append(tex(k.replace('_','')))
      L=[tex('CJ15'),tex('CDF')+r'$\ e$',tex('D\O')+r'$\ \mu$',tex('D\O')+r'$\ e$']

    ax.set_xlabel(r'$\eta_{\ell}$',size=25)
    ax.set_ylabel(r'$A_{\ell}$',size=25)

    ax.legend(H,L,frameon=0,loc=3,fontsize=22,numpoints=1)
   
    ###ax.text(0.5,0.8,tex('nrep=%d'%nrows),transform=ax.transAxes,size=20)
    py.tick_params(axis='both',labelsize=20)
    py.tight_layout()
    py.savefig('gallery/Lasy.pdf')
    py.close()

if __name__=='__main__':

  LASY()
